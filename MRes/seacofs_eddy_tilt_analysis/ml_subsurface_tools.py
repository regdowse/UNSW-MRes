"""Utilities for surface-to-subsurface eddy-tilt modelling.

The workflow predicts eastward and northward tilt components rather than
regressing directly on a circular bearing. All train/test and cross-validation
splits are grouped by eddy identity to prevent observations from the same eddy
appearing on both sides of a split.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, GroupShuffleSplit
from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


NUMERIC_FEATURES = ["lat", "w", "Omega", "Rc", "AR", "Age", "PV_gradient"]
CATEGORICAL_FEATURES = ["Cyc"]
FEATURES = NUMERIC_FEATURES + CATEGORICAL_FEATURES
TARGET_COMPONENTS = ["TiltEast", "TiltNorth"]


@dataclass
class ModelSplit:
    """Grouped train/test data for a single held-out-eddy experiment."""

    X_train: pd.DataFrame
    X_test: pd.DataFrame
    y_train: pd.DataFrame
    y_test: pd.DataFrame
    groups_train: pd.Series
    groups_test: pd.Series
    train_index: pd.Index
    test_index: pd.Index


def target_availability_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Summarise valid tilt coverage overall and by eddy polarity."""

    valid = df[["TiltDis", "TiltDir"]].notna().all(axis=1)
    rows = [
        {
            "group": "All",
            "rows": len(df),
            "valid_tilt_rows": int(valid.sum()),
            "valid_fraction": float(valid.mean()),
            "eddies": int(df["Eddy"].nunique()),
            "eddies_with_tilt": int(df.loc[valid, "Eddy"].nunique()),
        }
    ]
    for cyc, part in df.groupby("Cyc", dropna=False):
        part_valid = valid.loc[part.index]
        rows.append(
            {
                "group": str(cyc),
                "rows": len(part),
                "valid_tilt_rows": int(part_valid.sum()),
                "valid_fraction": float(part_valid.mean()),
                "eddies": int(part["Eddy"].nunique()),
                "eddies_with_tilt": int(part.loc[part_valid, "Eddy"].nunique()),
            }
        )
    return pd.DataFrame(rows).set_index("group")


def prepare_modelling_table(df: pd.DataFrame) -> pd.DataFrame:
    """Create a complete-case ML table and Cartesian tilt targets.

    ``TiltDir`` is assumed to be a compass bearing: 0 degrees is north and
    90 degrees is east, matching ``seacofs_tilt_tools.bearing_from_xy``.
    ``PV_gradient`` should be the magnitude produced from the core-mean PV
    calculation (``PV_grad_mag`` in ``seacofs_tilt_tools``).
    """

    required = {"Eddy", "TiltDis", "TiltDir", *FEATURES}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing required modelling columns: {missing}")

    out = df.copy()
    out["Cyc"] = out["Cyc"].astype("string")
    theta = np.deg2rad(out["TiltDir"].astype(float))
    out["TiltEast"] = out["TiltDis"] * np.sin(theta)
    out["TiltNorth"] = out["TiltDis"] * np.cos(theta)

    numeric_required = NUMERIC_FEATURES + ["TiltDis", "TiltDir", *TARGET_COMPONENTS]
    finite = np.isfinite(out[numeric_required].astype(float)).all(axis=1)
    complete = out[["Eddy", "Cyc"]].notna().all(axis=1)
    out = out.loc[finite & complete].copy()
    out["Eddy"] = out["Eddy"].astype(int)
    return out


def grouped_train_test_split(
    df: pd.DataFrame,
    *,
    test_size: float = 0.20,
    random_state: int = 42,
) -> ModelSplit:
    """Hold out complete eddies, never individual eddy-day rows."""

    splitter = GroupShuffleSplit(
        n_splits=1,
        test_size=test_size,
        random_state=random_state,
    )
    train_pos, test_pos = next(splitter.split(df, groups=df["Eddy"]))
    train = df.iloc[train_pos]
    test = df.iloc[test_pos]
    return ModelSplit(
        X_train=train[FEATURES].copy(),
        X_test=test[FEATURES].copy(),
        y_train=train[TARGET_COMPONENTS].copy(),
        y_test=test[TARGET_COMPONENTS].copy(),
        groups_train=train["Eddy"].copy(),
        groups_test=test["Eddy"].copy(),
        train_index=train.index,
        test_index=test.index,
    )


def _one_hot_encoder() -> OneHotEncoder:
    """Create a dense encoder across old and new scikit-learn versions."""

    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:  # scikit-learn < 1.2
        return OneHotEncoder(handle_unknown="ignore", sparse=False)


def make_preprocessor(*, scale_numeric: bool) -> ColumnTransformer:
    """Build leakage-safe preprocessing fitted only within each model fit."""

    numeric_steps = [("imputer", SimpleImputer(strategy="median"))]
    if scale_numeric:
        numeric_steps.append(("scaler", StandardScaler()))
    numeric = Pipeline(numeric_steps)
    categorical = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", _one_hot_encoder()),
        ]
    )
    return ColumnTransformer(
        [
            ("numeric", numeric, NUMERIC_FEATURES),
            ("categorical", categorical, CATEGORICAL_FEATURES),
        ],
        verbose_feature_names_out=False,
    )


def build_models(*, random_state: int = 42) -> dict[str, Pipeline]:
    """Return an interpretable linear model and a nonlinear tree model."""

    ridge = Pipeline(
        [
            ("preprocess", make_preprocessor(scale_numeric=True)),
            ("model", Ridge(alpha=10.0)),
        ]
    )
    gradient_boosting = Pipeline(
        [
            ("preprocess", make_preprocessor(scale_numeric=False)),
            (
                "model",
                MultiOutputRegressor(
                    HistGradientBoostingRegressor(
                        learning_rate=0.06,
                        max_iter=300,
                        max_leaf_nodes=31,
                        l2_regularization=1.0,
                        early_stopping=True,
                        random_state=random_state,
                    )
                ),
            ),
        ]
    )
    return {"Ridge": ridge, "Gradient boosting": gradient_boosting}


def fit_models(models: Mapping[str, Pipeline], split: ModelSplit) -> dict[str, Pipeline]:
    """Clone and fit each model on the grouped training partition."""

    fitted = {}
    for name, model in models.items():
        fitted[name] = clone(model).fit(split.X_train, split.y_train)
    return fitted


def components_to_polar(components) -> tuple[np.ndarray, np.ndarray]:
    """Convert east/north displacement components to distance and bearing."""

    values = np.asarray(components, dtype=float)
    east, north = values[:, 0], values[:, 1]
    distance = np.hypot(east, north)
    direction = (np.degrees(np.arctan2(east, north)) + 360.0) % 360.0
    return distance, direction


def angular_error_deg(observed, predicted) -> np.ndarray:
    """Smallest absolute difference between two compass bearings."""

    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    return np.abs((predicted - observed + 180.0) % 360.0 - 180.0)


def prediction_frame(
    observed_components,
    predicted_components,
    *,
    index=None,
) -> pd.DataFrame:
    """Combine component, distance, direction, and error diagnostics."""

    observed = np.asarray(observed_components, dtype=float)
    predicted = np.asarray(predicted_components, dtype=float)
    obs_dis, obs_dir = components_to_polar(observed)
    pred_dis, pred_dir = components_to_polar(predicted)
    return pd.DataFrame(
        {
            "observed_east": observed[:, 0],
            "observed_north": observed[:, 1],
            "predicted_east": predicted[:, 0],
            "predicted_north": predicted[:, 1],
            "observed_distance": obs_dis,
            "predicted_distance": pred_dis,
            "observed_direction": obs_dir,
            "predicted_direction": pred_dir,
            "vector_error": np.hypot(
                predicted[:, 0] - observed[:, 0],
                predicted[:, 1] - observed[:, 1],
            ),
            "angular_error": angular_error_deg(obs_dir, pred_dir),
        },
        index=index,
    )


def score_predictions(observed_components, predicted_components) -> dict[str, float]:
    """Return joint-vector, distance, component, and circular metrics."""

    pred = prediction_frame(observed_components, predicted_components)
    observed = np.asarray(observed_components, dtype=float)
    predicted = np.asarray(predicted_components, dtype=float)
    return {
        "vector_MAE_km": float(pred["vector_error"].mean()),
        "vector_RMSE_km": float(np.sqrt(np.mean(pred["vector_error"] ** 2))),
        "distance_MAE_km": float(
            mean_absolute_error(pred["observed_distance"], pred["predicted_distance"])
        ),
        "distance_RMSE_km": float(
            np.sqrt(
                mean_squared_error(
                    pred["observed_distance"], pred["predicted_distance"]
                )
            )
        ),
        "east_R2": float(r2_score(observed[:, 0], predicted[:, 0])),
        "north_R2": float(r2_score(observed[:, 1], predicted[:, 1])),
        "median_angular_error_deg": float(pred["angular_error"].median()),
        "within_30deg_fraction": float((pred["angular_error"] <= 30.0).mean()),
        "within_45deg_fraction": float((pred["angular_error"] <= 45.0).mean()),
    }


def mean_vector_baseline(y_train, n_predictions: int) -> np.ndarray:
    """Predict the training-set mean east/north tilt vector."""

    mean_vector = np.asarray(y_train, dtype=float).mean(axis=0)
    return np.tile(mean_vector, (n_predictions, 1))


def evaluate_models(
    fitted_models: Mapping[str, Pipeline],
    split: ModelSplit,
    *,
    include_baseline: bool = True,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Evaluate every model on eddies excluded from model fitting."""

    predictions = {}
    rows = []
    if include_baseline:
        baseline = mean_vector_baseline(split.y_train, len(split.y_test))
        predictions["Mean-vector baseline"] = prediction_frame(
            split.y_test, baseline, index=split.test_index
        )
        rows.append({"model": "Mean-vector baseline", **score_predictions(split.y_test, baseline)})

    for name, model in fitted_models.items():
        predicted = model.predict(split.X_test)
        predictions[name] = prediction_frame(split.y_test, predicted, index=split.test_index)
        rows.append({"model": name, **score_predictions(split.y_test, predicted)})
    scores = pd.DataFrame(rows).set_index("model").sort_values("vector_MAE_km")
    return scores, predictions


def grouped_cross_validation(
    models: Mapping[str, Pipeline],
    df: pd.DataFrame,
    *,
    n_splits: int = 5,
) -> pd.DataFrame:
    """Score models across folds containing mutually exclusive eddies."""

    groups = df["Eddy"]
    X = df[FEATURES]
    y = df[TARGET_COMPONENTS]
    splitter = GroupKFold(n_splits=n_splits)
    rows = []
    for fold, (train_pos, valid_pos) in enumerate(splitter.split(X, y, groups), start=1):
        for name, model in models.items():
            fitted = clone(model).fit(X.iloc[train_pos], y.iloc[train_pos])
            predicted = fitted.predict(X.iloc[valid_pos])
            rows.append(
                {
                    "model": name,
                    "fold": fold,
                    "validation_eddies": int(groups.iloc[valid_pos].nunique()),
                    **score_predictions(y.iloc[valid_pos], predicted),
                }
            )
    return pd.DataFrame(rows)


def performance_by_group(
    prediction: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    group_column: str = "Cyc",
) -> pd.DataFrame:
    """Report performance for groups such as AE and CE on held-out rows."""

    joined = prediction.join(metadata[[group_column]], how="left")
    rows = []
    for group, part in joined.groupby(group_column, dropna=False):
        rows.append(
            {
                group_column: group,
                "rows": len(part),
                "vector_MAE_km": part["vector_error"].mean(),
                "distance_MAE_km": np.abs(
                    part["predicted_distance"] - part["observed_distance"]
                ).mean(),
                "median_angular_error_deg": part["angular_error"].median(),
                "within_30deg_fraction": (part["angular_error"] <= 30).mean(),
            }
        )
    return pd.DataFrame(rows).set_index(group_column)


def direction_performance_by_tilt(
    prediction: pd.DataFrame,
    thresholds=(0.0, 5.0, 10.0, 20.0),
) -> pd.DataFrame:
    """Show angular skill after excluding near-zero, ill-defined bearings."""

    rows = []
    for threshold in thresholds:
        part = prediction[prediction["observed_distance"] >= threshold]
        rows.append(
            {
                "minimum_tilt_km": threshold,
                "rows": len(part),
                "median_angular_error_deg": part["angular_error"].median(),
                "within_30deg_fraction": (part["angular_error"] <= 30).mean(),
                "within_45deg_fraction": (part["angular_error"] <= 45).mean(),
            }
        )
    return pd.DataFrame(rows).set_index("minimum_tilt_km")


def raw_feature_permutation_importance(
    model: Pipeline,
    X_test: pd.DataFrame,
    y_test: pd.DataFrame,
    *,
    n_repeats: int = 10,
    random_state: int = 42,
) -> pd.DataFrame:
    """Permutation importance using negative mean displacement-vector error."""

    def vector_skill(estimator, X, y):
        predicted = estimator.predict(X)
        return -score_predictions(y, predicted)["vector_MAE_km"]

    result = permutation_importance(
        model,
        X_test,
        y_test,
        scoring=vector_skill,
        n_repeats=n_repeats,
        random_state=random_state,
        n_jobs=-1,
    )
    return pd.DataFrame(
        {
            "feature": X_test.columns,
            "importance_mean": result.importances_mean,
            "importance_std": result.importances_std,
        }
    ).sort_values("importance_mean", ascending=False)


def plot_model_comparison(scores: pd.DataFrame):
    """Plot primary held-out vector and angular errors."""

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ordered = scores.sort_values("vector_MAE_km")
    axes[0].barh(ordered.index, ordered["vector_MAE_km"], color="steelblue")
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Mean tilt-vector error (km)")
    axes[0].set_title("Held-out eddies")
    axes[1].barh(ordered.index, ordered["median_angular_error_deg"], color="darkorange")
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Median angular error (degrees)")
    axes[1].set_title("Direction")
    fig.tight_layout()
    return fig, axes


def plot_prediction_diagnostics(prediction: pd.DataFrame, *, title: str = "Model"):
    """Compare observed and predicted distances, components, and errors."""

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    axes[0].hexbin(
        prediction["observed_distance"],
        prediction["predicted_distance"],
        gridsize=45,
        mincnt=1,
        cmap="viridis",
    )
    limit = np.nanpercentile(
        np.r_[prediction["observed_distance"], prediction["predicted_distance"]], 99
    )
    axes[0].plot([0, limit], [0, limit], "k--", lw=1)
    axes[0].set(xlim=(0, limit), ylim=(0, limit), xlabel="Observed distance (km)", ylabel="Predicted distance (km)")

    axes[1].hexbin(
        prediction["observed_east"],
        prediction["observed_north"],
        C=prediction["vector_error"],
        reduce_C_function=np.mean,
        gridsize=40,
        mincnt=1,
        cmap="magma",
    )
    axes[1].set_aspect("equal", adjustable="box")
    axes[1].set(xlabel="Observed eastward tilt (km)", ylabel="Observed northward tilt (km)")
    axes[1].set_title("Mean vector error by observed tilt")

    axes[2].hist(prediction["angular_error"], bins=np.arange(0, 185, 5), color="slateblue", alpha=0.85)
    axes[2].axvline(90, color="k", ls="--", lw=1, label="random-direction median")
    axes[2].set(xlabel="Angular error (degrees)", ylabel="Count", xlim=(0, 180))
    axes[2].legend(frameon=False)
    fig.suptitle(title)
    fig.tight_layout()
    return fig, axes


def plot_permutation_importance(importance: pd.DataFrame):
    """Plot held-out raw-feature permutation importance."""

    ordered = importance.sort_values("importance_mean")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.barh(
        ordered["feature"],
        ordered["importance_mean"],
        xerr=ordered["importance_std"],
        color="teal",
        alpha=0.8,
    )
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Increase in vector MAE after permutation (km)")
    ax.set_title("Permutation importance on held-out eddies")
    fig.tight_layout()
    return fig, ax
