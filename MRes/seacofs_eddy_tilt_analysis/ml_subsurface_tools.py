"""Leakage-safe tools for predicting subsurface eddy tilt from surface data."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, GroupShuffleSplit
from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


STATIC_FEATURES = ["beta", "h"]
DYNAMIC_FEATURES = [
    "prop_east_km_day", "prop_north_km_day",
    "ellipse_major_cos2", "ellipse_major_sin2", "Omega", #"Rc", #"norm_time"
]
PV_MAG_FEATURES = ["PV_grad_mag"]
PV_COMPONENT_FEATURES = ["PV_grad_east", "PV_grad_north"]
FEATURES = STATIC_FEATURES + DYNAMIC_FEATURES + PV_MAG_FEATURES + PV_COMPONENT_FEATURES
TARGET_COMPONENTS = ["TiltEast", "TiltNorth"]
RELATIVE_TARGETS = ["LogTiltDis", "DeltaTiltEast", "DeltaTiltNorth"]
TARGET_MODES = ("cartesian", "pv_relative")

FEATURE_SETS = {
    "surface_without_PV": STATIC_FEATURES + DYNAMIC_FEATURES,
    "PV_magnitude_only": STATIC_FEATURES + DYNAMIC_FEATURES + PV_MAG_FEATURES,
    "PV_components_only": STATIC_FEATURES + DYNAMIC_FEATURES + PV_COMPONENT_FEATURES,
    "full_PV_vector": FEATURES,
    "full_without_dynamics": STATIC_FEATURES + PV_MAG_FEATURES + PV_COMPONENT_FEATURES,
}


@dataclass
class ModelSplit:
    """Training and untouched test observations split by complete eddies."""

    train: pd.DataFrame
    test: pd.DataFrame

    @property
    def groups_train(self) -> pd.Series:
        return self.train["Eddy"]

    @property
    def groups_test(self) -> pd.Series:
        return self.test["Eddy"]


def target_availability_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Summarise valid tilt coverage overall and by polarity."""

    valid = df[["TiltDis", "TiltDir"]].notna().all(axis=1)
    rows = []
    for label, part in [("All", df), *list(df.groupby("Cyc", dropna=False))]:
        part_valid = valid.loc[part.index]
        rows.append({
            "group": str(label), "rows": len(part),
            "valid_tilt_rows": int(part_valid.sum()),
            "valid_fraction": float(part_valid.mean()),
            "eddies": int(part["Eddy"].nunique()),
            "eddies_with_tilt": int(part.loc[part_valid, "Eddy"].nunique()),
        })
    return pd.DataFrame(rows).set_index("group")


def _major_axis_encoding(q11, q12, q22) -> tuple[np.ndarray, np.ndarray]:
    """Return sin(2 theta), cos(2 theta) for the ellipse's major axis."""

    q11, q12, q22 = map(lambda x: np.asarray(x, dtype=float), (q11, q12, q22))
    scale = np.hypot(q11 - q22, 2.0 * q12)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos2 = -(q11 - q22) / scale
        sin2 = -(2.0 * q12) / scale
    cos2[scale == 0] = np.nan
    sin2[scale == 0] = np.nan
    return sin2, cos2


def prepare_modelling_table(df: pd.DataFrame) -> pd.DataFrame:
    """Engineer predictors/targets; retain missing predictors for train-only imputation."""

    required = {
        "Eddy", "Day", "Cyc", "TiltDis", "TiltDir", "beta", "h",
        "Omega", "Rc", "xc", "yc", "q11", "q12", "q22",
        "PV_grad_mag", "PV_grad_theta",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing required modelling columns: {missing}")

    out = df.sort_values(["Eddy", "Day"]).copy()
    out["Cyc"] = out["Cyc"].astype("string")
    out["Day_idx"] = out.groupby("Eddy").cumcount()
    last_idx = out.groupby("Eddy")["Day_idx"].transform("max")
    out["norm_time"] = np.where(last_idx > 0, out["Day_idx"] / last_idx, 0.0)

    dt = out.groupby("Eddy")["Day"].diff().astype(float).where(lambda x: x > 0)
    out["prop_east_km_day"] = out.groupby("Eddy")["xc"].diff() / dt
    out["prop_north_km_day"] = out.groupby("Eddy")["yc"].diff() / dt
    # sin2, cos2 = _major_axis_encoding(out["q11"], out["q12"], out["q22"])
    # out["ellipse_major_sin2"] = sin2
    # out["ellipse_major_cos2"] = cos2
    # q11, q12 and q22 describe orientation relative to a model grid
    # rotated 20° clockwise from geographic north.
    grid_rotation_rad = np.deg2rad(20.0)
    sin2_grid, cos2_grid = _major_axis_encoding(
        out["q11"], out["q12"], out["q22"]
    )
    # Rotate the doubled-angle representation:
    # theta_geographic = theta_grid + 20°
    rotation_2theta = 2.0 * grid_rotation_rad
    out["ellipse_major_sin2"] = (
        sin2_grid * np.cos(rotation_2theta)
        + cos2_grid * np.sin(rotation_2theta)
    )
    out["ellipse_major_cos2"] = (
        cos2_grid * np.cos(rotation_2theta)
        - sin2_grid * np.sin(rotation_2theta)
    )

    pv_theta = np.deg2rad(out["PV_grad_theta"].astype(float))
    out["PV_grad_east"] = out["PV_grad_mag"] * np.sin(pv_theta)
    out["PV_grad_north"] = out["PV_grad_mag"] * np.cos(pv_theta)
    # Polarity-specific hypothesis: CE tilt is aligned with the PV gradient,
    # whereas AE tilt is aligned with the direction opposite the PV gradient.
    out["PV_reference_theta"] = np.where(
        out["Cyc"].eq("AE"),
        (out["PV_grad_theta"] + 180.0) % 360.0,
        out["PV_grad_theta"],
    )
    tilt_theta = np.deg2rad(out["TiltDir"].astype(float))
    out["TiltEast"] = out["TiltDis"] * np.sin(tilt_theta)
    out["TiltNorth"] = out["TiltDis"] * np.cos(tilt_theta)
    out["LogTiltDis"] = np.log1p(out["TiltDis"].clip(lower=0))
    delta = np.deg2rad((out["TiltDir"] - out["PV_reference_theta"] + 180.0) % 360.0 - 180.0)
    out["DeltaTiltEast"] = np.sin(delta)
    out["DeltaTiltNorth"] = np.cos(delta)

    targets = out[["TiltDis", "TiltDir", *TARGET_COMPONENTS, *RELATIVE_TARGETS]].astype(float)
    keep = np.isfinite(targets).all(axis=1) & out[["Eddy", "Cyc"]].notna().all(axis=1)
    out = out.loc[keep].copy()
    out["Eddy"] = out["Eddy"].astype(int)
    return out


def grouped_train_test_split(df, *, test_size=0.20, random_state=42) -> ModelSplit:
    """Reserve complete eddies for one final test evaluation."""

    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
    train_pos, test_pos = next(splitter.split(df, groups=df["Eddy"]))
    return ModelSplit(df.iloc[train_pos].copy(), df.iloc[test_pos].copy())


def eddy_equal_weights(groups: Sequence) -> np.ndarray:
    """Give every eddy approximately equal total weight."""

    groups = pd.Series(groups)
    weights = 1.0 / groups.map(groups.value_counts()).to_numpy(dtype=float)
    return weights / weights.mean()


def make_preprocessor(features, *, scale_numeric):
    steps = [("imputer", SimpleImputer(strategy="median", add_indicator=True))]
    if scale_numeric:
        steps.append(("scaler", StandardScaler()))
    return ColumnTransformer([("numeric", Pipeline(steps), list(features))], verbose_feature_names_out=False)


def build_model(family, features, params=None, *, random_state=42):
    """Construct one Ridge or histogram-gradient-boosting pipeline."""

    params = dict(params or {})
    if family == "Ridge":
        return Pipeline([
            ("preprocess", make_preprocessor(features, scale_numeric=True)),
            ("model", Ridge(**{"alpha": 10.0, **params})),
        ])
    if family == "Gradient boosting":
        defaults = {
            "learning_rate": 0.05, "max_iter": 400, "max_leaf_nodes": 31,
            "min_samples_leaf": 50, "l2_regularization": 1.0,
            "early_stopping": True, "random_state": random_state,
        }
        return Pipeline([
            ("preprocess", make_preprocessor(features, scale_numeric=False)),
            ("model", MultiOutputRegressor(HistGradientBoostingRegressor(**{**defaults, **params}))),
        ])
    raise ValueError(f"Unknown model family: {family}")


def default_candidates():
    """Compact hyperparameter and target-representation search space."""

    candidates = [
        {"family": "Ridge", "target_mode": mode, "params": {"alpha": alpha}}
        for mode, alpha in product(TARGET_MODES, [0.1, 1.0, 10.0, 100.0])
    ]
    boosting_grid = [
        {"learning_rate": 0.03, "max_leaf_nodes": 15, "min_samples_leaf": 50, "l2_regularization": 1.0},
        {"learning_rate": 0.05, "max_leaf_nodes": 31, "min_samples_leaf": 50, "l2_regularization": 1.0},
        {"learning_rate": 0.05, "max_leaf_nodes": 15, "min_samples_leaf": 100, "l2_regularization": 10.0},
        {"learning_rate": 0.10, "max_leaf_nodes": 15, "min_samples_leaf": 100, "l2_regularization": 1.0},
    ]
    candidates += [
        {"family": "Gradient boosting", "target_mode": mode, "params": params}
        for mode, params in product(TARGET_MODES, boosting_grid)
    ]
    return candidates


def target_columns(mode):
    if mode == "cartesian":
        return TARGET_COMPONENTS
    if mode == "pv_relative":
        return RELATIVE_TARGETS
    raise ValueError(f"Unknown target mode: {mode}")


def prediction_to_components(prediction, metadata, mode):
    predicted = np.asarray(prediction, dtype=float)
    if mode == "cartesian":
        return predicted
    distance = np.expm1(np.clip(predicted[:, 0], 0.0, None))
    delta = np.degrees(np.arctan2(predicted[:, 1], predicted[:, 2]))
    direction = (metadata["PV_reference_theta"].to_numpy(dtype=float) + delta) % 360.0
    theta = np.deg2rad(direction)
    return np.column_stack([distance * np.sin(theta), distance * np.cos(theta)])


def fit_pipeline(model, data, features, mode, *, equal_eddy_weight=True):
    kwargs = {"model__sample_weight": eddy_equal_weights(data["Eddy"])} if equal_eddy_weight else {}
    return model.fit(data[list(features)], data[target_columns(mode)], **kwargs)


def components_to_polar(components):
    values = np.asarray(components, dtype=float)
    east, north = values[:, 0], values[:, 1]
    return np.hypot(east, north), (np.degrees(np.arctan2(east, north)) + 360.0) % 360.0


def angular_error_deg(observed, predicted):
    return np.abs((np.asarray(predicted) - np.asarray(observed) + 180.0) % 360.0 - 180.0)


def prediction_frame(observed_components, predicted_components, *, index=None):
    observed, predicted = map(lambda x: np.asarray(x, dtype=float), (observed_components, predicted_components))
    obs_dis, obs_dir = components_to_polar(observed)
    pred_dis, pred_dir = components_to_polar(predicted)
    return pd.DataFrame({
        "observed_east": observed[:, 0], "observed_north": observed[:, 1],
        "predicted_east": predicted[:, 0], "predicted_north": predicted[:, 1],
        "observed_distance": obs_dis, "predicted_distance": pred_dis,
        "observed_direction": obs_dir, "predicted_direction": pred_dir,
        "vector_error": np.hypot(predicted[:, 0] - observed[:, 0], predicted[:, 1] - observed[:, 1]),
        "angular_error": angular_error_deg(obs_dir, pred_dir),
    }, index=index)


def score_predictions(observed_components, predicted_components, *, groups=None):
    pred = prediction_frame(observed_components, predicted_components)
    observed, predicted = map(lambda x: np.asarray(x, dtype=float), (observed_components, predicted_components))
    scores = {
        "vector_MAE_km": float(pred["vector_error"].mean()),
        "vector_RMSE_km": float(np.sqrt(np.mean(pred["vector_error"] ** 2))),
        "distance_MAE_km": float(np.abs(pred["predicted_distance"] - pred["observed_distance"]).mean()),
        "distance_RMSE_km": float(np.sqrt(mean_squared_error(pred["observed_distance"], pred["predicted_distance"]))),
        "east_R2": float(r2_score(observed[:, 0], predicted[:, 0])),
        "north_R2": float(r2_score(observed[:, 1], predicted[:, 1])),
        "median_angular_error_deg": float(pred["angular_error"].median()),
        "within_30deg_fraction": float((pred["angular_error"] <= 30).mean()),
    }
    if groups is not None:
        weights = eddy_equal_weights(groups)
        scores["eddy_weighted_vector_MAE_km"] = float(np.average(pred["vector_error"], weights=weights))
        scores["eddy_weighted_distance_MAE_km"] = float(np.average(np.abs(pred["predicted_distance"] - pred["observed_distance"]), weights=weights))
    return scores


def baseline_predictions(train, validation):
    # Compute climatology from per-eddy means so long-lived eddies do not
    # dominate the baseline against which equally weighted models are judged.
    mean_vector = train.groupby("Eddy")[TARGET_COMPONENTS].mean().mean().to_numpy(dtype=float)
    climatology = np.tile(mean_vector, (len(validation), 1))
    distance = float(train.groupby("Eddy")["TiltDis"].median().median())
    theta = np.deg2rad(validation["PV_reference_theta"].to_numpy(dtype=float))
    pv_direction = np.column_stack([distance * np.sin(theta), distance * np.cos(theta)])
    return {"Mean-vector baseline": climatology, "Polarity-aligned PV baseline": pv_direction}


def evaluate_candidate(candidate, train, validation, features=FEATURES, *, random_state=42):
    model = build_model(candidate["family"], features, candidate["params"], random_state=random_state)
    fitted = fit_pipeline(model, train, features, candidate["target_mode"])
    raw = fitted.predict(validation[list(features)])
    predicted = prediction_to_components(raw, validation, candidate["target_mode"])
    observed = validation[TARGET_COMPONENTS].to_numpy(dtype=float)
    return (
        score_predictions(observed, predicted, groups=validation["Eddy"]),
        fitted,
        prediction_frame(observed, predicted, index=validation.index),
    )


def tune_candidates_grouped(train, *, features=FEATURES, candidates=None, n_splits=5, random_state=42):
    """Select model family, targets, and hyperparameters using training eddies only."""

    candidates = list(candidates or default_candidates())
    splitter = GroupKFold(n_splits=n_splits)
    folds = list(splitter.split(train, groups=train["Eddy"]))
    rows = []
    for candidate_id, candidate in enumerate(candidates):
        for fold, (fit_pos, val_pos) in enumerate(folds, start=1):
            scores, _, _ = evaluate_candidate(candidate, train.iloc[fit_pos], train.iloc[val_pos], features, random_state=random_state + fold)
            rows.append({"candidate_id": candidate_id, "family": candidate["family"], "target_mode": candidate["target_mode"], "params": repr(candidate["params"]), "fold": fold, **scores})
    for fold, (fit_pos, val_pos) in enumerate(folds, start=1):
        fit_data, val_data = train.iloc[fit_pos], train.iloc[val_pos]
        observed = val_data[TARGET_COMPONENTS].to_numpy(dtype=float)
        for name, predicted in baseline_predictions(fit_data, val_data).items():
            rows.append({"candidate_id": name, "family": name, "target_mode": "baseline", "params": "{}", "fold": fold, **score_predictions(observed, predicted, groups=val_data["Eddy"])})
    return pd.DataFrame(rows)


def summarise_tuning(tuning):
    return tuning.groupby(["candidate_id", "family", "target_mode", "params"], dropna=False).agg(
        vector_MAE_mean=("vector_MAE_km", "mean"), vector_MAE_std=("vector_MAE_km", "std"),
        eddy_vector_MAE_mean=("eddy_weighted_vector_MAE_km", "mean"),
        distance_MAE_mean=("distance_MAE_km", "mean"),
        angular_error_median=("median_angular_error_deg", "median"),
        within_30deg_mean=("within_30deg_fraction", "mean"),
    ).sort_values("eddy_vector_MAE_mean")


def best_candidate(tuning, candidates=None):
    candidates = list(candidates or default_candidates())
    ids = pd.to_numeric(tuning["candidate_id"], errors="coerce")
    scores = tuning.loc[ids.notna()].assign(numeric_id=ids[ids.notna()].astype(int)).groupby("numeric_id")["eddy_weighted_vector_MAE_km"].mean()
    return dict(candidates[int(scores.idxmin())])


def fit_selected_and_test(candidate, split, *, features=FEATURES, random_state=42):
    scores, fitted, prediction = evaluate_candidate(candidate, split.train, split.test, features, random_state=random_state)
    rows = [{"model": "Selected model", **scores}]
    predictions = {"Selected model": prediction}
    observed = split.test[TARGET_COMPONENTS].to_numpy(dtype=float)
    for name, predicted in baseline_predictions(split.train, split.test).items():
        rows.append({"model": name, **score_predictions(observed, predicted, groups=split.test["Eddy"])})
        predictions[name] = prediction_frame(observed, predicted, index=split.test.index)
    scores_df = pd.DataFrame(rows).set_index("model").sort_values("eddy_weighted_vector_MAE_km")
    return fitted, prediction, scores_df, predictions


def feature_ablation_cv(train, candidate, *, feature_sets=FEATURE_SETS, n_splits=5, random_state=42):
    splitter = GroupKFold(n_splits=n_splits)
    folds = list(splitter.split(train, groups=train["Eddy"]))
    rows = []
    for name, features in feature_sets.items():
        for fold, (fit_pos, val_pos) in enumerate(folds, start=1):
            scores, _, _ = evaluate_candidate(candidate, train.iloc[fit_pos], train.iloc[val_pos], features, random_state=random_state + fold)
            rows.append({"feature_set": name, "fold": fold, **scores})
    return pd.DataFrame(rows)


def repeated_grouped_validation(data, candidate, *, features=FEATURES, n_repeats=5, test_size=0.20, random_state=100):
    rows = []
    for repeat in range(n_repeats):
        split = grouped_train_test_split(data, test_size=test_size, random_state=random_state + repeat)
        scores, _, _ = evaluate_candidate(candidate, split.train, split.test, features, random_state=random_state + repeat)
        rows.append({"repeat": repeat + 1, "model": "Selected model", **scores})
        observed = split.test[TARGET_COMPONENTS].to_numpy(dtype=float)
        for name, predicted in baseline_predictions(split.train, split.test).items():
            rows.append({"repeat": repeat + 1, "model": name, **score_predictions(observed, predicted, groups=split.test["Eddy"])})
    return pd.DataFrame(rows)


def direction_performance_by_tilt(prediction, thresholds=(0.0, 5.0, 10.0, 20.0)):
    rows = []
    for threshold in thresholds:
        part = prediction[prediction["observed_distance"] >= threshold]
        rows.append({"minimum_tilt_km": threshold, "rows": len(part), "median_angular_error_deg": part["angular_error"].median(), "within_30deg_fraction": (part["angular_error"] <= 30).mean()})
    return pd.DataFrame(rows).set_index("minimum_tilt_km")


def raw_feature_permutation_importance(fitted, test, candidate, *, features=FEATURES, n_repeats=10, random_state=42):
    y = test[target_columns(candidate["target_mode"])]
    def skill(estimator, X, raw_y):
        metadata = test.loc[X.index]
        predicted = prediction_to_components(estimator.predict(X), metadata, candidate["target_mode"])
        observed = metadata[TARGET_COMPONENTS].to_numpy(dtype=float)
        return -score_predictions(observed, predicted, groups=metadata["Eddy"])["eddy_weighted_vector_MAE_km"]
    result = permutation_importance(fitted, test[list(features)], y, scoring=skill, n_repeats=n_repeats, random_state=random_state, n_jobs=-1)
    return pd.DataFrame({"feature": features, "importance_mean": result.importances_mean, "importance_std": result.importances_std}).sort_values("importance_mean", ascending=False)


def plot_model_comparison(scores):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ordered = scores.sort_values("eddy_weighted_vector_MAE_km")
    axes[0].barh(ordered.index, ordered["eddy_weighted_vector_MAE_km"], color="steelblue")
    axes[0].invert_yaxis(); axes[0].set_xlabel("Equal-eddy-weighted vector MAE (km)")
    axes[1].barh(ordered.index, ordered["median_angular_error_deg"], color="darkorange")
    axes[1].invert_yaxis(); axes[1].set_xlabel("Median angular error (degrees)")
    fig.tight_layout(); return fig, axes


def plot_prediction_diagnostics(prediction, *, title):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    axes[0].hexbin(prediction["observed_distance"], prediction["predicted_distance"], gridsize=45, mincnt=1, cmap="viridis")
    limit = np.nanpercentile(np.r_[prediction["observed_distance"], prediction["predicted_distance"]], 99)
    axes[0].plot([0, limit], [0, limit], "k--", lw=1); axes[0].set(xlim=(0, limit), ylim=(0, limit), xlabel="Observed distance (km)", ylabel="Predicted distance (km)")
    axes[1].hexbin(prediction["observed_east"], prediction["observed_north"], C=prediction["vector_error"], reduce_C_function=np.mean, gridsize=40, mincnt=1, cmap="magma")
    axes[1].set_aspect("equal", adjustable="box"); axes[1].set(xlabel="Observed eastward tilt (km)", ylabel="Observed northward tilt (km)")
    axes[2].hist(prediction["angular_error"], bins=np.arange(0, 185, 5), color="slateblue", alpha=0.85)
    axes[2].axvline(90, color="k", ls="--", lw=1); axes[2].set(xlabel="Angular error (degrees)", ylabel="Count", xlim=(0, 180))
    fig.suptitle(title); fig.tight_layout(); return fig, axes


def plot_permutation_importance(importance, *, title):
    ordered = importance.sort_values("importance_mean")
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(ordered["feature"], ordered["importance_mean"], xerr=ordered["importance_std"], color="teal", alpha=0.8)
    ax.axvline(0, color="k", lw=0.8); ax.set_xlabel("Increase in equal-eddy vector MAE after permutation (km)"); ax.set_title(title)
    fig.tight_layout(); return fig, ax
