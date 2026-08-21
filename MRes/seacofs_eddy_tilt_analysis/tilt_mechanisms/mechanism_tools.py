"""Helpers unique to the EAC eddy-tilt mechanism notebooks.

Tilt geometry is never re-estimated here. ``TiltDis`` and ``TiltDir`` are
treated as authoritative measurements loaded by ``seacofs_tilt_tools``.
Existing project helpers remain the source of grid, PV, tracking, bootstrap,
and background-flow calculations.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path

import numpy as np
import pandas as pd


SECONDS_PER_DAY = 86400.0
METRES_PER_KM = 1000.0


def require_tilt_measurements(df: pd.DataFrame) -> pd.DataFrame:
    """Validate, but never calculate, the existing tilt measurements."""

    missing = {"Eddy", "Day", "Cyc", "TiltDis", "TiltDir"} - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")
    out = df.copy()
    out.loc[~out["TiltDir"].between(0, 360, inclusive="left"), "TiltDir"] = np.nan
    out.loc[out["TiltDis"] < 0, "TiltDis"] = np.nan
    return out


def add_tilt_components(df: pd.DataFrame) -> pd.DataFrame:
    """Encode measured tilt as east/north components (bearing is from north)."""

    out = require_tilt_measurements(df)
    theta = np.deg2rad(out["TiltDir"])
    out["tilt_east_km"] = out["TiltDis"] * np.sin(theta)
    out["tilt_north_km"] = out["TiltDis"] * np.cos(theta)
    return out


def bearing_from_east_north(east, north):
    return np.degrees(np.arctan2(east, north)) % 360.0


def signed_angle_difference(observed, reference):
    """Return observed minus reference in [-180, 180) degrees."""

    return (np.asarray(observed) - np.asarray(reference) + 180.0) % 360.0 - 180.0


def merge_one_to_one_or_many_to_one(left, right, keys=("Eddy", "Day")):
    """Merge diagnostic caches while guarding against duplicated eddy-days."""

    keys = list(keys)
    if right.duplicated(keys).any():
        duplicated = right.loc[right.duplicated(keys, keep=False), keys].head()
        raise ValueError(f"Right table contains duplicate keys:\n{duplicated}")
    return left.merge(right, on=keys, how="left", validate="many_to_one")


BACKGROUND_PAIRS = {
    "ann_200": ("ann_surface", "ann_200"),
    "ann_500": ("ann_surface", "ann_500"),
    "clim_200": ("clim_surface", "clim_200"),
    "clim_500": ("clim_surface", "clim_500"),
    "full_200": ("full_surface", "full_200"),
    "full_500": ("full_surface", "full_500"),
}


def add_background_shear(df: pd.DataFrame, pairs=BACKGROUND_PAIRS) -> pd.DataFrame:
    """Add surface-minus-depth shear and its alignment with measured tilt."""

    out = add_tilt_components(df)
    for label, (surface, deep) in pairs.items():
        east = f"{label}_shear_east_ms"
        north = f"{label}_shear_north_ms"
        required = [
            f"{surface}_east_ms", f"{surface}_north_ms",
            f"{deep}_east_ms", f"{deep}_north_ms",
        ]
        missing = set(required) - set(out.columns)
        if missing:
            raise KeyError(f"Cannot calculate {label}; missing {sorted(missing)}")
        out[east] = out[f"{surface}_east_ms"] - out[f"{deep}_east_ms"]
        out[north] = out[f"{surface}_north_ms"] - out[f"{deep}_north_ms"]
        out[f"{label}_shear_mag_ms"] = np.hypot(out[east], out[north])
        out[f"{label}_shear_dir"] = bearing_from_east_north(out[east], out[north])
        out[f"{label}_tilt_shear_offset"] = signed_angle_difference(
            out["TiltDir"], out[f"{label}_shear_dir"]
        )
    return out


def _rolling_vector_integral(part, east, north, window_days):
    """Integrate a vector velocity over trailing windows using trapezoids."""

    part = part.sort_values("Day")
    days = part["Day"].to_numpy(float)
    ue = part[east].to_numpy(float)
    vn = part[north].to_numpy(float)
    ie = np.full(len(part), np.nan)
    inn = np.full(len(part), np.nan)
    for stop in range(len(part)):
        use = (days >= days[stop] - window_days) & (days <= days[stop])
        use &= np.isfinite(days) & np.isfinite(ue) & np.isfinite(vn)
        idx = np.flatnonzero(use)
        if idx.size < 2:
            continue
        t = days[idx] * SECONDS_PER_DAY
        ie[stop] = np.trapz(ue[idx], t) / METRES_PER_KM
        inn[stop] = np.trapz(vn[idx], t) / METRES_PER_KM
    return pd.DataFrame({"east": ie, "north": inn}, index=part.index)


def add_accumulated_shear(df, method="ann_500", windows=(5, 10, 20, 30)):
    """Add trailing integrals of surface-minus-depth background velocity."""

    out = df.sort_values(["Eddy", "Day"]).copy()
    east = f"{method}_shear_east_ms"
    north = f"{method}_shear_north_ms"
    for window in windows:
        accum = pd.concat(
            [_rolling_vector_integral(part, east, north, window)
             for _, part in out.groupby("Eddy", sort=False)]
        ).sort_index()
        prefix = f"{method}_accum_{int(window)}d"
        out[f"{prefix}_east_km"] = accum["east"]
        out[f"{prefix}_north_km"] = accum["north"]
        out[f"{prefix}_mag_km"] = np.hypot(accum["east"], accum["north"])
        out[f"{prefix}_dir"] = bearing_from_east_north(accum["east"], accum["north"])
        out[f"{prefix}_offset"] = signed_angle_difference(out["TiltDir"], out[f"{prefix}_dir"])
    return out


def circular_offset_summary(df, columns, group=("Cyc",), eddy_equal=True):
    """Summarise signed angular offsets with optional equal eddy weighting."""

    records = []
    group = list(group)
    for keys, part in df.groupby(group, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        for column in columns:
            values = part[["Eddy", column]].dropna()
            if eddy_equal:
                from seacofs_tilt_tools import circular_mean_deg_true_north

                values = values.groupby("Eddy")[column].apply(circular_mean_deg_true_north).reset_index()
                angles = values[column].to_numpy(float)
            else:
                angles = values[column].to_numpy(float)
            if not len(angles):
                continue
            radians = np.deg2rad(angles)
            vector = np.nanmean(np.exp(1j * radians))
            record = dict(zip(group, keys))
            record.update({
                "metric": column,
                "eddies": int(values["Eddy"].nunique()),
                "mean_offset_deg": float(np.degrees(np.angle(vector))),
                "resultant_length": float(np.abs(vector)),
                "median_abs_offset_deg": float(np.nanmedian(np.abs(angles))),
            })
            records.append(record)
    return pd.DataFrame(records)


def add_stratification_proxies(df: pd.DataFrame) -> pd.DataFrame:
    """Add N2/f2 and transparent constant-N deformation/Burger proxies."""

    out = df.copy()
    if "f" not in out:
        raise KeyError("The table needs local Coriolis parameter `f`.")
    for depth in (200, 500):
        n2 = f"N2_{depth}m_s2"
        if n2 not in out:
            continue
        out[f"N2_over_f2_{depth}m"] = out[n2] / out["f"].pow(2)
        # Constant-N, flat-depth first-mode proxy: Rd ~ N H / (pi |f|).
        out[f"Rd_proxy_{depth}m_km"] = (
            np.sqrt(out[n2].clip(lower=0)) * depth
            / (np.pi * out["f"].abs()) / METRES_PER_KM
        )
        out[f"Bu_proxy_{depth}m"] = (
            out[f"Rd_proxy_{depth}m_km"] / out["Rc"]
        ).pow(2)
    return out


def _day_values(time):
    values = np.asarray(time.values)
    if np.issubdtype(values.dtype, np.datetime64):
        return ((values - np.datetime64("1990-01-01")) / np.timedelta64(1, "D")).astype(float)
    units = str(time.attrs.get("units", ""))
    scale = SECONDS_PER_DAY if "second" in units.lower() else 1.0
    return values.astype(float) / scale


def _depth_means(n2_columns, z_columns, depths):
    """Depth-average point columns; input arrays have shape point x vertical."""

    n2_columns = np.asarray(n2_columns, float)
    z_columns = np.asarray(z_columns, float)
    if n2_columns.shape != z_columns.shape or n2_columns.ndim != 2:
        raise ValueError("N2 and z must be matching point-by-vertical arrays.")
    output = {int(depth): np.full(n2_columns.shape[0], np.nan) for depth in depths}
    for point, (n2, z) in enumerate(zip(n2_columns, z_columns)):
        for depth in depths:
            use = np.isfinite(n2) & np.isfinite(z) & (z <= 0) & (z >= -float(depth))
            if use.sum() < 2:
                continue
            order = np.argsort(z[use])
            zz = z[use][order]
            nn = n2[use][order]
            span = zz[-1] - zz[0]
            if span > 0:
                output[int(depth)][point] = np.trapz(nn, zz) / span
    return output


@dataclass(frozen=True)
class N2CacheConfig:
    """Runtime and restart settings for the file-parallel N2 cache."""

    model_root: Path = Path("/srv/scratch/z3533156/26year_BRAN2020")
    output_path: Path = Path(
        "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/"
        "tilt_mechanisms/n2_eddy_day.parquet"
    )
    depths: tuple[int, ...] = (200, 500)
    rho0: float = 1025.0
    horizontal_chunk: int = 64
    skip_existing: bool = True

    @property
    def partition_root(self):
        return self.output_path.parent / "n2_file_partitions"


def _n2_partition_path(model_path, config):
    return config.partition_root / f"{Path(model_path).stem}.parquet"


def _atomic_parquet(table, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f".tmp-{os.getpid()}.parquet")
    table.to_parquet(temporary, index=False)
    os.replace(temporary, path)


def process_n2_model_file(model_path, file_rows, config=N2CacheConfig()):
    """Calculate all requested eddy columns in one model file and partition."""

    import xarray as xr
    import xroms

    model_path = Path(model_path)
    partition = _n2_partition_path(model_path, config)
    if config.skip_existing and partition.exists():
        return str(partition), "cached"

    rows = file_rows[["Eddy", "Day", "ic", "jc"]].copy()
    rows["Day"] = rows["Day"].round().astype(int)
    rows["ic"] = rows["ic"].astype(int)
    rows["jc"] = rows["jc"].astype(int)
    rows = rows.drop_duplicates(["Eddy", "Day"])
    output = []
    chunks = {
        "ocean_time": 1,
        "s_rho": -1,
        "eta_rho": config.horizontal_chunk,
        "xi_rho": config.horizontal_chunk,
    }
    with xr.open_dataset(model_path, chunks=chunks) as raw:
        ds, xgrid = xroms.roms_dataset(raw)
        rho = xroms.density(ds["temp"], ds["salt"])
        n2 = xroms.N2(rho, xgrid, rho0=config.rho0)
        z = n2.coords.get("z_w", ds.get("z_w"))
        if z is None:
            raise KeyError("xroms did not provide z_w for the N2 field.")
        vertical_dims = [dim for dim in n2.dims if dim.startswith("s_")]
        if len(vertical_dims) != 1:
            raise ValueError(f"Expected one N2 vertical dimension, found {vertical_dims}.")
        vertical_dim = vertical_dims[0]
        model_days = np.rint(_day_values(ds["ocean_time"])).astype(int)
        time_for_day = {int(day): t for t, day in enumerate(model_days)}

        for day, day_rows in rows.groupby("Day", sort=False):
            if int(day) not in time_for_day:
                raise KeyError(f"Day {day} was assigned to {model_path.name} but is absent from ocean_time.")
            # One vectorized point selection and one dask compute per timestep.
            xi = xr.DataArray(day_rows["ic"].to_numpy(), dims="point")
            eta = xr.DataArray(day_rows["jc"].to_numpy(), dims="point")
            selector = {"ocean_time": time_for_day[int(day)], "xi_rho": xi, "eta_rho": eta}
            selected = xr.Dataset({"N2": n2.isel(selector), "z": z.isel(selector)}).compute(
                scheduler="synchronous"
            )
            n2_values = selected["N2"].transpose("point", vertical_dim).values
            z_values = selected["z"].transpose("point", vertical_dim).values
            means = _depth_means(n2_values, z_values, config.depths)
            for position, row in enumerate(day_rows.itertuples(index=False)):
                record = {"Eddy": row.Eddy, "Day": int(day)}
                for depth in config.depths:
                    record[f"N2_{int(depth)}m_s2"] = means[int(depth)][position]
                output.append(record)

    columns = ["Eddy", "Day"] + [f"N2_{int(depth)}m_s2" for depth in config.depths]
    table = pd.DataFrame(output, columns=columns).sort_values(["Eddy", "Day"])
    _atomic_parquet(table, partition)
    return str(partition), "computed"


def build_n2_cache_xroms(
    eddies: pd.DataFrame,
    model_root: Path | str | None = None,
    output_path: Path | str | None = None,
    *,
    depths=(200, 500),
    rho0=1025.0,
    workers=4,
    horizontal_chunk=64,
    skip_existing=True,
):
    """Build a restartable, model-file-parallel eddy-centre N2 cache.

    The calculation follows the documented xroms sequence:
    ``roms_dataset -> density -> N2``. Results are written once as
    a parquet cache. ``ic`` indexes xi and ``jc`` indexes eta, matching the
    transposed arrays used by ``seacofs_tilt_tools``.
    """

    from joblib import Parallel, delayed
    from beta_effect_background_flow.background_flow_tools import model_file_for_day

    required = {"Eddy", "Day", "ic", "jc"}
    missing = required - set(eddies.columns)
    if missing:
        raise KeyError(f"Missing N2 extraction columns: {sorted(missing)}")
    defaults = N2CacheConfig()
    config = N2CacheConfig(
        model_root=Path(model_root) if model_root is not None else defaults.model_root,
        output_path=Path(output_path) if output_path is not None else defaults.output_path,
        depths=tuple(int(depth) for depth in depths),
        rho0=float(rho0),
        horizontal_chunk=int(horizontal_chunk),
        skip_existing=bool(skip_existing),
    )
    work = eddies[["Eddy", "Day", "ic", "jc"]].dropna().drop_duplicates(["Eddy", "Day"]).copy()
    work["model_file"] = work["Day"].map(
        lambda day: model_file_for_day(int(round(day)), config.model_root)
    )
    file_rows = {Path(path): part.drop(columns="model_file") for path, part in work.groupby("model_file")}
    missing_files = [path for path in file_rows if not path.exists()]
    if missing_files:
        raise FileNotFoundError(f"Missing {len(missing_files)} model files; first is {missing_files[0]}")

    results = Parallel(n_jobs=int(workers), prefer="processes", verbose=10)(
        delayed(process_n2_model_file)(path, rows, config)
        for path, rows in sorted(file_rows.items())
    )
    partition_paths = [path for path, _ in results]
    result = pd.concat([pd.read_parquet(path) for path in partition_paths], ignore_index=True)
    result = result.sort_values(["Eddy", "Day"]).reset_index(drop=True)
    if result.duplicated(["Eddy", "Day"]).any():
        raise ValueError("N2 output contains duplicate Eddy-Day rows.")
    expected = work[["Eddy", "Day"]].copy()
    expected["Day"] = expected["Day"].round().astype(int)
    key_check = expected.merge(
        result[["Eddy", "Day"]], on=["Eddy", "Day"], how="outer", indicator=True
    )
    if not key_check["_merge"].eq("both").all():
        counts = key_check["_merge"].value_counts().to_dict()
        raise ValueError(f"Reduced cache does not match requested Eddy-Day keys: {counts}")
    config.output_path.parent.mkdir(parents=True, exist_ok=True)
    _atomic_parquet(result, config.output_path)
    return result


def add_topographic_regimes(df: pd.DataFrame, shelf_depth=2000.0, dominance_ratio=1.0):
    """Attach readable shelf and planetary/topographic PV regimes."""

    out = df.copy()
    if "Region" in out:
        out["ShelfRegime"] = np.where(out["Region"].isin(["S1", "S2"]), "on_shelf", "off_shelf")
    else:
        out["ShelfRegime"] = np.where(out["h"] <= shelf_depth, "on_shelf", "off_shelf")
    ratio = out["PV_grad_topo_mag"] / out["PV_grad_plan_mag"]
    out["PVRegime"] = np.select(
        [ratio < 1 / dominance_ratio, ratio > dominance_ratio],
        ["planetary", "topographic"],
        default="mixed",
    )
    out["topo_plan_ratio_raw"] = ratio
    return out


def add_ekman_transport(df: pd.DataFrame, tau_east="tau_east_pa", tau_north="tau_north_pa", rho0=1025.0):
    """Calculate depth-integrated Ekman transport from geographic wind stress."""

    out = df.copy()
    required = {tau_east, tau_north, "f"}
    if missing := required - set(out.columns):
        raise KeyError(f"Missing wind columns: {sorted(missing)}")
    out["ekman_east_m2s"] = out[tau_north] / (rho0 * out["f"])
    out["ekman_north_m2s"] = -out[tau_east] / (rho0 * out["f"])
    out["ekman_dir"] = bearing_from_east_north(out["ekman_east_m2s"], out["ekman_north_m2s"])
    out["tilt_ekman_offset"] = signed_angle_difference(out["TiltDir"], out["ekman_dir"])
    return out


def standardise_within_between(df, columns, group="Eddy"):
    """Add eddy-mean and daily-anomaly versions of time-varying predictors."""

    out = df.copy()
    for column in columns:
        mean = out.groupby(group)[column].transform("mean")
        out[f"{column}_between"] = mean
        out[f"{column}_within"] = out[column] - mean
    return out
