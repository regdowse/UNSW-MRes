"""Helpers unique to the EAC eddy-tilt mechanism notebooks.

Tilt geometry is never re-estimated here. ``TiltDis`` and ``TiltDir`` are
treated as authoritative measurements loaded by ``seacofs_tilt_tools``.
Existing project helpers remain the source of grid, PV, tracking, bootstrap,
and background-flow calculations.
"""

from __future__ import annotations

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


def _column_depth_mean(n2_column, z_column, max_depth):
    n2 = np.asarray(n2_column, float)
    z = np.asarray(z_column, float)
    use = np.isfinite(n2) & np.isfinite(z) & (z <= 0) & (z >= -float(max_depth))
    if use.sum() < 2:
        return np.nan
    order = np.argsort(z[use])
    zz = z[use][order]
    nn = n2[use][order]
    span = zz[-1] - zz[0]
    return float(np.trapz(nn, zz) / span) if span > 0 else np.nan


def build_n2_cache_xroms(
    eddies: pd.DataFrame,
    model_root: Path | str,
    output_path: Path | str,
    *,
    depths=(200, 500),
    rho0=1025.0,
):
    """Extract depth-mean N2 at each eddy centre using xroms.

    The calculation follows the documented xroms sequence:
    ``roms_dataset -> density -> N2``. Results are written once as
    a parquet cache. ``ic`` indexes xi and ``jc`` indexes eta, matching the
    transposed arrays used by ``seacofs_tilt_tools``.
    """

    import xarray as xr
    import xroms
    from beta_effect_background_flow.background_flow_tools import model_file_for_day

    required = {"Eddy", "Day", "ic", "jc"}
    missing = required - set(eddies.columns)
    if missing:
        raise KeyError(f"Missing N2 extraction columns: {sorted(missing)}")
    work = eddies[list(required)].dropna().drop_duplicates(["Eddy", "Day"]).copy()
    work["model_file"] = work["Day"].map(lambda day: model_file_for_day(int(day), Path(model_root)))
    rows = []
    for model_file, file_rows in work.groupby("model_file", sort=True):
        with xr.open_dataset(model_file, chunks={"ocean_time": 1}) as raw:
            ds, xgrid = xroms.roms_dataset(raw)
            rho = xroms.density(ds["temp"], ds["salt"])
            n2 = xroms.N2(rho, xgrid, rho0=rho0)
            z = n2.coords.get("z_w", ds.get("z_w"))
            if z is None:
                raise KeyError("xroms did not provide z_w for the N2 field.")
            model_days = _day_values(ds["ocean_time"])
            for row in file_rows.itertuples(index=False):
                t = int(np.nanargmin(np.abs(model_days - float(row.Day))))
                selector = {"ocean_time": t, "xi_rho": int(row.ic), "eta_rho": int(row.jc)}
                n2_col = n2.isel(selector).compute().values
                z_col = z.isel(selector).compute().values
                record = {"Eddy": row.Eddy, "Day": row.Day}
                for depth in depths:
                    record[f"N2_{int(depth)}m_s2"] = _column_depth_mean(n2_col, z_col, depth)
                rows.append(record)
    result = pd.DataFrame(rows).sort_values(["Eddy", "Day"])
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path, index=False)
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
