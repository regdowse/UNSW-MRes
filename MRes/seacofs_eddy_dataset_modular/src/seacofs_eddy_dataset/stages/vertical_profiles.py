from __future__ import annotations

from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.core.doppio import find_directional_radii, nearest_ij, transect_indexer
from seacofs_eddy_dataset.core.esp import load_doppio_functions
from seacofs_eddy_dataset.core.grid import fnumber_from_outer_avg, read_reference_grid
from seacofs_eddy_dataset.core.velocity import rotate_uv
from seacofs_eddy_dataset.core.vertical import interp_3d_to_reference_depths
from seacofs_eddy_dataset.io import read_partitions, write_partition
from seacofs_eddy_dataset.stages.detection import find_model_files, reference_grid_path


PROFILE_COLUMNS = [
    "Eddy",
    "Day",
    "fnumber",
    "Depth",
    "z",
    "xc",
    "yc",
    "ic",
    "jc",
    "w",
    "q11",
    "q12",
    "q22",
    "Omega0",
    "Omega",
    "Rc",
    "psi0",
    "R",
]


def _processed_dataset(config: PipelineConfig) -> pd.DataFrame:
    path = config.output_root / "processed" / "eddy_dataset_processed.parquet"
    if not path.exists():
        raise FileNotFoundError(f"Missing processed eddy dataset: {path}")
    return pd.read_parquet(path)


def vertical_output_path(path: Path, config: PipelineConfig) -> Path:
    return config.output_root / "vertical_profiles" / f"fnumber={fnumber_from_outer_avg(path):05}.parquet"


def _q_from_row(row) -> np.ndarray:
    return np.array([[row.q11, row.q12], [row.q12, row.q22]], dtype=float)


def _surface_row(row) -> dict:
    return {
        "Eddy": int(row.Eddy),
        "Day": int(row.Day),
        "fnumber": int(row.fnumber),
        "Depth": 0.0,
        "z": 0,
        "xc": float(row.xc),
        "yc": float(row.yc),
        "ic": int(row.ic),
        "jc": int(row.jc),
        "w": float(row.w),
        "q11": float(row.q11),
        "q12": float(row.q12),
        "q22": float(row.q22),
        "Omega0": float(row.Omega0),
        "Omega": float(row.Omega),
        "Rc": float(row.Rc),
        "psi0": float(row.psi0),
        "R": float(row.R),
    }


def _fit_depth(
    row,
    depth_idx: int,
    depth: float,
    u2d,
    v2d,
    grid,
    doppio,
    out_core_param_fit,
    radius_km: float,
    rho_max: float,
    rho_min: float,
    max_jump_km: float,
) -> dict | None:
    q_prev = _q_from_row(row)
    try:
        ic, jc = nearest_ij(float(row.xc), float(row.yc), grid.X_grid, grid.Y_grid)
        x1, y1, x2, y2, ii, jj = transect_indexer(ic, jc, grid.X_grid, grid.Y_grid, r=radius_km)
        if len(ii) < 3 or len(jj) < 3:
            return None
        xc, yc, w, q, omega0 = doppio(
            x1,
            y1,
            u2d[ii, jc],
            v2d[ii, jc],
            x2,
            y2,
            u2d[ic, jj],
            v2d[ic, jj],
        )
        if np.hypot(xc - row.xc, yc - row.yc) > max_jump_km:
            return None
        if np.sign(w) != np.sign(row.w):
            return None

        radii = find_directional_radii(u2d, v2d, grid.X_grid, grid.Y_grid, xc, yc, q)
        radius = np.nanmean([radii["up"], radii["right"], radii["down"], radii["left"]])
        if not np.isfinite(radius) or radius < rho_min:
            return None

        rho_limit = min(radius * 1.5, rho_max)
        local = (grid.X_grid >= xc - rho_limit) & (grid.X_grid <= xc + rho_limit)
        local &= (grid.Y_grid >= yc - rho_limit) & (grid.Y_grid <= yc + rho_limit)
        if int(local.sum()) < 10:
            return None

        xloc = grid.X_grid[local]
        yloc = grid.Y_grid[local]
        uloc = u2d[local]
        vloc = v2d[local]
        rho = np.sqrt(np.maximum(q_prev[0, 0] * (xloc - xc) ** 2 + 2 * q_prev[0, 1] * (xloc - xc) * (yloc - yc) + q_prev[1, 1] * (yloc - yc) ** 2, 0))
        mask = np.isfinite(rho) & (rho <= rho_limit)
        if int(mask.sum()) < 10:
            return None
        rc, psi0, omega = out_core_param_fit(xloc[mask], yloc[mask], uloc[mask], vloc[mask], xc, yc, q, w)
    except Exception:
        return None

    return {
        "Eddy": int(row.Eddy),
        "Day": int(row.Day),
        "fnumber": int(row.fnumber),
        "Depth": float(depth),
        "z": int(depth_idx),
        "xc": float(xc),
        "yc": float(yc),
        "ic": int(ic),
        "jc": int(jc),
        "w": float(w),
        "q11": float(q[0, 0]),
        "q12": float(q[0, 1]),
        "q22": float(q[1, 1]),
        "Omega0": float(omega0),
        "Omega": float(omega),
        "Rc": float(rc),
        "psi0": float(psi0),
        "R": float(radius),
    }


def compute_profiles_for_file(path: Path, config: PipelineConfig) -> pd.DataFrame:
    fnumber = fnumber_from_outer_avg(path)
    df = _processed_dataset(config)
    df_file = df.loc[df["fnumber"].eq(fnumber)].copy()
    if df_file.empty:
        return pd.DataFrame(columns=PROFILE_COLUMNS)

    settings = config.raw.get("vertical_profiles", {})
    grid = read_reference_grid(reference_grid_path(config))
    doppio, out_core_param_fit = load_doppio_functions(config)
    z_r_path = Path(config.raw["paths"]["z_r"]).expanduser()
    z_r = np.load(z_r_path)
    if z_r.shape[:2] != grid.X_grid.shape and z_r.shape[-2:] == grid.X_grid.shape:
        z_r = np.moveaxis(z_r, 0, -1)

    configured_depths = settings.get("target_depths")
    if configured_depths is None:
        target_depths = np.nanmedian(np.abs(z_r), axis=(0, 1))
    else:
        target_depths = np.asarray(configured_depths, dtype=float)
    max_depth_levels = settings.get("max_depth_levels")
    if max_depth_levels is not None:
        target_depths = target_depths[: int(max_depth_levels)]

    rows: list[dict] = []
    with nc.Dataset(path) as dataset:
        days = np.asarray(dataset.variables["ocean_time"][:].data, dtype=float) / 86400
        for t, day_value in enumerate(days):
            day = int(round(day_value))
            df_day = df_file.loc[df_file["Day"].eq(day)]
            if df_day.empty:
                continue
            u3d = np.transpose(dataset["u_eastward"][t, :, :, :], (2, 1, 0))
            v3d = np.transpose(dataset["v_northward"][t, :, :, :], (2, 1, 0))
            u_depth = interp_3d_to_reference_depths(u3d, z_r, target_depths)
            v_depth = interp_3d_to_reference_depths(v3d, z_r, target_depths)

            for row in df_day.itertuples(index=False):
                rows.append(_surface_row(row))
                for depth_idx, depth in enumerate(target_depths):
                    if depth <= 0:
                        continue
                    u2d, v2d = rotate_uv(u_depth[:, :, depth_idx], v_depth[:, :, depth_idx], grid.angle)
                    fitted = _fit_depth(
                        row,
                        depth_idx,
                        depth,
                        u2d,
                        v2d,
                        grid,
                        doppio,
                        out_core_param_fit,
                        radius_km=float(settings.get("transect_radius_km", 30.0)),
                        rho_max=float(settings.get("rho_max_km", 200.0)),
                        rho_min=float(settings.get("rho_min_km", 30.0)),
                        max_jump_km=float(settings.get("max_jump_km", 100.0)),
                    )
                    if fitted is None:
                        break
                    rows.append(fitted)

    return pd.DataFrame(rows, columns=PROFILE_COLUMNS)


def compute_and_write_file(path: Path, config: PipelineConfig) -> Path:
    out_path = vertical_output_path(path, config)
    if config.skip_existing and out_path.exists():
        return out_path
    return write_partition(compute_profiles_for_file(path, config), out_path)


def run(config: PipelineConfig) -> None:
    from joblib import Parallel, delayed

    files = find_model_files(config)
    backend = config.raw.get("parallel", {}).get("backend", "process")
    prefer = "threads" if backend == "thread" else "processes"
    written = Parallel(n_jobs=config.workers, prefer=prefer)(
        delayed(compute_and_write_file)(path, config) for path in files
    )
    for path in written:
        print(path)


def run_qc(config: PipelineConfig) -> None:
    source_paths = sorted((config.output_root / "vertical_profiles").glob("fnumber=*.parquet"))
    profiles = read_partitions(source_paths)
    out_path = config.output_root / "vertical_profiles_confirmed" / "profiles.parquet"
    if profiles.empty:
        write_partition(profiles, out_path)
        print(out_path)
        return

    settings = config.raw.get("tilt", {})
    min_count = int(settings.get("min_profile_depth_count", 3))
    counts = profiles.groupby(["Eddy", "Day"])["Depth"].transform("count")
    confirmed = profiles.loc[counts >= min_count].copy()
    write_partition(confirmed, out_path)
    print(out_path)
