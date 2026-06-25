from __future__ import annotations

from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.core.doppio import bad_doppio_row, find_directional_radii, transect_indexer
from seacofs_eddy_dataset.core.esp import load_doppio_functions
from seacofs_eddy_dataset.core.grid import fnumber_from_outer_avg, read_reference_grid
from seacofs_eddy_dataset.core.velocity import rotate_uv
from seacofs_eddy_dataset.io import partition_path, write_partition
from seacofs_eddy_dataset.stages.detection import find_model_files, reference_grid_path


SURFACE_COLUMNS = [
    "Day",
    "fnumber",
    "nxc",
    "nyc",
    "nCyc",
    "nic",
    "njc",
    "xc",
    "yc",
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


def detection_path_for_file(path: Path, config: PipelineConfig) -> Path:
    return partition_path(config.output_root, "detections", f"fnumber={fnumber_from_outer_avg(path):05}")


def surface_output_path(path: Path, config: PipelineConfig) -> Path:
    return partition_path(config.output_root, "surface_eddies", f"fnumber={fnumber_from_outer_avg(path):05}")


def _row_with_q(row: dict, q) -> dict:
    out = dict(row)
    if np.asarray(q).shape == (2, 2):
        out["q11"] = float(q[0, 0])
        out["q12"] = float(q[0, 1])
        out["q22"] = float(q[1, 1])
    else:
        out["q11"] = np.nan
        out["q12"] = np.nan
        out["q22"] = np.nan
    out.pop("Q", None)
    return out


def _bad_row(row, fnumber: int) -> dict:
    base = bad_doppio_row(row, fnumber)
    base["nic"] = getattr(row, "nic", np.nan)
    base["njc"] = getattr(row, "njc", np.nan)
    return _row_with_q(base, np.nan)


def _fit_detection_row(row, fnumber: int, u_rot, v_rot, grid, doppio, out_core_param_fit, radius_km: float) -> dict:
    try:
        x1, y1, x2, y2, ii, jj = transect_indexer(int(row.nic), int(row.njc), grid.X_grid, grid.Y_grid, r=radius_km)
        if len(ii) < 3 or len(jj) < 3:
            return _bad_row(row, fnumber)
        xc, yc, w, q, omega0 = doppio(
            x1,
            y1,
            u_rot[ii, int(row.njc)],
            v_rot[ii, int(row.njc)],
            x2,
            y2,
            u_rot[int(row.nic), jj],
            v_rot[int(row.nic), jj],
        )
        radii = find_directional_radii(u_rot, v_rot, grid.X_grid, grid.Y_grid, xc, yc, q)
        radius = np.nanmean([radii["up"], radii["right"], radii["down"], radii["left"]])
        if not np.isfinite(radius):
            return _bad_row(row, fnumber)

        rho_limit = radius * 1.5
        local = (grid.X_grid >= xc - rho_limit) & (grid.X_grid <= xc + rho_limit)
        local &= (grid.Y_grid >= yc - rho_limit) & (grid.Y_grid <= yc + rho_limit)
        if int(local.sum()) < 10:
            return _bad_row(row, fnumber)

        xloc = grid.X_grid[local]
        yloc = grid.Y_grid[local]
        uloc = u_rot[local]
        vloc = v_rot[local]
        rho = np.sqrt(np.maximum(q[0, 0] * (xloc - xc) ** 2 + 2 * q[0, 1] * (xloc - xc) * (yloc - yc) + q[1, 1] * (yloc - yc) ** 2, 0))
        mask = np.isfinite(rho) & (rho <= rho_limit)
        if int(mask.sum()) < 10:
            return _bad_row(row, fnumber)
        rc, psi0, omega = out_core_param_fit(xloc[mask], yloc[mask], uloc[mask], vloc[mask], xc, yc, q, w)
    except Exception:
        return _bad_row(row, fnumber)

    return _row_with_q(
        {
            "Day": int(row.Day),
            "fnumber": fnumber,
            "nxc": float(row.nxc),
            "nyc": float(row.nyc),
            "nCyc": row.Cyc,
            "nic": int(row.nic),
            "njc": int(row.njc),
            "xc": float(xc),
            "yc": float(yc),
            "w": float(w),
            "Q": q,
            "Omega0": float(omega0),
            "Omega": float(omega),
            "Rc": float(rc),
            "psi0": float(psi0),
            "R": float(radius),
        },
        q,
    )


def fit_surface_file(path: Path, config: PipelineConfig) -> pd.DataFrame:
    fnumber = fnumber_from_outer_avg(path)
    detection_path = detection_path_for_file(path, config)
    if not detection_path.exists():
        raise FileNotFoundError(f"Missing detection partition for {path.name}: {detection_path}")

    detections = pd.read_parquet(detection_path)
    if detections.empty:
        return pd.DataFrame(columns=SURFACE_COLUMNS)

    grid = read_reference_grid(reference_grid_path(config))
    doppio, out_core_param_fit = load_doppio_functions(config)
    radius_km = float(config.raw.get("surface_fit", {}).get("transect_radius_km", 30.0))
    rows = []

    with nc.Dataset(path) as dataset:
        for t, day_value in enumerate(np.asarray(dataset.variables["ocean_time"][:].data, dtype=float) / 86400):
            day = int(round(day_value))
            day_detections = detections.loc[detections["Day"].eq(day)]
            if day_detections.empty:
                continue
            u = dataset["u_eastward"][t, -1, :, :].T
            v = dataset["v_northward"][t, -1, :, :].T
            u_rot, v_rot = rotate_uv(u, v, grid.angle)
            rows.extend(
                _fit_detection_row(row, fnumber, u_rot, v_rot, grid, doppio, out_core_param_fit, radius_km)
                for row in day_detections.itertuples(index=False)
            )

    return pd.DataFrame(rows, columns=SURFACE_COLUMNS)


def fit_and_write_file(path: Path, config: PipelineConfig) -> Path:
    out_path = surface_output_path(path, config)
    if config.skip_existing and out_path.exists():
        return out_path
    return write_partition(fit_surface_file(path, config), out_path)


def run(config: PipelineConfig) -> None:
    from joblib import Parallel, delayed

    files = find_model_files(config)
    backend = config.raw.get("parallel", {}).get("backend", "process")
    prefer = "threads" if backend == "thread" else "processes"
    written = Parallel(n_jobs=config.workers, prefer=prefer)(
        delayed(fit_and_write_file)(path, config) for path in files
    )
    for path in written:
        print(path)
