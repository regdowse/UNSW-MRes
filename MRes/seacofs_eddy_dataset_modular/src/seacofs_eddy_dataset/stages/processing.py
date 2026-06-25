from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.core.grid import read_reference_grid
from seacofs_eddy_dataset.core.doppio import nearest_ij
from seacofs_eddy_dataset.io import write_partition
from seacofs_eddy_dataset.stages.detection import reference_grid_path


def _add_lon_lat(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    grid = read_reference_grid(reference_grid_path(config))
    out = df.copy()
    points = np.column_stack((out["xc"].to_numpy(float), out["yc"].to_numpy(float)))
    lon_interp = RegularGridInterpolator((grid.x_grid, grid.y_grid), grid.lon_rho, bounds_error=False, fill_value=np.nan)
    lat_interp = RegularGridInterpolator((grid.x_grid, grid.y_grid), grid.lat_rho, bounds_error=False, fill_value=np.nan)
    out["lon"] = lon_interp(points)
    out["lat"] = lat_interp(points)
    ij = [nearest_ij(xc, yc, grid.X_grid, grid.Y_grid) for xc, yc in zip(out["xc"], out["yc"], strict=False)]
    out["ic"] = [i for i, _ in ij]
    out["jc"] = [j for _, j in ij]
    return out


def run(config: PipelineConfig) -> None:
    out_path = config.output_root / "processed" / "eddy_dataset_processed.parquet"
    if config.skip_existing and out_path.exists():
        print(out_path)
        return

    source = config.output_root / "tracked" / "eddy_tracks.parquet"
    if not source.exists():
        raise FileNotFoundError(f"Missing tracked dataset: {source}")

    df = pd.read_parquet(source)
    if df.empty:
        write_partition(df, out_path)
        print(out_path)
        return

    out = _add_lon_lat(df, config)
    out["fname"] = out["fnumber"].map(lambda n: f"outer_avg_{int(n):05}.nc")
    out = out.sort_values(["Eddy", "Day"]).reset_index(drop=True)
    write_partition(out, out_path)
    print(out_path)
