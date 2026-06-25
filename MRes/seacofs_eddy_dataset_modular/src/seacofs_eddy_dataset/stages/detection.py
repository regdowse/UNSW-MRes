from __future__ import annotations

from pathlib import Path

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.core.grid import fnumber_from_outer_avg
from seacofs_eddy_dataset.core.nencioli import nencioli
from seacofs_eddy_dataset.core.velocity import interpolate_uv_to_grid


def find_model_files(config: PipelineConfig) -> list[Path]:
    pattern = config.raw.get("files", {}).get("model_pattern", "outer_avg_*.nc")
    return sorted(config.model_root.glob(pattern))


def detect_file(path: Path, config: PipelineConfig):
    """Detect Nencioli eddy candidates for one model file.

    This stage should load one `outer_avg_*.nc`, use `core.velocity` to rotate
    and interpolate surface velocities, then call `core.nencioli.nencioli`.

    Expected output columns include Day, nxc, nyc, Cyc, nic, njc, and fnumber.
    """
    _ = (fnumber_from_outer_avg, interpolate_uv_to_grid, nencioli)
    raise NotImplementedError("Move Nencioli detection logic here from nencioli_dataset.ipynb")


def run(config: PipelineConfig) -> None:
    files = find_model_files(config)
    if not files:
        raise FileNotFoundError(f"No model files found in {config.model_root}")
    raise NotImplementedError("Parallel detection runner is scaffolded but not implemented yet")
