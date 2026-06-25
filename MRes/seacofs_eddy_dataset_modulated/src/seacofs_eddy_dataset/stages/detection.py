from __future__ import annotations

from pathlib import Path

from seacofs_eddy_dataset.config import PipelineConfig


def find_model_files(config: PipelineConfig) -> list[Path]:
    pattern = config.raw.get("files", {}).get("model_pattern", "outer_avg_*.nc")
    return sorted(config.model_root.glob(pattern))


def detect_file(path: Path, config: PipelineConfig):
    """Detect Nencioli eddy candidates for one model file.

    This will be extracted from `nencioli_dataset.ipynb`.
    Expected output columns include Day, nxc, nyc, Cyc, nic, njc, and fnumber.
    """
    raise NotImplementedError("Move Nencioli detection logic here from nencioli_dataset.ipynb")


def run(config: PipelineConfig) -> None:
    files = find_model_files(config)
    if not files:
        raise FileNotFoundError(f"No model files found in {config.model_root}")
    raise NotImplementedError("Parallel detection runner is scaffolded but not implemented yet")
