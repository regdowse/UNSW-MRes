from __future__ import annotations

from pathlib import Path

from seacofs_eddy_dataset.config import PipelineConfig


def compute_profiles_for_file(path: Path, config: PipelineConfig):
    """Compute vertical profiles for all tracked eddy-days in one model file.

    This will be extracted from `vert_dataset.ipynb` and should write one atomic
    partition per work unit rather than appending to one large pickle.
    """
    raise NotImplementedError("Move vertical profile logic here from vert_dataset.ipynb")


def run(config: PipelineConfig) -> None:
    raise NotImplementedError("Vertical profile runner is scaffolded but not implemented yet")


def run_qc(config: PipelineConfig) -> None:
    """Confirm usable vertical profile outputs.

    This replaces the current manual `vert_check_parallised` and
    `vert_confirmed_vert_dataset_parallised` notebook workflow.
    """
    raise NotImplementedError("Move vertical profile QC and confirmation logic here")
