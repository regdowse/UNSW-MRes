from __future__ import annotations

from seacofs_eddy_dataset.config import PipelineConfig


def clean_surface_eddies(config: PipelineConfig):
    """Apply QC before tracking.

    This will be extracted from `eddy_tracking.ipynb`, especially `clean_doppio`.
    """
    raise NotImplementedError("Move surface eddy QC logic here from eddy_tracking.ipynb")


def track_eddies(config: PipelineConfig):
    """Link cleaned surface eddies through time.

    This will wrap the existing clim_functions tracking helpers in a deterministic stage.
    """
    raise NotImplementedError("Move tracking logic here from eddy_tracking.ipynb")


def run(config: PipelineConfig) -> None:
    raise NotImplementedError("Tracking runner is scaffolded but not implemented yet")
