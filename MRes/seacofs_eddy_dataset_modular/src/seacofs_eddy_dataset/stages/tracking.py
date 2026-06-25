from __future__ import annotations

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.core.tracking import clean_surface_eddies as clean_surface_eddies_core


def clean_surface_eddies(config: PipelineConfig):
    """Apply QC before tracking.

    This will be extracted from `eddy_tracking.ipynb`, especially `clean_doppio`.
    """
    _ = clean_surface_eddies_core
    raise NotImplementedError("Move surface eddy QC logic here from eddy_tracking.ipynb")


def track_eddies(config: PipelineConfig):
    """Link cleaned surface eddies through time.

    This will wrap the existing clim_functions tracking helpers in a deterministic stage.
    """
    raise NotImplementedError("Move tracking logic here from eddy_tracking.ipynb")


def run(config: PipelineConfig) -> None:
    raise NotImplementedError("Tracking runner is scaffolded but not implemented yet")
