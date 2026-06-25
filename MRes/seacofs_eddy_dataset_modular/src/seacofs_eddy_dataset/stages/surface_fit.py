from __future__ import annotations

from pathlib import Path

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.core.doppio import bad_doppio_row, find_directional_radii, transect_indexer
from seacofs_eddy_dataset.core.velocity import rotate_uv


def fit_surface_file(path: Path, config: PipelineConfig):
    """Fit DOPPIO surface eddy parameters for one model file.

    This stage should read one detection partition, rotate surface velocities,
    extract transects with `core.doppio`, and call the external DOPPIO fit.

    Expected output columns include Day, fnumber, nxc, nyc, nCyc, xc, yc, w, Q,
    Omega0, Omega, Rc, psi0, and R.
    """
    _ = (bad_doppio_row, find_directional_radii, rotate_uv, transect_indexer)
    raise NotImplementedError("Move DOPPIO surface fitting logic here from doppio_dataset.ipynb")


def run(config: PipelineConfig) -> None:
    raise NotImplementedError("Parallel DOPPIO surface fit runner is scaffolded but not implemented yet")
