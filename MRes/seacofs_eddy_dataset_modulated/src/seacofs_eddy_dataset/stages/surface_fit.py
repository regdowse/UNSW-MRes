from __future__ import annotations

from pathlib import Path

from seacofs_eddy_dataset.config import PipelineConfig


def fit_surface_file(path: Path, config: PipelineConfig):
    """Fit DOPPIO surface eddy parameters for one model file.

    This will be extracted from `doppio_dataset.ipynb`.
    Expected output columns include Day, fnumber, nxc, nyc, nCyc, xc, yc, w, Q,
    Omega0, Omega, Rc, psi0, and R.
    """
    raise NotImplementedError("Move DOPPIO surface fitting logic here from doppio_dataset.ipynb")


def run(config: PipelineConfig) -> None:
    raise NotImplementedError("Parallel DOPPIO surface fit runner is scaffolded but not implemented yet")
