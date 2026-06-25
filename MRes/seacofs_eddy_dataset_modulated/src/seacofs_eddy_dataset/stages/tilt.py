from __future__ import annotations

from seacofs_eddy_dataset.config import PipelineConfig


def run(config: PipelineConfig) -> None:
    """Compute tilt metrics from confirmed vertical profile outputs.

    This will be extracted from `tilt_dataset.ipynb` and related tilt notebooks.
    """
    raise NotImplementedError("Move tilt dataset calculation logic here from tilt_dataset.ipynb")


def run_analysis(config: PipelineConfig) -> None:
    """Run final tilt analysis products.

    Analysis notebooks can remain notebooks, but reusable metrics and plotting data
    preparation should live here.
    """
    raise NotImplementedError("Move reusable tilt analysis logic here from tilt analysis notebooks")
