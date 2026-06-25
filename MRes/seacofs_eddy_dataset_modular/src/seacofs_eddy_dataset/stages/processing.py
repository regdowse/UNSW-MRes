from __future__ import annotations

from seacofs_eddy_dataset.config import PipelineConfig


def run(config: PipelineConfig) -> None:
    """Clean and normalise tracked eddy data for vertical-profile extraction.

    This will be extracted from `dataset_processing.ipynb`.
    """
    raise NotImplementedError("Move tracked dataset processing logic here from dataset_processing.ipynb")
