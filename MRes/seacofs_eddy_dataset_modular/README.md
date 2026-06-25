# SEACOFS Eddy Dataset Modular Pipeline

This folder is a fresh modular version of the SEACOFS eddy dataset workflow.

The goal is to replace the current notebook-driven workflow with a reproducible, parallel pipeline that can run each stage in order, cache intermediate outputs, and resume safely after failures.

## Pipeline Stages

1. `detect_nencioli`
   - Reads `outer_avg_*.nc` model files.
   - Detects surface eddy candidates with the Nencioli method.
   - Writes partitioned detection tables by `fnumber`.

2. `fit_doppio_surface`
   - Reads Nencioli detections and the matching model file.
   - Fits DOPPIO surface eddy geometry and parameters.
   - Writes partitioned surface eddy tables by `fnumber`.

3. `track_eddies`
   - Reads all surface eddy partitions.
   - Applies QC and links eddies through time.
   - Writes the tracked eddy dataset.

4. `process_tracked_dataset`
   - Cleans and normalises the tracked dataset for downstream work.
   - Writes the processed eddy dataset.

5. `compute_vertical_profiles`
   - Reads the processed eddy dataset and 3D velocity fields.
   - Computes vertical profiles for each eddy-day.
   - Writes partitioned vertical profile outputs.

6. `qc_vertical_profiles`
   - Checks and confirms usable vertical profiles.
   - Writes the confirmed vertical profile dataset.

7. `compute_tilt`
   - Reads confirmed vertical profiles.
   - Computes tilt metrics and supporting diagnostics.
   - Writes the tilt dataset.

8. `analyse_tilt`
   - Reads the tilt dataset.
   - Produces figures, tables, and summary statistics.

## Design Principles

- Keep science kernels in importable Python modules, not notebooks.
- Keep notebooks for exploration, QC, and figures only.
- Write one output partition per independent work unit, usually `fnumber`, `day`, or `eddy_id`.
- Avoid repeatedly rewriting one large pickle during long jobs.
- Make every stage resumable by checking whether expected partition outputs already exist.
- Store paths, thresholds, and worker counts in config files.

## First Extraction Targets

The current notebooks map onto this folder roughly as follows:

- `nencioli_dataset.ipynb` -> `src/seacofs_eddy_dataset/stages/detection.py`
- `doppio_dataset.ipynb` -> `src/seacofs_eddy_dataset/stages/surface_fit.py`
- `eddy_tracking.ipynb` -> `src/seacofs_eddy_dataset/stages/tracking.py`
- `dataset_processing.ipynb` -> `src/seacofs_eddy_dataset/stages/processing.py`
- `vert_dataset.ipynb` -> `src/seacofs_eddy_dataset/stages/vertical_profiles.py`
- `tilt_dataset.ipynb` -> `src/seacofs_eddy_dataset/stages/tilt.py`

## Function Modules

Reusable functions live under `src/seacofs_eddy_dataset/core/`:

- `grid.py`: reference grid loading, Cartesian grid construction, filename helpers.
- `velocity.py`: fill-value cleaning, velocity rotation, interpolation.
- `nencioli.py`: Nencioli eddy detection kernel.
- `doppio.py`: transect, radius, local eddy geometry, and DOPPIO table helpers.
- `tracking.py`: surface eddy QC before temporal linking.
- `vertical.py`: 3D-to-depth interpolation and profile table helpers.
- `tilt.py`: tilt displacement and summary metrics.

Stage modules under `src/seacofs_eddy_dataset/stages/` should orchestrate I/O,
parallel execution, and resume behaviour. They should call `core/` functions
rather than accumulating science logic directly.

## Example Usage

```bash
python -m seacofs_eddy_dataset.cli run-stage detect_nencioli --config config/example.yaml
python -m seacofs_eddy_dataset.cli run-stage fit_doppio_surface --config config/example.yaml
python -m seacofs_eddy_dataset.cli run-all --config config/example.yaml
```

There is also a concise notebook control panel at
`notebooks/run_pipeline.ipynb`. It loads `config/local.yaml`, runs selected
pipeline stages, and shows quick output summaries.

`detect_nencioli` is implemented first. It writes one parquet partition per
`outer_avg_*.nc` file under `output_root/detections/`, using names like
`fnumber=01461.parquet`.

The remaining stages are scaffolded and will be implemented by moving logic out
of the existing notebooks in order.
