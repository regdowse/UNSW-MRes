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

## Example Usage

```bash
python -m seacofs_eddy_dataset.cli run-stage detect_nencioli --config config/example.yaml
python -m seacofs_eddy_dataset.cli run-stage fit_doppio_surface --config config/example.yaml
python -m seacofs_eddy_dataset.cli run-all --config config/example.yaml
```

The commands above are scaffolded first. The next step is to move the existing notebook logic into the corresponding stage modules.
