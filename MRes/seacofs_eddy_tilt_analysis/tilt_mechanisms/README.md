# Eddy tilt mechanism analysis

These notebooks test why the already-measured EAC eddies tilt. `TiltDis` and `TiltDir` are loaded as authoritative measurements and are never recomputed.

Run order:

1. `00_build_stratification_cache.ipynb`
2. `01_vertical_shear.ipynb`
3. `02_beta_stratification_burger.ipynb`
4. `03_pv_topographic_steering.ipynb`
5. `04_depth_dependent_propagation.ipynb`
6. `05_wind_ekman_sensitivity.ipynb` (optional; requires a wind-stress cache)
7. `06_joint_mechanism_comparison.ipynb`

The stratification cache is file-parallel and restartable. For every eddy-day,
it calculates depth-mean N2 at each ocean grid column inside the same
elliptical maximum-tangential-velocity contour used by `compute_core_mean`,
then takes the unweighted core mean. Centre-column values and core coverage
diagnostics are retained for sensitivity and quality control. It can be run from
the notebook or from Katana with:

```bash
python build_stratification_cache.py --workers 5 --point-batch-size 128
```

Start with 4–8 workers and check job memory before increasing concurrency.
The point batch size controls vectorized core-column reads and is not a Dask grid
chunk; the calculation never constructs full-domain xgcm metrics.
Each completed ROMS file is written to a cache-versioned partition folder;
reruns skip those partitions unless `--overwrite-partitions` is supplied.

`mechanism_tools.py` contains only helpers unique to this workflow. The notebooks reuse `seacofs_tilt_tools.py`, `ml_subsurface_tools.py`, and `beta_effect_background_flow/*` for existing functionality.

The notebooks assume they are launched with this directory as the working directory on Katana. Cache paths are declared near the top of the relevant notebooks.
