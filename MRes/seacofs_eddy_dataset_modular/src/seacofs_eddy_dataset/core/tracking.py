from __future__ import annotations

import numpy as np
import pandas as pd


def clean_surface_eddies(df: pd.DataFrame, X_grid, Y_grid, threshold: float = 4.0, omega0_thresh: float = 5e-5):
    required = ["nCyc", "xc", "yc", "nxc", "nyc", "Omega0"]
    clean = df.dropna(subset=required).copy()

    nenc_cyc = np.where(clean["nCyc"].eq("AE"), 1, -1)
    doppio_cyc = np.sign(clean["Omega0"].values)
    clean = clean.loc[nenc_cyc == doppio_cyc].copy()

    err = np.hypot(clean["xc"] - clean["nxc"], clean["yc"] - clean["nyc"])
    clean["err"] = err
    boundary_mask = clean["xc"].between(0, X_grid.max(), inclusive="neither") & clean["yc"].between(
        0, Y_grid.max(), inclusive="neither"
    )

    log_err = np.log1p(err)
    err_med = np.median(log_err)
    err_mad = np.median(np.abs(log_err - err_med))
    upper_err = np.inf if err_mad == 0 else np.expm1(err_med + threshold * err_mad / 0.6745)

    clean["Omega0_abs"] = clean["Omega0"].abs()
    final_mask = (
        boundary_mask
        & (err <= upper_err)
        & (clean["Omega0_abs"] <= omega0_thresh)
        & np.isfinite(err)
        & np.isfinite(clean["Omega0_abs"])
        & (clean["Omega0_abs"] > 0)
    )
    return clean.loc[final_mask].reset_index(drop=True)

