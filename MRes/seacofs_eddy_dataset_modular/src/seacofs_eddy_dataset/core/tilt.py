from __future__ import annotations

import numpy as np
import pandas as pd


def add_tilt_displacement(profile: pd.DataFrame) -> pd.DataFrame:
    out = profile.sort_values("Depth").copy()
    if out.empty:
        return out
    x0 = out.iloc[0]["xc"]
    y0 = out.iloc[0]["yc"]
    out["tilt_dx"] = out["xc"] - x0
    out["tilt_dy"] = out["yc"] - y0
    out["tilt_distance"] = np.hypot(out["tilt_dx"], out["tilt_dy"])
    out["tilt_angle_deg"] = np.degrees(np.arctan2(out["tilt_dy"], out["tilt_dx"]))
    return out


def summarise_tilt(profile: pd.DataFrame) -> dict:
    profile = add_tilt_displacement(profile)
    valid = profile.dropna(subset=["Depth", "tilt_distance"])
    if len(valid) < 2:
        return {"tilt_slope": np.nan, "max_tilt_distance": np.nan, "max_depth": np.nan}
    depth = valid["Depth"].to_numpy(dtype=float)
    distance = valid["tilt_distance"].to_numpy(dtype=float)
    slope = np.polyfit(depth, distance, deg=1)[0]
    return {
        "tilt_slope": slope,
        "max_tilt_distance": np.nanmax(distance),
        "max_depth": np.nanmax(depth),
    }

