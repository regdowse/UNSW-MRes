"""Shared utilities for the SEACOFS eddy tilt analysis notebooks.

The notebooks in this folder are intentionally thin: they define the question,
load the common data, then call the reusable helpers collected here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import linregress


DEFAULT_EDDY_PATH = Path(
    "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/"
    "DOPPIO_SEACOFS_26yr_50m_vert_check/df_eddies_50m_vert_checked_processed.pkl"
)
DEFAULT_TILT_PATH = Path(
    "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/"
    "DOPPIO_SEACOFS_26yr_50m_vert_check/df_tilt_vert_checked.pkl"
)
DEFAULT_VERT_PATH = Path(
    "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/"
    "DOPPIO_SEACOFS_26yr_50m_vert_check/dic_vert_doppio_50m_vert_checked.pkl"
)
DEFAULT_OLD_EDDY_PATH = Path(
    "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/df_eddies_processed.pkl"
)
DEFAULT_BETA_EDDY_PATH = Path(
    "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/df_eddies_beta_data_w.pkl"
)
DEFAULT_GRID_PATH = Path("/srv/scratch/z3533156/26year_BRAN2020/outer_avg_01461.nc")
DEFAULT_ZR_PATH = Path("/srv/scratch/z5297792/z_r.npy")
DEFAULT_CLIM_FUNCTIONS_DIR = Path("/home/z5297792/UNSW-MRes/MRes/SEACOFS_dataset")

KM_PER_DAY_TO_M_PER_S = 1000.0 / 86400.0
LEVELS_LAT = [-40, -35, -30, -25]
LEVELS_LON = [150, 155, 160]


@dataclass
class Paths:
    """Centralised paths used by the notebooks."""

    eddies: Path = DEFAULT_EDDY_PATH
    tilt: Path = DEFAULT_TILT_PATH
    vertical_dictionary: Path = DEFAULT_VERT_PATH
    old_eddies: Path = DEFAULT_OLD_EDDY_PATH
    beta_eddies: Path = DEFAULT_BETA_EDDY_PATH
    grid: Path = DEFAULT_GRID_PATH
    z_r: Path = DEFAULT_ZR_PATH
    clim_functions_dir: Path = DEFAULT_CLIM_FUNCTIONS_DIR


@dataclass
class Grid:
    lon_rho: np.ndarray
    lat_rho: np.ndarray
    mask_rho: np.ndarray
    h: np.ndarray
    f: np.ndarray
    angle: float
    z_r: np.ndarray
    x_grid: np.ndarray
    y_grid: np.ndarray
    X_grid: np.ndarray
    Y_grid: np.ndarray


def configure_clim_functions(path: Path | str = DEFAULT_CLIM_FUNCTIONS_DIR) -> None:
    """Add the existing SEACOFS helper folder to ``sys.path``."""

    path = str(path)
    if path not in sys.path:
        sys.path.append(path)


def clim_function(name: str):
    """Import one function from the existing ``clim_functions`` module."""

    configure_clim_functions()
    from clim_functions import __dict__ as clim_namespace

    return clim_namespace[name]


def distance_km(lat1, lon1, lat2, lon2, earth_radius_km: float = 6357.0):
    """Great-circle distance in kilometres."""

    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return earth_radius_km * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def load_grid(grid_path: Path | str = DEFAULT_GRID_PATH, z_r_path: Path | str = DEFAULT_ZR_PATH) -> Grid:
    """Load the SEACOFS grid and build kilometre x/y coordinates."""

    dataset = nc.Dataset(grid_path)
    lon_rho = np.transpose(dataset.variables["lon_rho"], axes=(1, 0))
    lat_rho = np.transpose(dataset.variables["lat_rho"], axes=(1, 0))
    mask_rho = np.transpose(dataset.variables["mask_rho"], axes=(1, 0))
    h = np.transpose(dataset.variables["h"], axes=(1, 0))
    f = np.transpose(dataset.variables["f"], axes=(1, 0))
    angle = float(dataset.variables["angle"][0, 0])
    z_r = np.transpose(np.load(z_r_path), (1, 2, 0))

    j_mid = lon_rho.shape[1] // 2
    i_mid = lon_rho.shape[0] // 2
    dx = distance_km(lat_rho[:-1, j_mid], lon_rho[:-1, j_mid], lat_rho[1:, j_mid], lon_rho[1:, j_mid])
    dy = distance_km(lat_rho[i_mid, :-1], lon_rho[i_mid, :-1], lat_rho[i_mid, 1:], lon_rho[i_mid, 1:])
    x_grid = np.insert(np.cumsum(dx), 0, 0)
    y_grid = np.insert(np.cumsum(dy), 0, 0)
    X_grid, Y_grid = np.meshgrid(x_grid, y_grid, indexing="ij")

    return Grid(lon_rho, lat_rho, mask_rho, h, f, angle, z_r, x_grid, y_grid, X_grid, Y_grid)


def load_tilt_tables(paths: Paths = Paths(), *, add_regions: bool = False, grid: Grid | None = None):
    """Load eddy and tilt tables, then merge ``TiltDis`` and ``TiltDir``."""

    df_eddies = pd.read_pickle(paths.eddies)
    df_tilt = pd.read_pickle(paths.tilt)
    df_eddies = df_eddies.merge(
        df_tilt[["Eddy", "Day", "TiltDis", "TiltDir"]],
        how="left",
        on=["Eddy", "Day"],
    )

    if add_regions:
        if grid is None:
            raise ValueError("Pass grid when add_regions=True.")
        df_eddies = add_region_labels(df_eddies, grid)

    return df_eddies, df_tilt


def load_vertical_dictionary(paths: Paths = Paths()):
    return pd.read_pickle(paths.vertical_dictionary)


def add_region_labels(df: pd.DataFrame, grid: Grid) -> pd.DataFrame:
    """Attach the six SEACOFS region labels used in the original notebooks."""

    add_region_column = clim_function("add_region_column")
    return add_region_column(df, grid.X_grid, grid.Y_grid, grid.lon_rho, grid.lat_rho, grid.h, grid.mask_rho)


def add_time_coordinates(df: pd.DataFrame) -> pd.DataFrame:
    """Add day index and normalised lifetime coordinates per eddy."""

    out = df.copy()
    out["Day_idx"] = out.groupby("Eddy").cumcount()
    max_day = out.groupby("Eddy")["Day_idx"].transform("max")
    out["norm_time"] = np.where(max_day > 0, out["Day_idx"] / max_day, np.nan)
    return out


def bearing_from_xy(x, y):
    """Bearing in degrees, where 0 is north and 90 is east."""

    return (90.0 - np.degrees(np.arctan2(y, x))) % 360.0


def angle_diff_180(a, b):
    """Absolute angular difference between two bearings in degrees."""

    return np.abs((a - b + 180.0) % 360.0 - 180.0)


def add_pv_gradient_terms(df: pd.DataFrame, grid: Grid) -> pd.DataFrame:
    """Compute planetary, topographic, and total shallow-water PV-gradient terms."""

    phys_grad = clim_function("phys_grad")
    out = df.copy()
    out["f"] = grid.f[out.ic, out.jc]
    out["h"] = grid.h[out.ic, out.jc]

    dhdx, dhdy = phys_grad(grid.h, grid.X_grid * 1e3, grid.Y_grid * 1e3, grid.mask_rho)
    dh_dN = -(np.sin(grid.angle) * dhdx + np.cos(grid.angle) * dhdy)
    dh_dE = -(np.cos(grid.angle) * dhdx - np.sin(grid.angle) * dhdy)
    out["dhdx"] = dh_dE[out.ic, out.jc]
    out["dhdy"] = dh_dN[out.ic, out.jc]

    dfdx, dfdy = phys_grad(grid.f, grid.X_grid * 1e3, grid.Y_grid * 1e3, grid.mask_rho)
    df_dN = -(np.sin(grid.angle) * dfdx + np.cos(grid.angle) * dfdy)
    out["beta"] = df_dN[out.ic, out.jc]

    omega_f = out["w"] + out["f"]
    out["PV_grad_plan_x"] = 0.0
    out["PV_grad_plan_y"] = out["beta"] / out["h"]
    out["PV_grad_topo_x"] = -omega_f * out["dhdx"] / out["h"] ** 2
    out["PV_grad_topo_y"] = -omega_f * out["dhdy"] / out["h"] ** 2
    out["PV_grad_x"] = out["PV_grad_plan_x"] + out["PV_grad_topo_x"]
    out["PV_grad_y"] = out["PV_grad_plan_y"] + out["PV_grad_topo_y"]

    for prefix in ["PV_grad_plan", "PV_grad_topo", "PV_grad"]:
        out[f"{prefix}_mag"] = np.hypot(out[f"{prefix}_x"], out[f"{prefix}_y"])
        out[f"{prefix}_theta"] = bearing_from_xy(out[f"{prefix}_x"], out[f"{prefix}_y"])

    out["dtheta_PV_grad"] = angle_diff_180(out["TiltDir"], out["PV_grad_theta"])
    out["dtheta_PV_grad_topo"] = angle_diff_180(out["TiltDir"], out["PV_grad_topo_theta"])
    out["dtheta_PV_grad_plan"] = angle_diff_180(out["TiltDir"], out["PV_grad_plan_theta"])
    out["Ro"] = np.abs(out["w"] / out["f"])
    out["topo_plan_ratio"] = np.log(out["PV_grad_topo_mag"] / out["PV_grad_plan_mag"])
    return out


def add_top_bottom_speeds(df: pd.DataFrame, dic_vert: dict, zmax: float = 1000.0) -> pd.DataFrame:
    """Add top-centre, bottom-centre, and surface-bottom propagation diagnostics."""

    out = df.copy()
    dt = out.groupby("Eddy")["Day"].diff()

    out["dx_top"] = out.groupby("Eddy")["xc"].diff()
    out["dy_top"] = out.groupby("Eddy")["yc"].diff()
    out["EddyProp"] = np.hypot(out.dx_top, out.dy_top) / dt * KM_PER_DAY_TO_M_PER_S

    x_btm, y_btm, z_btm = [], [], []
    for row in out.itertuples():
        try:
            profile = dic_vert[f"Eddy{int(row.Eddy)}"][f"Day{int(row.Day)}"]
        except KeyError:
            profile = None

        if profile is None or len(profile) == 0:
            x_btm.append(np.nan)
            y_btm.append(np.nan)
            z_btm.append(np.nan)
            continue

        deep = profile[profile.Depth.abs() <= zmax]
        bottom = deep.iloc[-1] if len(deep) else profile.iloc[-1]
        x_btm.append(bottom.xc)
        y_btm.append(bottom.yc)
        z_btm.append(bottom.Depth)

    out["x_btm"] = x_btm
    out["y_btm"] = y_btm
    out["z_btm"] = z_btm
    out["dx_btm"] = out.groupby("Eddy")["x_btm"].diff()
    out["dy_btm"] = out.groupby("Eddy")["y_btm"].diff()
    out["btm_prop"] = np.hypot(out.dx_btm, out.dy_btm) / dt * KM_PER_DAY_TO_M_PER_S
    out["sep_km"] = np.hypot(out.x_btm - out.xc, out.y_btm - out.yc)
    out["sep_rate_ms"] = out.groupby("Eddy")["sep_km"].diff() / dt * KM_PER_DAY_TO_M_PER_S
    out["top_btm_diff"] = np.hypot(out.dx_btm - out.dx_top, out.dy_btm - out.dy_top) / dt * KM_PER_DAY_TO_M_PER_S
    return out


def binned_tilt_panel(
    ax,
    df: pd.DataFrame,
    xcol: str,
    xlabel: str,
    *,
    ylabel: str | None = None,
    xlim: tuple[float, float] | None = None,
    percentile_xlim: tuple[float, float] = (10, 90),
    styles: dict | None = None,
    scatter: bool = False,
    linfit: bool = False,
    bins: int | None = None,
):
    """Median/IQR tilt-distance panel split by AE/CE."""

    if styles is None:
        styles = {
            "AE": {"line": "darkred", "fill": "red"},
            "CE": {"line": "navy", "fill": "blue"},
        }

    data = df.copy()
    data[xcol] = np.ma.filled(np.ma.asarray(data[xcol]), np.nan)
    data["TiltDis"] = np.ma.filled(np.ma.asarray(data["TiltDis"]), np.nan)
    finite_x = data[xcol].to_numpy(float)
    finite_x = finite_x[np.isfinite(finite_x)]
    if finite_x.size == 0:
        ax.set_axis_off()
        return ax

    if bins is None:
        bins = min(30, max(8, int(np.sqrt(finite_x.size))))
    edges = np.unique(np.nanquantile(finite_x, np.linspace(0, 1, bins + 1)))
    if len(edges) < 3:
        edges = np.linspace(np.nanmin(finite_x), np.nanmax(finite_x), 9)
    centers = 0.5 * (edges[:-1] + edges[1:])

    for cyc in ["AE", "CE"]:
        sub = data.loc[data.Cyc == cyc].dropna(subset=[xcol, "TiltDis"])
        x = sub[xcol].to_numpy(float)
        y = sub["TiltDis"].to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        idx = np.digitize(x, edges)
        med = np.array([np.nanmedian(y[idx == i]) if np.any(idx == i) else np.nan for i in range(1, len(edges))])
        q25 = np.array([np.nanquantile(y[idx == i], 0.25) if np.any(idx == i) else np.nan for i in range(1, len(edges))])
        q75 = np.array([np.nanquantile(y[idx == i], 0.75) if np.any(idx == i) else np.nan for i in range(1, len(edges))])

        if scatter:
            ax.scatter(x, y, s=1, alpha=0.08, color=styles[cyc]["fill"])
        mask = np.isfinite(med)
        ax.plot(centers[mask], med[mask], lw=2.5, color=styles[cyc]["line"], label=cyc)
        ax.fill_between(centers[mask], q25[mask], q75[mask], color=styles[cyc]["fill"], alpha=0.10)

        if linfit and ok.sum() > 3:
            lo, hi = np.nanpercentile(x, percentile_xlim)
            reg = (x >= lo) & (x <= hi)
            if reg.sum() > 2:
                m, c, *_ = linregress(x[reg], y[reg])
                xf = np.linspace(lo, hi, 200)
                ax.plot(xf, m * xf + c, "--", lw=2, color=styles[cyc]["fill"])

    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlim is None:
        ax.set_xlim(*np.nanpercentile(finite_x, percentile_xlim))
    else:
        ax.set_xlim(*xlim)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", alpha=0.2)
    return ax


def circular_mean_deg_true_north(deg):
    """Circular mean for true-north bearings in degrees."""

    x = np.asarray(deg, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    s = np.nanmean(np.sin(np.deg2rad(x)))
    c = np.nanmean(np.cos(np.deg2rad(x)))
    return np.rad2deg(np.arctan2(s, c)) % 360


def shared_bins(*arrays, min_bins: int = 12, max_bins: int = 40):
    """Freedman-Diaconis-like bins shared across one or more arrays."""

    vals = np.concatenate([np.asarray(a, float)[np.isfinite(a)] for a in arrays if len(np.asarray(a))])
    if vals.size == 0:
        return np.linspace(0, 1, min_bins + 1)
    q25, q75 = np.nanpercentile(vals, [25, 75])
    iqr = q75 - q25
    width = 2 * iqr / np.cbrt(vals.size) if iqr > 0 else (np.nanmax(vals) - np.nanmin(vals)) / min_bins
    if width <= 0 or not np.isfinite(width):
        width = 1.0
    n = int(np.clip(np.ceil((np.nanmax(vals) - np.nanmin(vals)) / width), min_bins, max_bins))
    return np.linspace(np.nanmin(vals), np.nanmax(vals), n + 1)


def mirrored_hist(ax, ae, ce, bins, xlabel, *, ylabel="Count", colors=("darkred", "navy")):
    """Plot AE above zero and CE below zero using shared bins."""

    ae_counts, edges = np.histogram(np.asarray(ae, float)[np.isfinite(ae)], bins=bins)
    ce_counts, _ = np.histogram(np.asarray(ce, float)[np.isfinite(ce)], bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    ax.bar(centers, ae_counts, width=widths, align="center", color=colors[0], alpha=0.35, label="AE")
    ax.bar(centers, -ce_counts, width=widths, align="center", color=colors[1], alpha=0.35, label="CE")
    ax.axhline(0, color="0.3", lw=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return ax


def choose_dir_bins_cardinal(*dfs, col: str = "TiltDir", min_bins: int = 8, max_bins: int = 36, min_avg_per_sector: int = 8):
    """Choose circular direction bins while keeping cardinal directions aligned."""

    candidates = np.array([4, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 30, 36, 45, 60, 72])
    candidates = candidates[(candidates >= min_bins) & (candidates <= max_bins)]
    candidates = candidates[candidates % 4 == 0]
    counts = []
    for df in dfs:
        x = df[col].to_numpy(float)
        counts.append(np.isfinite(x).sum())
    n = min(counts) if counts else 0
    if n == 0 or candidates.size == 0:
        k = 16
    else:
        k0 = int(np.ceil(2 * n ** (1 / 3)))
        k = candidates[np.argmin(np.abs(candidates - k0))]
        while k > candidates.min() and (n / k) < min_avg_per_sector:
            k = candidates[candidates < k].max()
    return np.linspace(0.0, 360.0, k + 1), 180.0 / k


def add_season(df: pd.DataFrame, time_col: str = "Date", season_col: str = "Season") -> pd.DataFrame:
    """Add austral seasons from a datetime-like column."""

    out = df.copy()
    month = pd.to_datetime(out[time_col]).dt.month
    out[season_col] = np.select(
        [month.isin([12, 1, 2]), month.isin([3, 4, 5]), month.isin([6, 7, 8]), month.isin([9, 10, 11])],
        ["DJF", "MAM", "JJA", "SON"],
        default=np.nan,
    )
    return out


def windrose_counts(directions_deg, magnitudes, *, mag_bins, dir_bins=None, dir_shift=None):
    """Count tilt directions and magnitudes for polar bar plots."""

    directions = np.asarray(directions_deg, float)
    magnitudes = np.asarray(magnitudes, float)
    ok = np.isfinite(directions) & np.isfinite(magnitudes)
    directions = directions[ok]
    magnitudes = magnitudes[ok]
    if directions.size == 0:
        return None
    if dir_bins is None or dir_shift is None:
        tmp = pd.DataFrame({"TiltDir": directions})
        dir_bins, dir_shift = choose_dir_bins_cardinal(tmp)
    k = len(dir_bins) - 1
    binw_deg = 360.0 / k
    angles = np.deg2rad(np.arange(k) * binw_deg)
    width = np.deg2rad(binw_deg)
    directions = np.mod(directions + dir_shift, 360.0)
    dir_idx = np.digitize(directions, dir_bins, right=False) - 1
    mag_idx = np.digitize(magnitudes, mag_bins, right=False) - 1
    counts = np.zeros((len(mag_bins) - 1, k), float)
    ok = (dir_idx >= 0) & (dir_idx < k) & (mag_idx >= 0) & (mag_idx < len(mag_bins) - 1)
    for d_i, m_i in zip(dir_idx[ok], mag_idx[ok]):
        counts[m_i, d_i] += 1
    return counts, angles, width


def plot_windrose(ax, df: pd.DataFrame, *, title: str = "", mag_bins=(0, 10, 20, 30, 40, np.inf), colors=None):
    """Draw one stacked tilt windrose."""

    if colors is None:
        cmap = "Reds" if df.Cyc.iloc[0] == "AE" else "Blues"
        colors = getattr(plt.cm, cmap)(np.linspace(0.15, 1, len(mag_bins) - 1))
    data = windrose_counts(df.TiltDir, df.TiltDis, mag_bins=mag_bins)
    if data is None:
        ax.set_axis_off()
        return ax
    counts, angles, width = data
    bottom = np.zeros(counts.shape[1])
    for i in range(counts.shape[0]):
        hi = "inf" if np.isinf(mag_bins[i + 1]) else f"{mag_bins[i + 1]:g}"
        ax.bar(angles, counts[i], width=width, bottom=bottom, color=colors[i], edgecolor=(0, 0, 0, 0.2), label=f"{mag_bins[i]:g}-{hi}")
        bottom += counts[i]
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_title(title)
    return ax


def assign_six_regions(df: pd.DataFrame, grid: Grid, *, lon_split: float = 157.0, lat_split: float = -33.5) -> tuple[pd.DataFrame, np.ndarray]:
    """Assign eddy days to the six map bins used by the windrose map."""

    region_mask = (grid.h < 4e3) & (grid.X_grid < 400) & (grid.lon_rho < 154.85) & (grid.mask_rho == 1)
    bin_grid = np.full(grid.X_grid.shape, np.nan)
    bin_grid[region_mask & (grid.lat_rho >= lat_split)] = 1
    bin_grid[region_mask & (grid.lat_rho < lat_split)] = 2
    bin_grid[(~region_mask) & (grid.lon_rho < lon_split) & (grid.mask_rho == 1) & (grid.lat_rho >= lat_split)] = 3
    bin_grid[(~region_mask) & (grid.lon_rho < lon_split) & (grid.mask_rho == 1) & (grid.lat_rho < lat_split)] = 4
    bin_grid[(grid.lon_rho >= lon_split) & (grid.mask_rho == 1) & (grid.lat_rho >= lat_split)] = 5
    bin_grid[(grid.lon_rho >= lon_split) & (grid.mask_rho == 1) & (grid.lat_rho < lat_split)] = 6

    out = df.copy()
    tree = cKDTree(np.column_stack([grid.X_grid.ravel(), grid.Y_grid.ravel()]))
    _, idx = tree.query(np.column_stack([out.xc, out.yc]))
    out["bin_id"] = bin_grid.ravel()[idx]
    out = out.dropna(subset=["bin_id"])
    out["bin_id"] = out["bin_id"].astype(int)
    return out, bin_grid


def point_from_bearing(origin, distance, bearing_deg):
    """Endpoint from an origin, distance, and true-north bearing."""

    theta_rad = np.radians(bearing_deg)
    dx = distance * np.sin(theta_rad)
    dy = distance * np.cos(theta_rad)
    return origin[0] - dx, origin[1] - dy


def match_old_eddies(sample_eddies_old, df_eddies_old, df_eddies, min_overlap_frac: float = 0.5, max_mean_dist: float = np.inf):
    """Match old sample eddy IDs to the vertically checked eddy IDs."""

    matches = []
    for eddy_old in sample_eddies_old:
        old = (
            df_eddies_old.loc[df_eddies_old.Eddy == eddy_old, ["Day", "xc", "yc"]]
            .drop_duplicates("Day")
            .sort_values("Day")
        )
        old_days = set(old.Day)
        candidates = df_eddies.loc[df_eddies.Day.isin(old_days), "Eddy"].unique()
        scores = []
        for eddy_new in candidates:
            new = (
                df_eddies.loc[df_eddies.Eddy == eddy_new, ["Day", "xc", "yc"]]
                .drop_duplicates("Day")
                .sort_values("Day")
            )
            merged = old.merge(new, on="Day", suffixes=("_old", "_new"))
            if merged.empty:
                continue
            overlap_frac = len(merged) / len(old)
            mean_dist = np.hypot(merged.xc_old - merged.xc_new, merged.yc_old - merged.yc_new).mean()
            if overlap_frac >= min_overlap_frac and mean_dist <= max_mean_dist:
                scores.append((eddy_new, overlap_frac, mean_dist, len(merged)))
        if scores:
            best = sorted(scores, key=lambda x: (-x[1], x[2]))[0]
            matches.append({"old_eddy": eddy_old, "new_eddy": best[0], "overlap_frac": best[1], "mean_dist_km": best[2], "n_overlap": best[3]})
        else:
            matches.append({"old_eddy": eddy_old, "new_eddy": np.nan, "overlap_frac": 0.0, "mean_dist_km": np.nan, "n_overlap": 0})
    return pd.DataFrame(matches)
