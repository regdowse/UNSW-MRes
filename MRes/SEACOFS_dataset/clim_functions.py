import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import time
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import Normalize
from matplotlib.ticker import FormatStrFormatter

# import sys
# sys.path.append('/home/z5297792/UNSW-MRes/MRes/SEACOFS_dataset') 

import sys
sys.path.append("/home/z5297792/ESP_zonodo")
from functions import doppio, out_core_param_fit

fname = f'/srv/scratch/z3533156/26year_BRAN2020/outer_avg_01461.nc'
dataset = nc.Dataset(fname)
lon_rho = np.transpose(dataset.variables['lon_rho'], axes=(1, 0))
lat_rho = np.transpose(dataset.variables['lat_rho'], axes=(1, 0))
mask_rho = np.transpose(dataset.variables['mask_rho'], axes=(1, 0))
h = np.transpose(dataset.variables['h'], axes=(1, 0))
f = np.transpose(dataset.variables['f'], axes=(1, 0))
angle = dataset.variables['angle'][0, 0]
z_r = np.load('/srv/scratch/z5297792/z_r.npy')
z_r = np.transpose(z_r, (1, 2, 0))
def distance(lat1, lon1, lat2, lon2):
    EARTH_RADIUS = 6357
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    return EARTH_RADIUS * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
j_mid, i_mid = lon_rho.shape[1] // 2, lon_rho.shape[0] // 2
dx = distance(lat_rho[:-1, j_mid], lon_rho[:-1, j_mid],
              lat_rho[1:, j_mid], lon_rho[1:, j_mid])
dy = distance(lat_rho[i_mid, :-1], lon_rho[i_mid, :-1],
              lat_rho[i_mid, 1:], lon_rho[i_mid, 1:])
x_grid = np.insert(np.cumsum(dx), 0, 0)
y_grid = np.insert(np.cumsum(dy), 0, 0)
X_grid, Y_grid = np.meshgrid(x_grid, y_grid, indexing='ij')



def nencioli(u, v, X, Y, a, b, flip_sign=False):
    """
    Nencioli eddy detection on a Cartesian grid (X, Y in km)

    Returns:
    - eddy_uv: points passing constraints 1–2
    - eddy_c: points passing constraints 1–3
    - eddy: final eddy centres (x, y, type)
    
    type (cyclonic=-1, anticyclonic=1) (personal preference Reg Dowse)
    """

    borders = max(a, b) + 1  # ensures all stencil operations stay inside domain
    vel = np.sqrt(u**2 + v**2)
    bound = vel.shape

    eddy_uv = np.empty((0, 2))
    eddy_c = np.empty((0, 2))
    eddy = np.empty((0, 3))

    for i in range(borders, bound[0] - borders + 1):

        wrk = v[i, :]  # zonal slice of v (detect zero-crossings along X)
        s = np.sign(wrk)
        ds = np.diff(s)

        # --- CONSTRAINT 1: zero-crossing in v + monotonicity away from crossing ---
        indx = np.where((ds != 0) & ~np.isnan(ds))[0]  # valid sign changes only
        indx = indx[(indx >= borders) & (indx <= bound[1] - borders - 1)]

        for ii in indx:
            var = 0  # 1 = cyclonic, -1 = anticyclonic

            if wrk[ii] >= 0:  # candidate anticyclone
                if wrk[ii - a] > wrk[ii] and wrk[ii + 1 + a] < wrk[ii + 1]:
                    var = -1  # satisfies constraint 1
            else:  # candidate cyclone
                if wrk[ii - a] < wrk[ii] and wrk[ii + 1 + a] > wrk[ii + 1]:
                    var = 1

            # --- CONSTRAINT 2: u reversal across the zero-crossing (ii OR ii+1) ---
            if var == -1:
                ok1 = (u[i - a, ii] <= 0 and u[i - a, ii] <= u[i - 1, ii] and
                       u[i + a, ii] >= 0 and u[i + a, ii] >= u[i + 1, ii])
                ok2 = (u[i - a, ii + 1] <= 0 and u[i - a, ii + 1] <= u[i - 1, ii + 1] and
                       u[i + a, ii + 1] >= 0 and u[i + a, ii + 1] >= u[i + 1, ii + 1])

                if ok1 or ok2:  # must satisfy reversal at at least one side
                    eddy_uv = np.vstack([
                        eddy_uv,
                        [X[i, ii], Y[i, ii]],
                        [X[i, ii + 1], Y[i, ii + 1]]
                    ])
                else:
                    var = 0  # reject if no reversal

            elif var == 1:
                ok1 = (u[i - a, ii] >= 0 and u[i - a, ii] >= u[i - 1, ii] and
                       u[i + a, ii] <= 0 and u[i + a, ii] <= u[i + 1, ii])
                ok2 = (u[i - a, ii + 1] >= 0 and u[i - a, ii + 1] >= u[i - 1, ii + 1] and
                       u[i + a, ii + 1] <= 0 and u[i + a, ii + 1] <= u[i + 1, ii + 1])

                if ok1 or ok2:
                    eddy_uv = np.vstack([
                        eddy_uv,
                        [X[i, ii], Y[i, ii]],
                        [X[i, ii + 1], Y[i, ii + 1]]
                    ])
                else:
                    var = 0

            # --- CONSTRAINT 3: local minimum of velocity magnitude ---
            if var != 0:
                srch = vel[i - b:i + b + 1, ii - b:ii + b + 2]  # MATLAB-consistent window
                sx = X[i - b:i + b + 1, ii - b:ii + b + 2]
                sy = Y[i - b:i + b + 1, ii - b:ii + b + 2]

                Xind, Yind = np.unravel_index(np.nanargmin(srch), srch.shape)
                imin = i - b + Xind
                jmin = ii - b + Yind

                # second check: ensure this minimum is locally consistent
                srch2 = vel[
                    max(imin - b, 0):min(imin + b + 1, bound[0]),
                    max(jmin - b, 0):min(jmin + b + 1, bound[1])
                ]

                if np.nanmin(srch2) == np.nanmin(srch):
                    eddy_c = np.vstack([eddy_c, [sx[Xind, Yind], sy[Xind, Yind]]])
                else:
                    var = 0  # reject if not a true local minimum

            # --- CONSTRAINT 4: consistent rotation of velocity vectors around centre ---
            if var != 0:
                d = a - 1
                i1, i2 = imin, jmin

                u_small = u[
                    max(i1 - d, 0):min(i1 + d + 1, bound[0]),
                    max(i2 - d, 0):min(i2 + d + 1, bound[1])
                ]
                v_small = v[
                    max(i1 - d, 0):min(i1 + d + 1, bound[0]),
                    max(i2 - d, 0):min(i2 + d + 1, bound[1])
                ]

                if not np.isnan(u_small).any() and not np.isnan(v_small).any():  # require full stencil
                    # extract boundary vectors in anticlockwise order
                    u_bound = np.concatenate([
                        u_small[0, :],
                        u_small[1:, -1],
                        u_small[-1, -2::-1],
                        u_small[-2:0:-1, 0]
                    ])
                    v_bound = np.concatenate([
                        v_small[0, :],
                        v_small[1:, -1],
                        v_small[-1, -2::-1],
                        v_small[-2:0:-1, 0]
                    ])

                    # assign each boundary vector to a quadrant (1→4)
                    quadrants = np.zeros_like(u_bound, dtype=int)
                    quadrants[(u_bound >= 0) & (v_bound >= 0)] = 1
                    quadrants[(u_bound < 0) & (v_bound >= 0)] = 2
                    quadrants[(u_bound < 0) & (v_bound < 0)] = 3
                    quadrants[(u_bound >= 0) & (v_bound < 0)] = 4

                    spin = np.where(quadrants == 4)[0]

                    # require full rotation and avoid trivial cases
                    if spin.size > 0 and spin.size != quadrants.size:
                        if spin[0] == 0:
                            spin = np.where(quadrants != 4)[0]
                            spin = np.array([spin[0] - 1])

                        quadrants[spin[-1] + 1:] += 4  # unwrap sequence to enforce monotonicity

                        # check no jumps >1 quadrant and no backward rotation
                        if not np.any(np.diff(quadrants) > 1) and not np.any(np.diff(quadrants) < 0):
                            eddy = np.vstack([eddy, [sx[Xind, Yind], sy[Xind, Yind], var]])

    eddy_uv = np.unique(eddy_uv, axis=0)
    eddy_c = np.unique(eddy_c, axis=0)
    eddy = np.unique(eddy, axis=0)

    if flip_sign and eddy.size:
        eddy[:, 2] *= -1  # optional convention flip only (no hemisphere logic in Cartesian)

    return eddy_uv, eddy_c, eddy


def doppio_pipeliner(ic, jc, ut, vt, X, Y, r=30_000.0):
    '''
    Return orthogonal transects centered at grid index (ic, jc) with radius r.
    '''
    x = np.asarray(X[:, 0], float)
    y = np.asarray(Y[0, :], float)
    nxc = x[ic]
    nyc = y[jc]

    # x-transect: y = y[jc]
    i0 = np.searchsorted(x, nxc - r, side="left")
    i1 = np.searchsorted(x, nxc + r, side="right")

    x1 = x[i0:i1]
    y1 = np.full(x1.size, nyc)
    u1 = ut[i0:i1, jc]
    v1 = vt[i0:i1, jc]

    # y-transect: x = x[ic]
    j0 = np.searchsorted(y, nyc - r, side="left")
    j1 = np.searchsorted(y, nyc + r, side="right")

    y2 = y[j0:j1]
    x2 = np.full(y2.size, nxc)
    u2 = ut[ic, j0:j1]
    v2 = vt[ic, j0:j1]

    return x1, y1, u1, v1, x2, y2, u2, v2

def tangential_velocity(xp, yp, up, vp, xc, yc, Q, det1=False):
    Q = np.asarray(Q, float)
    if Q.shape == (3,):
        q11, q12, q22 = Q
        Q = np.array([[q11, q12], [q12, q22]], float)
    if det1:
        d = np.linalg.det(Q)
        if d != 0:
            Q /= np.sqrt(d)

    xp, yp, up, vp = (np.asarray(a, float) for a in (xp, yp, up, vp))
    r   = np.stack((xp - xc, yp - yc), axis=-1)
    g   = 2.0 * (r @ Q.T)                    # ∇F
    J   = np.array([[0., -1.], [1., 0.]])    # +90° rot
    tau = g @ J.T                             # tangent
    nrm = np.linalg.norm(tau, axis=-1, keepdims=True)
    t_hat = np.divide(tau, nrm, out=np.zeros_like(tau), where=nrm > 0)

    vel = np.stack((up, vp), axis=-1)
    vt  = np.sum(vel * t_hat, axis=-1)
    vt  = np.where(nrm.squeeze() > 0, vt, np.nan)
    return vt

def find_directional_radii(u, v, x, y, xc, yc, Q, return_index=False):
    """
    Returns dict of 'up','right','down','left' where each value is either:
      - steps from (nic,njc) where |v_theta| stops growing (return_index=True), or
      - Euclidean distance from (xc,yc) to that stopping point (return_index=False).
    Stops early if a NaN is encountered.
    """
    dis = np.hypot(x - xc, y - yc)
    dis_f = np.where(np.isfinite(dis), dis, np.inf)
    nic, njc = np.unravel_index(np.argmin(dis_f), dis.shape)

    def walk(di, dj, max_r):
        v_old = 0.0
        steps = 0
        for r in range(1, max_r + 1):
            i, j = nic + di * r, njc + dj * r
            vt = np.abs(tangential_velocity(x[i, j], y[i, j], u[i, j], v[i, j], xc, yc, Q))
            if np.isnan(vt) or vt < v_old:
                break
            v_old = vt
            steps = r
        return steps

    steps = {
        'up':    walk(-1,  0, nic),
        'right': walk( 0,  1, u.shape[1] - njc - 1),
        'down':  walk( 1,  0, u.shape[0] - nic - 1),
        'left':  walk( 0, -1, njc),
    }

    if return_index:
        return steps

    dists = {}
    for direction, r in steps.items():
        if r == 0:
            dists[direction] = 0.0
        else:
            if direction == 'up':
                i0, j0 = nic - r, njc
            elif direction == 'right':
                i0, j0 = nic, njc + r
            elif direction == 'down':
                i0, j0 = nic + r, njc
            else:
                i0, j0 = nic, njc - r
            dists[direction] = float(np.hypot(x[i0, j0] - xc, y[i0, j0] - yc))
    return dists


def compute_AR_from_Q(Q):
    q11 = Q[:, 0, 0]
    q12 = Q[:, 1, 0]
    q22 = Q[:, 1, 1]

    tr  = q11 + q22
    rad = np.sqrt((q11 - q22)**2 + 4*q12**2)

    lam_min = 0.5 * (tr - rad)
    lam_max = 0.5 * (tr + rad)

    AR = np.sqrt(lam_max / np.maximum(lam_min, 1e-12))
    AR[lam_min <= 0] = np.nan

    return AR


def day_plot(day, df_eddies, out_core_flag=False):

    fnumber = 1461 + ((day - 1462) // 30)*30
    fname = f'/srv/scratch/z3533156/26year_BRAN2020/outer_avg_{fnumber:05}.nc'
    with nc.Dataset(fname) as ds:
        ocean_time = ds['ocean_time'][:] / 86400
        t = np.where(ocean_time == day)[0][0]
        ut = ds['u_eastward'][t, -1, :, :].T
        vt = ds['v_northward'][t, -1, :, :].T

    df_day = df_eddies.loc[df_eddies.Day.eq(day)].copy()

    cs = np.hypot(ut, vt)

    fig, ax = plt.subplots(figsize=(8, 10))
    im = ax.pcolor(X_grid, Y_grid, cs, shading='nearest', vmin=0, vmax=2.5, cmap='Blues_r')
    fig.colorbar(im, ax=ax, label=r'Current speed (ms$^{-1}$)')

    clrs = np.where(df_day.Cyc.eq('CE'), 'c', 'r')
    ax.scatter(df_day.xc, df_day.yc, c=clrs, edgecolors='k', linewidths=0.8, s=60, zorder=10)

    if 'Q' not in df_day.columns:
        df_day['Q'] = list(
            np.stack([
                np.stack([df_day.q11.values, df_day.q12.values], axis=1),
                np.stack([df_day.q12.values, df_day.q22.values], axis=1)
            ], axis=1)
        )

    for xc, yc, e, Q, Rc, R, cyc in zip(df_day.xc, df_day.yc, df_day.Eddy, df_day.Q, df_day.Rc, df_day.R, df_day.Cyc):

        # ----- Where I plot the eddy's maximum tangenital velocity contour -----
        dx_ell, dy_ell = X_grid - xc, Y_grid - yc
        rho2_ell = Q[0,0]*dx_ell**2 + 2*Q[1,0]*dx_ell*dy_ell + Q[1,1]*dy_ell**2 # rho^2
        ax.contour(X_grid, Y_grid, rho2_ell, levels=[Rc**2/2], colors='r' if cyc=='AE' else 'c')
        if out_core_flag:
            ax.contour(X_grid, Y_grid, rho2_ell, levels=[(1.75*R)**2], linestyles='--', colors='r' if cyc=='AE' else 'c')
        # -----------------------------------------------------------------------

        ax.annotate(
            str(e), (xc, yc),
            textcoords='offset points', xytext=(3, 3),
            fontsize=12, color='w', weight='bold',
            path_effects=[pe.withStroke(linewidth=2, foreground='k')],
            zorder=11
        )

    c1 = ax.contour(X_grid, Y_grid, lat_rho, levels=[-40, -35, -30, -25], colors='k', linewidths=.5)
    ax.clabel(c1, fmt=lambda v: f"{np.abs(v):.0f}°S", inline=True, colors='k')
    c2 = ax.contour(X_grid, Y_grid, lon_rho, levels=[150, 155, 160], colors='k', linewidths=.5)
    ax.clabel(c2, fmt=lambda v: f"{v:.0f}°E", inline=True, colors='k')
                
    ax.set_title(f'Day {day} | {pd.Timestamp("1990-01-01") + pd.Timedelta(days=day)}')
    ax.set_aspect('equal', adjustable='datalim')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_xlim(x_grid.min(), x_grid.max())
    ax.set_ylim(y_grid.min(), y_grid.max())


# 50m depth check parallisation

def rotate_uv(u, v, angle):
    u = np.where(np.abs(u) > 1e30, np.nan, u).astype(float)
    v = np.where(np.abs(v) > 1e30, np.nan, v).astype(float)
    u_rot = v * np.sin(angle) + u * np.cos(angle)
    v_rot = v * np.cos(angle) - u * np.sin(angle)

    return u_rot, v_rot

def transect_indexer(ic, jc, X, Y, r=30.0):
    x = np.asarray(X[:, 0])
    y = np.asarray(Y[0, :])

    x0, y0 = x[ic], y[jc]

    i0 = np.searchsorted(x, x0 - r, side='left')
    i1 = np.searchsorted(x, x0 + r, side='right')
    j0 = np.searchsorted(y, y0 - r, side='left')
    j1 = np.searchsorted(y, y0 + r, side='right')

    ii = np.arange(i0, i1)
    jj = np.arange(j0, j1)

    return x[ii], np.full(ii.size, y0), np.full(jj.size, x0), y[jj], ii, jj

def nearest_ij(xc, yc, X, Y):
    x = np.asarray(X[:, 0])
    y = np.asarray(Y[0, :])

    ic = np.abs(x - xc).argmin()
    jc = np.abs(y - yc).argmin()

    return int(ic), int(jc)

def interp_3d_to_reference_depths(var3d, z3d, target_depths):
    target_depths = np.asarray(target_depths)
    nz = target_depths.size

    out = np.full(var3d.shape[:2] + (nz,), np.nan)

    for i in range(var3d.shape[0]):
        for j in range(var3d.shape[1]):
            out[i, j, :] = np.interp(
                target_depths,
                z3d[i, j, :],
                var3d[i, j, :],
                left=np.nan,
                right=np.nan
            )

    return out 

def local_eddy_arrays(X, Y, u2d, v2d, xc, yc, Q, rho_max):
    local = (
        (X >= xc - rho_max)
        & (X <= xc + rho_max)
        & (Y >= yc - rho_max)
        & (Y <= yc + rho_max)
    )

    Xloc = X[local]
    Yloc = Y[local]
    uloc = u2d[local]
    vloc = v2d[local]

    dx = Xloc - xc
    dy = Yloc - yc

    rho2 = (
        Q[0, 0] * dx**2
        + 2 * Q[0, 1] * dx * dy
        + Q[1, 1] * dy**2
    )

    rho = np.sqrt(np.where(rho2 >= 0, rho2, np.nan))

    return Xloc, Yloc, uloc, vloc, rho

def add_reaches_50m_flag_fname(
    df,
    fname,
    X_grid, Y_grid, angle, z_r, target_depths,
    r=30.0, max_jump=100,
    depth_thresh=50.0,
    flag_col='reaches_50m'
):

    z_r = np.abs(z_r)

    if flag_col not in df.columns:
        df[flag_col] = np.nan

    if 'ic' not in df.columns or 'jc' not in df.columns:
        ij = np.array([
            nearest_ij(xc, yc, X_grid, Y_grid)
            for xc, yc in zip(df.xc, df.yc)
        ])
        df['ic'], df['jc'] = ij[:, 0], ij[:, 1]

    df_file = df.loc[df.fname.eq(fname)]

    if df_file.empty:
        return df

    with nc.Dataset(fname) as ds:

        ocean_time = ds['ocean_time'][:] / 86400
        time_lookup = {int(d): i for i, d in enumerate(ocean_time)}

        for day, df_day in df_file.groupby('Day', sort=False):

            t = time_lookup.get(int(day))

            if t is None:
                continue

            u3d = np.flip(
                ds['u_eastward'][t].T.astype(float),
                axis=2
            )

            v3d = np.flip(
                ds['v_northward'][t].T.astype(float),
                axis=2
            )

            u3d[np.abs(u3d) > 1e30] = np.nan
            v3d[np.abs(v3d) > 1e30] = np.nan

            u_depth = interp_3d_to_reference_depths(
                u3d, z_r, target_depths
            )

            v_depth = interp_3d_to_reference_depths(
                v3d, z_r, target_depths
            )

            for row in df_day.itertuples():

                xc_prev, yc_prev = row.xc, row.yc
                ic, jc = int(row.ic), int(row.jc)
                w_surf = row.w

                reached = False

                for k, target_depth in enumerate(target_depths):

                    if (
                        (xc_prev < r)
                        or (xc_prev > X_grid.max() - r)
                        or (yc_prev < r)
                        or (yc_prev > Y_grid.max() - r)
                    ):
                        # print(1)
                        break

                    u2d = u_depth[:, :, k]
                    v2d = v_depth[:, :, k]

                    x1, y1, x2, y2, ii, jj = transect_indexer(
                        ic, jc, X_grid, Y_grid, r=r
                    )

                    u1, v1 = rotate_uv(
                        u2d[ii, jc],
                        v2d[ii, jc],
                        angle
                    )

                    u2, v2 = rotate_uv(
                        u2d[ic, jj],
                        v2d[ic, jj],
                        angle
                    )

                    if any(
                        np.all(np.isnan(a))
                        for a in [u1, v1, u2, v2]
                    ):
                        # print(2)
                        break

                    try:
                        xc, yc, w, Q, Omega0 = doppio(
                            x1, y1, u1, v1,
                            x2, y2, u2, v2
                        )
                    except Exception:
                        # print(3)
                        break

                    if (
                        np.sign(w) != np.sign(w_surf)
                        or np.hypot(
                            xc - xc_prev,
                            yc - yc_prev
                        ) > max_jump
                    ):
                        # print(4)
                        break

                    if target_depth >= depth_thresh:
                        reached = True
                        break

                    xc_prev = xc
                    yc_prev = yc

                df.loc[row.Index, flag_col] = reached

    return df

# vert dataset
def vert_doppio_dataset(
    dic, df_eddies, fname, X_grid, Y_grid, angle, z_r, 
    r=30.0,
    max_jump=100,
    max_depth_levels=None,
    rho_max=200,
    rho_min=30,
    target_depths=None
):
    z_r = np.abs(z_r)

    df_file = df_eddies.loc[df_eddies.fname.eq(fname)]

    if df_file.empty:
        return dic

    with nc.Dataset(fname) as ds:

        ocean_time = ds['ocean_time'][:] / 86400
        time_lookup = {int(d): i for i, d in enumerate(ocean_time)}

        for day, df_day in df_file.groupby('Day', sort=False):

            t = time_lookup.get(int(day))
            if t is None:
                continue

            u3d = ds['u_eastward'][t, :, :, :].T.astype(float)
            u3d = np.flip(u3d, axis=2)
            v3d = ds['v_northward'][t, :, :, :].T.astype(float)
            v3d = np.flip(v3d, axis=2)

            u3d[np.abs(u3d) > 1e30] = np.nan
            v3d[np.abs(v3d) > 1e30] = np.nan

            # nz = u3d.shape[-1]
            nz = len(target_depths)
            if max_depth_levels is not None:
                nz = min(nz, max_depth_levels)

            u_depth = interp_3d_to_reference_depths(
                u3d, z_r, target_depths
            )
            
            v_depth = interp_3d_to_reference_depths(
                v3d, z_r, target_depths
            )

            for row in df_day.itertuples():

                eddy_key = f'Eddy{row.Eddy}'
                day_key = f'Day{int(row.Day)}'

                dic.setdefault(eddy_key, {})

                out = []

                out.append({
                    'z': 0,
                    'Depth': 0,
                    'xc': row.xc,
                    'yc': row.yc,
                    'ic': row.ic,
                    'jc': row.jc,
                    'w': row.w,
                    'Q': np.array([[row.q11, row.q12], [row.q12, row.q22]]),
                    'Omega0': np.nan,
                    'Omega': row.Omega,
                    'Rc': row.Rc,
                    'psi0': row.psi0,
                    'R': row.R
                })

                xc_prev = row.xc
                yc_prev = row.yc
                w_surf = row.w

                ic = int(row.ic)
                jc = int(row.jc)
                # print(eddy_key)
                for k in range(nz):
                    # print(f'k = {k}')
            
                    if (
                        (xc_prev < r)
                        or (xc_prev > X_grid.max() - r)
                        or (yc_prev < r)
                        or (yc_prev > Y_grid.max() - r)
                    ):
                        # print(1)
                        break

                    target_depth = target_depths[k]

                    u2d = u_depth[:, :, k]
                    v2d = v_depth[:, :, k]

                    x1, y1, x2, y2, ii, jj = transect_indexer(
                        ic, jc, X_grid, Y_grid, r=r
                    )

                    u1 = u2d[ii, jc]
                    v1 = v2d[ii, jc]
                    u2 = u2d[ic, jj]
                    v2 = v2d[ic, jj]

                    u1, v1 = rotate_uv(u1, v1, angle)
                    u2, v2 = rotate_uv(u2, v2, angle)

                    if any(np.all(np.isnan(a)) for a in [u1, v1, u2, v2]):
                        # print(2)
                        break

                    try:
                        xc, yc, w, Q, Omega0 = doppio(
                            x1, y1, u1, v1,
                            x2, y2, u2, v2
                        )
                    except Exception:
                        print(3)
                        break

                    if (
                        np.sign(w) != np.sign(w_surf)
                        or np.hypot(xc - xc_prev, yc - yc_prev) > max_jump
                    ):
                        # print(4)
                        break

                    ic, jc = nearest_ij(xc, yc, X_grid, Y_grid)

                    u2d_rot, v2d_rot = rotate_uv(u2d, v2d, angle)

                    Xloc, Yloc, uloc, vloc, rho = local_eddy_arrays(
                        X_grid, Y_grid,
                        u2d_rot, v2d_rot,
                        xc, yc, Q,
                        rho_max=rho_max
                    )

                    mask0 = rho < rho_max

                    if mask0.sum() < 10:
                        break

                    radii = find_directional_radii(
                        u2d_rot, v2d_rot, X_grid, Y_grid, xc, yc, Q
                    )

                    R = np.nanmean([
                        radii['up'],
                        radii['right'],
                        radii['down'],
                        radii['left'],
                    ])

                    if not np.isfinite(R):
                        # print(5)
                        break

                    rho_lim = max(min(R * 1.75, rho_max), rho_min)
                    mask = rho < rho_lim

                    if mask.sum() < 10:
                        break

                    try:
                        Rc, psi0, Omega = out_core_param_fit(
                            Xloc[mask],
                            Yloc[mask],
                            uloc[mask],
                            vloc[mask],
                            xc, yc, Q,
                            Omega0=Omega0
                        )
                    except Exception:
                        break

                    out.append({
                        'z': k+1,
                        'Depth': target_depth,
                        'xc': xc,
                        'yc': yc,
                        'ic': ic,
                        'jc': jc,
                        'w': w*1e-3,
                        'Q': Q,
                        'Omega0': Omega0*1e-3,
                        'Omega': Omega*1e-3,
                        'Rc': Rc,
                        'psi0': psi0,
                        'R': R
                    })

                    xc_prev = xc
                    yc_prev = yc

                dic[eddy_key][day_key] = pd.DataFrame(out)

    return dic



# Tracking

def collect_tracking_R_with_Omega0(
    df_data,
    L_SCALE=50,
    Omega0_SCALE=1e-5,
    LOOKBACK=4,
    K_NEIGH=5
):
    df = df_data.copy()
    unique_days = np.sort(df['Day'].unique())

    daily_groups = {
        d: (
            df.loc[(df['Day'] == d) & df['xc'].notna()]
              .groupby('eddy_idx', as_index=False, sort=False)
              .first()
        )
        for d in unique_days
    }

    rows = []

    for day in unique_days[1:]:

        pres = daily_groups.get(day)
        if pres is None or len(pres) == 0:
            continue

        for _, pres_eddy in pres.iterrows():

            for delta in range(1, LOOKBACK + 1):

                prev_day = day - delta
                prev = daily_groups.get(prev_day)

                if prev is None or len(prev) == 0:
                    continue

                prev = prev[prev['Cyc'].eq(pres_eddy['Cyc'])]

                if len(prev) == 0:
                    continue

                coords = np.column_stack([
                    prev['xc'].values / L_SCALE,
                    prev['yc'].values / L_SCALE,
                    prev['Omega0'].values / Omega0_SCALE
                ])

                tree = cKDTree(coords)

                query = np.array([
                    pres_eddy['xc'] / L_SCALE,
                    pres_eddy['yc'] / L_SCALE,
                    pres_eddy['Omega0'] / Omega0_SCALE
                ])

                dist, idx = tree.query(query, k=min(K_NEIGH, len(prev)))

                dist = np.atleast_1d(dist)
                idx = np.atleast_1d(idx)

                for dR, j in zip(dist, idx):

                    prev_eddy = prev.iloc[j]

                    rows.append({
                        'Day': day,
                        'prev_day': prev_day,
                        'delta': delta,
                        'eddy_idx': pres_eddy['eddy_idx'],
                        'prev_eddy_idx': prev_eddy['eddy_idx'],
                        'R': dR,
                        'dx': pres_eddy['xc'] - prev_eddy['xc'],
                        'dy': pres_eddy['yc'] - prev_eddy['yc'],
                        'dO': pres_eddy['Omega0'] - prev_eddy['Omega0'],
                        'Cyc': pres_eddy['Cyc'],
                    })

                # only diagnose most recent available day
                break

    return pd.DataFrame(rows)


def tracking_kdtree_with_Omega0(
    df_data,
    start_ID,
    next_num,
    L_SCALE=50,
    Omega0_SCALE=1e-5,
    R_THRESH=1,
    LOOKBACK=4
):
    tic = time.perf_counter()
    df = df_data.dropna(subset=['xc', 'yc', 'Omega0']).copy()

    min_day = df['Day'].min()
    df['Eddy'] = -1
    df.loc[df['Day'] == min_day, 'Eddy'] = start_ID
    df['Eddy'] = df['Eddy'].astype('Int64')

    unique_days = np.sort(df['Day'].unique())

    for day in unique_days[1:]:

        pres = (
            df.loc[df['Day'].eq(day)]
              .groupby('eddy_idx', as_index=False, sort=False)
              .first()
        )

        assigned = set()

        for _, pres_eddy in pres.iterrows():

            best_id = None

            query = np.array([
                pres_eddy['xc'] / L_SCALE,
                pres_eddy['yc'] / L_SCALE,
                pres_eddy['Omega0'] / Omega0_SCALE
            ])

            for delta in range(1, LOOKBACK + 1):

                prev_day = day - delta

                prev = (
                    df.loc[df['Day'].eq(prev_day) & df['Eddy'].ne(-1)]
                      .groupby('eddy_idx', as_index=False, sort=False)
                      .first()
                )

                if len(prev) == 0:
                    continue

                coords = np.column_stack([
                    prev['xc'].values / L_SCALE,
                    prev['yc'].values / L_SCALE,
                    prev['Omega0'].values / Omega0_SCALE
                ])

                tree = cKDTree(coords)
                idxs = tree.query_ball_point(query, r=R_THRESH)

                if len(idxs) == 0:
                    continue

                dists = np.linalg.norm(coords[idxs] - query, axis=1)

                for j in np.array(idxs)[np.argsort(dists)]:
                    prev_eddy = prev.iloc[j]

                    if (
                        pres_eddy['Cyc'] == prev_eddy['Cyc']
                        and prev_eddy['Eddy'] not in assigned
                    ):
                        best_id = prev_eddy['Eddy']
                        break

                if best_id is not None:
                    break

            mask = (
                df['Day'].eq(day)
                & df['eddy_idx'].eq(pres_eddy['eddy_idx'])
            )

            if best_id is not None:
                df.loc[mask, 'Eddy'] = best_id
                assigned.add(best_id)
            else:
                df.loc[mask, 'Eddy'] = next_num
                assigned.add(next_num)
                next_num += 1

        if day % 200 == 0:
            print(f"Day {day}, elapsed: {time.perf_counter() - tic:.2f}s")

    df = df.loc[df['Eddy'].ne(-1)].copy()
    df['next_num'] = next_num

    assert not df.duplicated(subset=['Eddy', 'Day']).any(), \
        "Duplicate (Eddy, Day) pairs found!"

    return df


# Tilt measurement
def compute_tilt_data(
    dic_all,
    eddy,
    num=6,
    depth_int=10,
    max_depth=1000,
    min_depth_range=200,
    var='Depth',
    eps=1e-10,
    min_points=5
):
    df_tilt_data = pd.DataFrame(columns=['Eddy', 'Day', 'TiltDis', 'TiltDir'])

    diffs_xc = {}
    diffs_yc = {}

    target_depths = np.arange(0, max_depth + depth_int, depth_int)
    full_idx = target_depths[:-1]

    dic = dic_all[f'Eddy{eddy}']

    day_nums = sorted(int(k[3:]) for k in dic.keys())
    all_day_nums = np.arange(min(day_nums), max(day_nums) + 1)
    days = [f'Day{d}' for d in all_day_nums]

    for d, day in enumerate(days):

        key = f't_{d}'

        if day not in dic:
            diffs_xc[key] = pd.Series(np.nan, index=full_idx)
            diffs_yc[key] = pd.Series(np.nan, index=full_idx)
            continue

        df = dic[day].copy()

        if len(df) < 2:
            diffs_xc[key] = pd.Series(np.nan, index=full_idx)
            diffs_yc[key] = pd.Series(np.nan, index=full_idx)
            continue

        if var == 'Depth':
            df[var] = np.abs(df[var])

        df = df[df[var] <= max_depth]
        df = df.set_index(var).sort_index()

        depths = df.index.values

        valid_depths = target_depths[
            (target_depths >= depths.min())
            & (target_depths <= depths.max())
        ]

        if len(valid_depths) < 2:
            diffs_xc[key] = pd.Series(np.nan, index=full_idx)
            diffs_yc[key] = pd.Series(np.nan, index=full_idx)
            continue

        xc_col = 'xc' if 'xc' in df.columns else 'x'
        yc_col = 'yc' if 'yc' in df.columns else 'y'

        xc_i = np.interp(valid_depths, depths, df[xc_col].values, left=np.nan, right=np.nan)
        yc_i = np.interp(valid_depths, depths, df[yc_col].values, left=np.nan, right=np.nan)

        diffs_xc[key] = pd.Series(np.diff(xc_i), index=valid_depths[:-1])
        diffs_yc[key] = pd.Series(np.diff(yc_i), index=valid_depths[:-1])

    df_X_all = pd.DataFrame(diffs_xc)
    df_Y_all = pd.DataFrame(diffs_yc)

    for ref_day in range(num // 2, len(days) - num // 2):

        df_X = df_X_all.iloc[:, ref_day - num // 2:ref_day + num // 2 + 1]
        df_Y = df_Y_all.iloc[:, ref_day - num // 2:ref_day + num // 2 + 1]

        if var == 'rho':
            df_X = df_X.dropna()
            df_Y = df_Y.dropna()

        df_data = pd.DataFrame(index=df_X.index)

        df_data['dxc'] = df_X.mean(axis=1)
        df_data['dyc'] = df_Y.mean(axis=1)

        df_data['sum_dxc'] = df_data['dxc'].cumsum()
        df_data['sum_dyc'] = df_data['dyc'].cumsum()

        df_data['var_dxc'] = df_X.var(axis=1)
        df_data['var_dyc'] = df_Y.var(axis=1)
        df_data['total_var'] = df_data['var_dxc'] + df_data['var_dyc']

        df_data['weight'] = 1 / (df_data['total_var'] + eps)
        df_data[var] = df_data.index

        df_data = df_data.replace([np.inf, -np.inf], np.nan)
        df_data = df_data.dropna(subset=['sum_dxc', 'sum_dyc', var, 'weight'])

        depth_range = df_data[var].max() - df_data[var].min() if len(df_data) else 0

        if depth_range < min_depth_range or len(df_data) < min_points:
            df_tilt_data.loc[len(df_tilt_data)] = {
                'Eddy': eddy,
                'Day': int(days[ref_day][3:]),
                'TiltDis': np.nan,
                'TiltDir': np.nan
            }
            continue

        xc = df_data['sum_dxc'].values
        yc = df_data['sum_dyc'].values
        z = df_data[var].values
        weights = df_data['weight'].values

        if not np.all(np.isfinite(weights)):
            df_tilt_data.loc[len(df_tilt_data)] = {
                'Eddy': eddy,
                'Day': int(days[ref_day][3:]),
                'TiltDis': np.nan,
                'TiltDir': np.nan
            }
            continue

        weights = weights / np.nanmax(weights)
        W = np.sum(weights)

        if not np.isfinite(W) or W <= 0:
            df_tilt_data.loc[len(df_tilt_data)] = {
                'Eddy': eddy,
                'Day': int(days[ref_day][3:]),
                'TiltDis': np.nan,
                'TiltDir': np.nan
            }
            continue

        mean = np.array([
            np.dot(weights, xc),
            np.dot(weights, yc),
            np.dot(weights, z)
        ]) / W

        X = np.vstack((xc, yc, z)).T
        Xc = X - mean
        Xw = Xc * np.sqrt(weights)[:, None]

        try:
            _, _, Vt = np.linalg.svd(Xw, full_matrices=False)
            direction = Vt[0]
        except Exception:
            df_tilt_data.loc[len(df_tilt_data)] = {
                'Eddy': eddy,
                'Day': int(days[ref_day][3:]),
                'TiltDis': np.nan,
                'TiltDir': np.nan
            }
            continue

        if not np.all(np.isfinite(direction)) or np.abs(direction[2]) < eps:
            df_tilt_data.loc[len(df_tilt_data)] = {
                'Eddy': eddy,
                'Day': int(days[ref_day][3:]),
                'TiltDis': np.nan,
                'TiltDir': np.nan
            }
            continue

        z_top = np.nanmin(z)
        z_btm = np.nanmax(z)

        t_top = (z_top - mean[2]) / direction[2]
        t_btm = (z_btm - mean[2]) / direction[2]

        p_top = mean + t_top * direction
        p_btm = mean + t_btm * direction

        tilt_dist = np.hypot(
            p_top[0] - p_btm[0],
            p_top[1] - p_btm[1]
        )

        tilt_dir = (bearing(p_btm, p_top) + 20) % 360

        df_tilt_data.loc[len(df_tilt_data)] = {
            'Eddy': eddy,
            'Day': int(days[ref_day][3:]),
            'TiltDis': tilt_dist,
            'TiltDir': tilt_dir
        }

    df_tilt_data['Day'] = df_tilt_data['Day'].astype(int)

    return df_tilt_data
    

def bearing(a, b):
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    angle_rad = np.arctan2(dx, dy)  # note the order: dx, dy
    angle_deg = np.degrees(angle_rad)
    bearing = (angle_deg + 360) % 360
    return bearing
    

def plot_tilt_method(
    dic_all, eddy, ref_day_idx,
    ax_x=None, ax_y=None,
    num=6,
    depth_int=10,
    max_depth=1000,
    min_depth_range=200,
    var='Depth',
    eps=1e-10,
    min_points=5,
    color='tab:blue',
    show=True
):
    dic = dic_all[f'Eddy{eddy}']

    day_nums = sorted(int(k[3:]) for k in dic.keys())
    all_day_nums = np.arange(min(day_nums), max(day_nums) + 1)
    days = [f'Day{d}' for d in all_day_nums]

    target_depths = np.arange(0, max_depth + depth_int, depth_int)
    full_idx = target_depths[:-1]

    diffs_xc = {}
    diffs_yc = {}
    interp_xc = {}
    interp_yc = {}

    for d, day in enumerate(days):

        key = f't_{d}'

        if day not in dic:
            diffs_xc[key] = pd.Series(np.nan, index=full_idx)
            diffs_yc[key] = pd.Series(np.nan, index=full_idx)
            interp_xc[key] = pd.Series(np.nan, index=target_depths)
            interp_yc[key] = pd.Series(np.nan, index=target_depths)
            continue

        df = dic[day].copy()

        if len(df) < 2:
            diffs_xc[key] = pd.Series(np.nan, index=full_idx)
            diffs_yc[key] = pd.Series(np.nan, index=full_idx)
            interp_xc[key] = pd.Series(np.nan, index=target_depths)
            interp_yc[key] = pd.Series(np.nan, index=target_depths)
            continue

        if var == 'Depth':
            df[var] = np.abs(df[var])

        df = df[df[var] <= max_depth]
        df = df.set_index(var).sort_index()

        depths = df.index.values

        valid_depths = target_depths[
            (target_depths >= depths.min())
            & (target_depths <= depths.max())
        ]

        if len(valid_depths) < 2:
            diffs_xc[key] = pd.Series(np.nan, index=full_idx)
            diffs_yc[key] = pd.Series(np.nan, index=full_idx)
            interp_xc[key] = pd.Series(np.nan, index=target_depths)
            interp_yc[key] = pd.Series(np.nan, index=target_depths)
            continue

        xc_col = 'xc' if 'xc' in df.columns else 'x'
        yc_col = 'yc' if 'yc' in df.columns else 'y'

        xc_i = np.interp(
            valid_depths,
            depths,
            df[xc_col].values,
            left=np.nan,
            right=np.nan
        )

        yc_i = np.interp(
            valid_depths,
            depths,
            df[yc_col].values,
            left=np.nan,
            right=np.nan
        )

        diffs_xc[key] = pd.Series(np.diff(xc_i), index=valid_depths[:-1])
        diffs_yc[key] = pd.Series(np.diff(yc_i), index=valid_depths[:-1])

        interp_xc[key] = pd.Series(
            np.r_[0, np.cumsum(np.diff(xc_i))],
            index=valid_depths
        )

        interp_yc[key] = pd.Series(
            np.r_[0, np.cumsum(np.diff(yc_i))],
            index=valid_depths
        )

    df_X_all = pd.DataFrame(diffs_xc)
    df_Y_all = pd.DataFrame(diffs_yc)

    i0 = ref_day_idx - num // 2
    i1 = ref_day_idx + num // 2 + 1

    if i0 < 0 or i1 > len(days):
        print(f'Eddy {eddy}: ref_day_idx outside valid compute_tilt_data range.')
        return None, None, pd.DataFrame()

    df_X = df_X_all.iloc[:, i0:i1]
    df_Y = df_Y_all.iloc[:, i0:i1]

    if var == 'rho':
        df_X = df_X.dropna()
        df_Y = df_Y.dropna()

    df_data = pd.DataFrame(index=df_X.index)

    df_data['dxc'] = df_X.mean(axis=1)
    df_data['dyc'] = df_Y.mean(axis=1)

    df_data['sum_dxc'] = df_data['dxc'].cumsum()
    df_data['sum_dyc'] = df_data['dyc'].cumsum()

    df_data['var_dxc'] = df_X.var(axis=1)
    df_data['var_dyc'] = df_Y.var(axis=1)
    df_data['total_var'] = df_data['var_dxc'] + df_data['var_dyc']

    df_data['weight'] = 1 / (df_data['total_var'] + eps)
    df_data[var] = df_data.index

    df_data = df_data.replace([np.inf, -np.inf], np.nan)
    df_data = df_data.dropna(subset=['sum_dxc', 'sum_dyc', var, 'weight'])

    depth_range = df_data[var].max() - df_data[var].min() if len(df_data) else 0

    if depth_range < min_depth_range:
        print(f'Eddy {eddy}: tilt not measurable, depth range < {min_depth_range} m.')
        return None, None, df_data

    if len(df_data) < min_points:
        print(f'Eddy {eddy}: tilt not measurable, too few valid points.')
        return None, None, df_data

    xc = df_data['sum_dxc'].values
    yc = df_data['sum_dyc'].values
    z = df_data[var].values
    weights = df_data['weight'].values

    if not np.all(np.isfinite(weights)):
        print(f'Eddy {eddy}: tilt not measurable, invalid weights.')
        return None, None, df_data

    weights = weights / np.nanmax(weights)
    W = np.sum(weights)

    if not np.isfinite(W) or W <= 0:
        print(f'Eddy {eddy}: tilt not measurable, invalid total weight.')
        return None, None, df_data

    mean = np.array([
        np.dot(weights, xc),
        np.dot(weights, yc),
        np.dot(weights, z)
    ]) / W

    X = np.vstack((xc, yc, z)).T
    Xc = X - mean
    Xw = Xc * np.sqrt(weights)[:, None]

    try:
        _, _, Vt = np.linalg.svd(Xw, full_matrices=False)
        direction = Vt[0]
    except Exception:
        print(f'Eddy {eddy}: tilt not measurable, SVD failed.')
        return None, None, df_data

    if not np.all(np.isfinite(direction)) or np.abs(direction[2]) < eps:
        print(f'Eddy {eddy}: tilt not measurable, invalid SVD direction.')
        return None, None, df_data

    z_top = np.nanmin(z)
    z_btm = np.nanmax(z)

    t_top = (z_top - mean[2]) / direction[2]
    t_btm = (z_btm - mean[2]) / direction[2]

    p_top = mean + t_top * direction
    p_btm = mean + t_btm * direction

    tilt_dist = np.hypot(
        p_top[0] - p_btm[0],
        p_top[1] - p_btm[1]
    )

    tilt_dir = (bearing(p_btm, p_top) + 20) % 360

    if ax_x is None or ax_y is None:
        fig, axs = plt.subplots(1, 2, figsize=(7, 5), sharey=True)
        ax_x, ax_y = axs
    else:
        fig = ax_x.figure

    df_xi = pd.DataFrame(interp_xc).iloc[:, i0:i1]
    df_yi = pd.DataFrame(interp_yc).iloc[:, i0:i1]

    for col in df_xi.columns:
        ax_x.plot(df_xi[col], df_xi.index, color=color, alpha=0.25, lw=0.8)
        ax_y.plot(df_yi[col], df_yi.index, color=color, alpha=0.25, lw=0.8)

    spread_x = df_data['var_dxc'].values
    spread_y = df_data['var_dyc'].values

    ax_x.fill_betweenx(
        z,
        xc - spread_x,
        xc + spread_x,
        color=color,
        alpha=0.25
    )

    ax_y.fill_betweenx(
        z,
        yc - spread_y,
        yc + spread_y,
        color=color,
        alpha=0.25
    )

    ax_x.plot(xc, z, color=color, lw=2, label='mean cumulative displacement')
    ax_y.plot(yc, z, color=color, lw=2)

    ax_x.plot(
        [p_top[0], p_btm[0]],
        [p_top[2], p_btm[2]],
        'k--',
        lw=2,
        label='weighted SVD tilt'
    )

    ax_y.plot(
        [p_top[1], p_btm[1]],
        [p_top[2], p_btm[2]],
        'k--',
        lw=2
    )

    ax_x.invert_yaxis()

    ax_x.set_xlabel('x displacement')
    ax_y.set_xlabel('y displacement')
    ax_x.set_ylabel(var)

    day_actual = int(days[ref_day_idx][3:])

    ax_x.set_title(f'Eddy {eddy}, Day {day_actual}: x-depth')
    ax_y.set_title(f'Tilt = {tilt_dist:.1f} m, Dir = {tilt_dir:.1f}°')

    ax_x.legend(frameon=False)

    if show:
        plt.tight_layout()
        plt.show()

    return fig, (ax_x, ax_y), df_data

def phys_grad(F, X, Y, mask=None):
    # index-space gradients
    x_i, x_j = np.gradient(X)
    y_i, y_j = np.gradient(Y)
    F_i, F_j = np.gradient(F)
    # Jacobian
    J = x_i*y_j - x_j*y_i
    # physical gradients
    dFdx = ( F_i*y_j - F_j*y_i) / J
    dFdy = (-F_i*x_j + F_j*x_i) / J
    # handle singular cells
    bad = np.isclose(J, 0)
    dFdx[bad] = np.nan
    dFdy[bad] = np.nan
    # apply mask (1=ocean, 0=land)
    if mask is not None:
        m = mask.astype(bool)
        dFdx = np.where(m, dFdx, np.nan)
        dFdy = np.where(m, dFdy, np.nan)
    return dFdx, dFdy

def compute_core_mean(
    df_eddies,
    X_grid,
    Y_grid,
    mask_rho,
    base_path=None,
    varname=None,
    fixed_field=None,
    colname=None,
    circle_region_flag=False
):
    """
    Core-mean of either
      - a 3D field (x,y,t) loaded as <varname>_<fnumber>.npy, or
      - a fixed 2D field (x,y) passed in as fixed_field.
    """
    if fixed_field is None and (base_path is None or varname is None):
        raise ValueError("Either fixed_field OR (base_path and varname) must be provided.")
    mode_2d = fixed_field is not None
    if colname is None:
        if mode_2d:
            colname = "field_core"
        else:
            colname = f"{varname}"
    df = df_eddies[~df_eddies["TiltDis"].isna()].copy()
    chunks = []
    if mode_2d:
        field2d = np.where(mask_rho, fixed_field, np.nan)
    for fname, df_loc in df.groupby("fname"):
        if not mode_2d:
            fnumber  = int(fname[-8:-3])
            base_day = fnumber + 1
            data3d = np.load(f"{base_path}/{varname}_{fnumber:05}.npy")
            data3d = np.where(mask_rho[:, :, None], data3d, np.nan)
        df_loc = df_loc.copy().reset_index(drop=False)
        core_vals = np.full(len(df_loc), np.nan)
        for idx, row in enumerate(df_loc.itertuples(index=False)):
            dx = X_grid - row.xc
            dy = Y_grid - row.yc
            if circle_region_flag:
                if hasattr(row, 'rmax') and np.isfinite(row.rmax):
                    rho2 = (dx**2 + dy**2)
                    core_mask = rho2 <= row.rmax**2
                else:
                    rho2 = np.full_like(dx, np.nan, dtype=float)
            else:
                if hasattr(row, 'q11') and np.isfinite(row.q11):
                    rho2 = (
                        row.q11 * dx**2
                        + 2 * row.q12 * dx * dy
                        + row.q22 * dy**2
                    )
                elif isinstance(row.Q, np.ndarray) and row.Q.shape == (2, 2) and np.isfinite(row.Q).all():
                    rho2 = (
                        row.Q[0, 0] * dx**2
                        + 2 * row.Q[1, 0] * dx * dy
                        + row.Q[1, 1] * dy**2
                    )
                else:
                    rho2 = np.full_like(dx, np.nan, dtype=float)
                core_mask = rho2 <= row.Rc**2 / 2
            if not core_mask.any():
                continue
            if mode_2d:
                vals = field2d[core_mask]
            else:
                t_idx = int(row.Day - base_day)
                vals = data3d[:, :, t_idx][core_mask]
            core_vals[idx] = np.nanmean(vals)
        chunks.append(pd.DataFrame({
            "Eddy": df_loc["Eddy"].to_numpy(),
            "Day":  df_loc["Day"].to_numpy(),
            colname: core_vals
        }))
    df_core = pd.concat(chunks, ignore_index=True)
    df_out = df_eddies.merge(
        df_core[["Eddy", "Day", colname]],
        how="left",
        on=["Eddy", "Day"]
    )
    return df_out

def _nice_step(h, base):
    s = h / base
    for k in [1, 2, 2.5, 5, 10]:
        if s <= k: return k * base
    return np.ceil(s) * base

def _grid_step(G):
    gx = np.diff(np.sort(np.unique(G.ravel())))
    return np.nanmedian(gx[gx > 0])

def bin_edges_fd(x, xgrid, rule='fd'): # Freedman-Diaconis (fg) rationale
    n = len(x)
    if n < 2: return np.array([np.min(x), np.max(x)])
    rng = np.ptp(x)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    std = np.std(x, ddof=1)

    # raw width (km)
    if rule.lower() == 'fd':
        h = 2 * (iqr if iqr > 0 else 1.349*std) / (n ** (1/3))
    else:  # 'scott'
        h = 3.5 * std / (n ** (1/3))

    # fallback if degenerate
    if not np.isfinite(h) or h <= 0:
        h = rng / max(10, np.sqrt(n))

    # snap to grid spacing
    base = _grid_step(xgrid)
    h = _nice_step(h, base)

    lo = np.floor(np.min(x) / h) * h
    hi = np.ceil(np.max(x) / h) * h
    return np.arange(lo, hi + h, h)

def binned_median(x, y, v, xbins, ybins):
    ix = np.digitize(x, xbins) - 1
    iy = np.digitize(y, ybins) - 1

    nx, ny = len(xbins) - 1, len(ybins) - 1
    ok = (
        (ix >= 0) & (ix < nx) &
        (iy >= 0) & (iy < ny) &
        np.isfinite(v)
    )

    bins = {}
    for i, j, val in zip(ix[ok], iy[ok], v[ok]):
        bins.setdefault((j, i), []).append(val)

    out = np.full((ny, nx), np.nan)
    for (j, i), vals in bins.items():
        out[j, i] = np.nanmedian(vals)

    return out

def plot_binned_median_map(
    df_eddies,
    metric='Rc',
    X_grid=X_grid,
    Y_grid=Y_grid,
    h=h,
    mask_rho=mask_rho,
    lat_rho=lat_rho,
    lon_rho=lon_rho,
    vmin=0,
    vmax=120,
    rule='fd',
    levels_lat=[-40, -35, -30, -25],
    levels_lon=[150, 155, 160],
    cmaps={'AE': 'Reds', 'CE': 'Blues'},
    units='km',
    figsize=(9, 8)
):

    xbins = bin_edges_fd(pd.to_numeric(df_eddies.xc, errors='coerce').to_numpy(dtype=float), X_grid, rule=rule)
    ybins = bin_edges_fd(pd.to_numeric(df_eddies.yc, errors='coerce').to_numpy(dtype=float), Y_grid, rule=rule)

    norm = Normalize(vmin=vmin, vmax=vmax)

    fig, axs = plt.subplots(1, 2, figsize=figsize, sharey=True)

    for ax, cyc in zip(axs, ['AE', 'CE']):

        df = (
            df_eddies[df_eddies.Cyc.eq(cyc)]
            .dropna(subset=['xc', 'yc', metric])
            .sort_values(metric, kind='mergesort', ignore_index=True)
        )

        ax.contour(X_grid, Y_grid, h, levels=[4000], colors='k')

        H = binned_median(
            df.xc.to_numpy(dtype=float),
            df.yc.to_numpy(dtype=float),
            df[metric].to_numpy(dtype=float),
            xbins,
            ybins
        )

        m = ax.pcolormesh(
            xbins, ybins, H,
            cmap=cmaps[cyc],
            norm=norm,
            shading='auto',
            rasterized=True
        )

        cb = fig.colorbar(m, ax=ax, location='top', shrink=0.9, pad=0.02)
        cb.set_label(fr'{cyc} median surface ${metric}$ ({units})', fontsize=12)
        cb.set_ticks(np.linspace(vmin, vmax, 5))

        ax.contourf(
            X_grid, Y_grid, np.where(mask_rho == 0, 1, np.nan),
            levels=[0.5, 1.5], colors=['k'], alpha=0.5
        )

        c1 = ax.contour(X_grid, Y_grid, lat_rho, levels=levels_lat, colors='k', linewidths=0.5)
        ax.clabel(c1, fmt=lambda v: f"{-v:.0f}°S", inline=True, colors='k')

        c2 = ax.contour(X_grid, Y_grid, lon_rho, levels=levels_lon, colors='k', linewidths=0.5)
        ax.clabel(c2, fmt=lambda v: f"{v:.0f}°E", inline=True, colors='k')

        ax.axis('equal')
        ax.set_xlim(15, X_grid.max())
        ax.set_ylim(Y_grid.min(), Y_grid.max())
        ax.set_xlabel('x (km)', fontsize=11)

    axs[0].set_ylabel('y (km)', fontsize=11)

    plt.tight_layout()
    plt.show()

    return fig, axs

def tilt_distance_Li(x, y, z, zmin=None, zmax=None):
    x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if zmin is not None: m &= (z >= zmin)
    if zmax is not None: m &= (z <= zmax)
    x, y = x[m], y[m]
    if x.size == 0: return np.nan, np.nan, (np.nan, np.nan)

    TDx = np.nanmax(x) - np.nanmin(x)
    TDy = np.nanmax(y) - np.nanmin(y)
    TD  = np.hypot(TDx, TDy)
    theta_deg = np.degrees(np.arctan2(TDy, TDx))
    return TD, theta_deg, (TDx, TDy)


##### Composite eddy ######
def smooth_esp_params_depth(
    df,
    window=5,
    min_periods=2,
    centre=True,
    remove_outliers=True,
    z_col='Depth'
):
    """
    Smooth ESP parameters with depth for one eddy-day dataframe.
    Smooths xc, yc, q11, q12, q22, Omega, Rc.
    """

    df = df.copy().sort_values(z_col)

    # unpack Q into scalar columns
    df['q11'] = df['Q'].apply(lambda Q: Q[0, 0] if Q is not None else np.nan)
    df['q12'] = df['Q'].apply(lambda Q: Q[0, 1] if Q is not None else np.nan)
    df['q22'] = df['Q'].apply(lambda Q: Q[1, 1] if Q is not None else np.nan)

    cols = ['xc', 'yc', 'q11', 'q12', 'q22', 'Omega', 'Rc']

    for col in cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    if remove_outliers:
        for col in cols:
            med = df[col].rolling(window, center=centre, min_periods=min_periods).median()
            mad = (df[col] - med).abs().rolling(window, center=centre, min_periods=min_periods).median()

            outlier = (df[col] - med).abs() > 4 * 1.4826 * mad
            df.loc[outlier, col] = np.nan

    # interpolate then smooth
    for col in cols:
        df[col] = (
            df[col]
            .interpolate(limit_direction='both')
            .rolling(window, center=centre, min_periods=min_periods)
            .mean()
        )

    return df

def composite_eddy_velocity(
    eddy,
    dic_vert,
    z_r,
    xlim_km=150,
    nx=100,
    ny=100,
    zmax=1.2e3,
    smooth_params=True,
    smooth_window=5,
    plot=False,
    plot_panel='both',      # 'both', 'zonal', or 'meridional'
    ax=None,                # pass external axis for single-panel plotting
    add_cbar=True,
    zlim=1000,
    levels=20,
    cmap='RdBu_r',
):

    x = np.linspace(-xlim_km, xlim_km, nx) * 1e3
    y = np.linspace(-xlim_km, xlim_km, ny) * 1e3

    z_grid = np.insert(np.abs(z_r[150, 150, 1:]), 0, 0)
    z_grid = z_grid[z_grid < zmax]

    X, Y = np.meshgrid(x, y, indexing='ij')
    nx, ny, nz = len(x), len(y), len(z_grid)

    U_days = []
    V_days = []

    last_w = np.nan

    for day, df in dic_vert[f'Eddy{eddy}'].items():

        df = df.copy().sort_values('Depth')

        if smooth_params:
            df = smooth_esp_params_depth(
                df,
                window=smooth_window,
                z_col='Depth'
            )
        else:
            df['q11'] = df['Q'].apply(lambda Q: Q[0, 0] if Q is not None else np.nan)
            df['q12'] = df['Q'].apply(lambda Q: Q[0, 1] if Q is not None else np.nan)
            df['q22'] = df['Q'].apply(lambda Q: Q[1, 1] if Q is not None else np.nan)

        if 'w' in df.columns:
            last_w = df['w'].iloc[0]

        df['xc'] -= df['xc'].iloc[0]
        df['yc'] -= df['yc'].iloc[0]

        U = np.zeros((nx, ny, nz))
        V = np.zeros((nx, ny, nz))

        for _, data in df.iterrows():

            if not np.all(np.isfinite([
                data.xc, data.yc,
                data.q11, data.q12, data.q22,
                data.Omega, data.Rc
            ])):
                continue

            if data.Rc <= 0:
                continue

            k = np.argmin(np.abs(z_grid - abs(data.Depth)))

            dx = X - data.xc * 1e3
            dy = Y - data.yc * 1e3

            rho2 = (
                data.q11 * dx**2
                + 2 * data.q12 * dx * dy
                + data.q22 * dy**2
            )

            fac = data.Omega * np.exp(-rho2 / (data.Rc * 1e3)**2)

            U[:, :, k] = -fac * (data.q12*dx + data.q22*dy)
            V[:, :, k] =  fac * (data.q11*dx + data.q12*dy)

        U = np.nan_to_num(U, nan=0.0, posinf=0.0, neginf=0.0)
        V = np.nan_to_num(V, nan=0.0, posinf=0.0, neginf=0.0)

        U_days.append(U)
        V_days.append(V)

    U_comp = np.mean(U_days, axis=0)
    V_comp = np.mean(V_days, axis=0)

    if not plot:
        return X, Y, z_grid, U_comp, V_comp

    ix0 = np.argmin(np.abs(x))
    iy0 = np.argmin(np.abs(y))

    vmax = np.nanmax(np.abs([U_comp, V_comp]))
    clevels = np.linspace(-vmax, vmax, levels + 1)

    cyc = 'AE' if np.sign(last_w) > 0 else 'CE'

    def _plot_zonal(ax):
        m = ax.contourf(
            x / 1e3,
            z_grid,
            V_comp[:, iy0, :].T,
            cmap=cmap,
            levels=clevels,
            extend='both'
        )

        ax.contour(
            x / 1e3,
            z_grid,
            V_comp[:, iy0, :].T,
            levels=[0],
            colors='k',
            linewidths=2,
            alpha=.7
        )

        ax.set_title('Zonal: $v$')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('Depth (m)')
        ax.axvline(0, color='k', lw=0.8, ls='--', alpha=.8)
        ax.set_ylim(zlim, 0)

        return m

    def _plot_meridional(ax):
        m = ax.contourf(
            y / 1e3,
            z_grid,
            U_comp[ix0, :, :].T,
            cmap=cmap,
            levels=clevels,
            extend='both'
        )

        ax.contour(
            y / 1e3,
            z_grid,
            U_comp[ix0, :, :].T,
            levels=[0],
            colors='k',
            linewidths=2,
            alpha=.7
        )

        ax.set_title('Meridional: $u$')
        ax.set_xlabel('y (km)')
        ax.set_ylabel('Depth (m)')
        ax.axvline(0, color='k', lw=0.8, ls='--', alpha=.8)
        ax.set_ylim(zlim, 0)

        return m

    if plot_panel == 'zonal':

        if ax is None:
            fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)
        else:
            fig = ax.figure

        m = _plot_zonal(ax)
        ax.set_title(f'{cyc}{eddy}: zonal $v$')

        if add_cbar:
            cbar = fig.colorbar(m, ax=ax, shrink=0.9)
            cbar.set_label('v (m s$^{-1}$)')
            cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        return X, Y, z_grid, U_comp, V_comp, fig, ax

    if plot_panel == 'meridional':

        if ax is None:
            fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)
        else:
            fig = ax.figure

        m = _plot_meridional(ax)
        ax.set_title(f'{cyc}{eddy}: meridional $u$')

        if add_cbar:
            cbar = fig.colorbar(m, ax=ax, shrink=0.9)
            cbar.set_label('u (m s$^{-1}$)')
            cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        return X, Y, z_grid, U_comp, V_comp, fig, ax

    fig, axs = plt.subplots(
        1, 2,
        figsize=(11, 4),
        sharey=True,
        constrained_layout=True
    )

    m0 = _plot_zonal(axs[0])
    _plot_meridional(axs[1])

    axs[1].set_ylabel('')

    if add_cbar:
        cbar = fig.colorbar(m0, ax=axs, location='right', shrink=0.9)
        cbar.set_label('(m s$^{-1}$)')

    fig.suptitle(f'{cyc}{eddy}')

    return X, Y, z_grid, U_comp, V_comp, fig, axs

# Regioning 
def make_region_grids(
    X_grid,
    Y_grid,
    lon_rho,
    lat_rho,
    h,
    mask_rho,
    lon_split=157,
    lat_split=-33,
    shelf_hmax=4000,
    shelf_xmax=400,
    shelf_lonmax=154.85,
):
    """
    Create region_mask_grid and bin_grid for regions:
    S1, S2, U1, D1, U2, D2.
    """

    region_mask_grid = (
        (h < shelf_hmax)
        & (X_grid < shelf_xmax)
        & (lon_rho < shelf_lonmax)
        & (mask_rho == 1)
    )

    bin_grid = np.full(X_grid.shape, np.nan)

    # Shelf
    bin_grid[region_mask_grid & (lat_rho >= lat_split)] = 1  # S1
    bin_grid[region_mask_grid & (lat_rho <  lat_split)] = 2  # S2

    # Offshore west
    bin_grid[
        (~region_mask_grid)
        & (lon_rho < lon_split)
        & (mask_rho == 1)
        & (lat_rho >= lat_split)
    ] = 3  # U1

    bin_grid[
        (~region_mask_grid)
        & (lon_rho < lon_split)
        & (mask_rho == 1)
        & (lat_rho < lat_split)
    ] = 4  # D1

    # Offshore east
    bin_grid[
        (lon_rho >= lon_split)
        & (mask_rho == 1)
        & (lat_rho >= lat_split)
    ] = 5  # U2

    bin_grid[
        (lon_rho >= lon_split)
        & (mask_rho == 1)
        & (lat_rho < lat_split)
    ] = 6  # D2

    return region_mask_grid, bin_grid

def plot_region_map(
    ax,
    X_grid,
    Y_grid,
    lon_rho,
    lat_rho,
    h,
    mask_rho,
    levels_lat,
    levels_lon,
    lon_split=157,
    lat_split=-33,
    title=None,
    borders_only=False
):
    """
    Plot the six-region map onto an existing axis.
    """

    region_mask_grid, bin_grid = make_region_grids(
        X_grid,
        Y_grid,
        lon_rho,
        lat_rho,
        h,
        mask_rho,
        lon_split=lon_split,
        lat_split=lat_split
    )

    if not borders_only:
        ax.contourf(
            X_grid, Y_grid,
            np.where(mask_rho == 0, 1, np.nan),
            levels=[0.5, 1.5],
            colors=['k'],
            alpha=0.5
        )

    
        ax.contourf(
            X_grid, Y_grid,
            bin_grid,
            levels=np.arange(0.5, 7.5, 1),
            alpha=0.25,
            cmap='gist_rainbow'
        )

        c1 = ax.contour(
            X_grid, Y_grid,
            lat_rho,
            levels=levels_lat,
            colors='k',
            linewidths=0.5
        )
    
        ax.clabel(
            c1,
            fmt=lambda v: f"{np.abs(v):.0f}°S",
            inline=True,
            colors='k'
        )
    
        c2 = ax.contour(
            X_grid, Y_grid,
            lon_rho,
            levels=levels_lon,
            colors='k',
            linewidths=0.5
        )
    
        ax.clabel(
            c2,
            fmt=lambda v: f"{v:.0f}°E",
            inline=True,
            colors='k'
        )
    
        ax.contour(
            X_grid, Y_grid,
            h,
            levels=[4000],
            colors='k',
            linewidths=1
        )

    ax.contour(
        X_grid, Y_grid,
        region_mask_grid.astype(float),
        levels=[0.5],
        colors='magenta',
        linewidths=2,
        linestyles='-'
    )

    ax.contour(
        X_grid, Y_grid,
        lon_rho,
        levels=[lon_split],
        colors='magenta',
        linewidths=2,
        linestyles='-'
    )

    ax.contour(
        X_grid, Y_grid,
        np.where(mask_rho, lat_rho, np.nan),
        levels=[lat_split],
        colors='magenta',
        linewidths=2,
        linestyles='-'
    )

    labels = [
        ('S1', 220, 1300),
        ('S2', 120, 50),
        ('U1', 400, 1450),
        ('U2', 800, 1450),
        ('D1', 400, 700),
        ('D2', 800, 700),
    ]

    for txt, x, y in labels:
        ax.text(
            x, y, txt,
            ha='center',
            va='center',
            fontsize=11,
            fontweight='bold',
            color='m' if borders_only else 'k'
        )

    ax.set_aspect('equal')
    ax.set_xlim(X_grid.min(), X_grid.max())
    ax.set_ylim(Y_grid.min(), Y_grid.max())

    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')

    if title is not None:
        ax.set_title(title)

    return ax

def add_region_column(
    df,
    X_grid,
    Y_grid,
    lon_rho,
    lat_rho,
    h,
    mask_rho,
    lon_split=157,
    lat_split=-33,
):
    """
    Add a Region column (S1, S2, U1, D1, U2, D2)
    to a dataframe containing xc and yc coordinates.
    """

    region_mask_grid, bin_grid = make_region_grids(
        X_grid,
        Y_grid,
        lon_rho,
        lat_rho,
        h,
        mask_rho,
        lon_split=lon_split,
        lat_split=lat_split
    )

    tree = cKDTree(
        np.column_stack([X_grid.ravel(), Y_grid.ravel()])
    )

    _, idx = tree.query(
        np.column_stack([df.xc, df.yc])
    )

    region_map = {
        1: 'S1',
        2: 'S2',
        3: 'U1',
        4: 'D1',
        5: 'U2',
        6: 'D2'
    }

    df = df.copy()

    df['Region'] = (
        pd.Series(
            bin_grid.ravel()[idx],
            index=df.index
        )
        .map(region_map)
    )

    return df

# def plot_region_map(
#     ax,
#     X_grid,
#     Y_grid,
#     lon_rho,
#     lat_rho,
#     h,
#     mask_rho,
#     bin_grid,
#     region_mask_grid,
#     levels_lat,
#     levels_lon,
#     lon_split=157,
#     lat_split=-33,
#     title=None,
#     borders_only=False
# ):
#     """
#     Plot the six-region map onto an existing axis.
#     """

#     ax.contourf(
#         X_grid, Y_grid,
#         np.where(mask_rho == 0, 1, np.nan),
#         levels=[0.5, 1.5],
#         colors=['k'],
#         alpha=0.5
#     )

#     ax.contourf(
#         X_grid, Y_grid,
#         bin_grid,
#         levels=np.arange(0.5, 7.5, 1),
#         alpha=0.25,
#         cmap='gist_rainbow'
#     )

#     c1 = ax.contour(
#         X_grid, Y_grid,
#         lat_rho,
#         levels=levels_lat,
#         colors='k',
#         linewidths=0.5
#     )

#     ax.clabel(
#         c1,
#         fmt=lambda v: f"{np.abs(v):.0f}°S",
#         inline=True,
#         colors='k'
#     )

#     c2 = ax.contour(
#         X_grid, Y_grid,
#         lon_rho,
#         levels=levels_lon,
#         colors='k',
#         linewidths=0.5
#     )

#     ax.clabel(
#         c2,
#         fmt=lambda v: f"{v:.0f}°E",
#         inline=True,
#         colors='k'
#     )

#     ax.contour(
#         X_grid, Y_grid,
#         h,
#         levels=[4000],
#         colors='k',
#         linewidths=1
#     )

#     ax.contour(
#         X_grid, Y_grid,
#         region_mask_grid.astype(float),
#         levels=[0.5],
#         colors='magenta',
#         linewidths=2,
#         linestyles='-'
#     )

#     ax.contour(
#         X_grid, Y_grid,
#         lon_rho,
#         levels=[lon_split],
#         colors='magenta',
#         linewidths=2,
#         linestyles='-'
#     )

#     ax.contour(
#         X_grid, Y_grid,
#         np.where(mask_rho, lat_rho, np.nan),
#         levels=[lat_split],
#         colors='magenta',
#         linewidths=2,
#         linestyles='-'
#     )

#     labels = [
#         ('S1', 220, 1300),
#         ('S2', 120, 50),
#         ('U1', 400, 1450),
#         ('U2', 800, 1450),
#         ('D1', 400, 700),
#         ('D2', 800, 700),
#     ]

#     for txt, x, y in labels:
#         ax.text(
#             x, y, txt,
#             ha='center',
#             va='center',
#             fontsize=11,
#             fontweight='bold'
#         )

#     ax.set_aspect('equal')
#     ax.set_xlim(X_grid.min(), X_grid.max())
#     ax.set_ylim(Y_grid.min(), Y_grid.max())

#     ax.set_xlabel('x (km)')
#     ax.set_ylabel('y (km)')

#     if title is not None:
#         ax.set_title(title)

#     return ax

