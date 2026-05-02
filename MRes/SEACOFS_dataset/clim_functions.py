import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import time
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

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


# Tracking

def collect_tracking_R(
    df_data,
    L_SCALE=50,
    W_SCALE=1e-5,
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
                    prev['w'].values / W_SCALE
                ])

                tree = cKDTree(coords)

                query = np.array([
                    pres_eddy['xc'] / L_SCALE,
                    pres_eddy['yc'] / L_SCALE,
                    pres_eddy['w'] / W_SCALE
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
                        'dw': pres_eddy['w'] - prev_eddy['w'],
                        'Cyc': pres_eddy['Cyc'],
                    })

                # only diagnose most recent available day
                break

    return pd.DataFrame(rows)


def tracking_kdtree(
    df_data,
    start_ID,
    next_num,
    L_SCALE=50,
    W_SCALE=1e-5,
    R_THRESH=1,
    LOOKBACK=4
):
    tic = time.perf_counter()
    df = df_data.dropna(subset=['xc', 'yc', 'w']).copy()

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
                pres_eddy['w'] / W_SCALE
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
                    prev['w'].values / W_SCALE
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
