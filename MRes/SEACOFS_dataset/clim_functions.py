import numpy as np

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