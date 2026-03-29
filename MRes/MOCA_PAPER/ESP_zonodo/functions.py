import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

# SOLO 
def solo(l, VT, VN, core_thresh=30_000, plot=False, ax=None):
    from scipy.optimize import curve_fit
    '''
    l: The along track distance, 0 at the beginning
    VT: The tangential velocity i.e., the along-track velocity 
    VN: The normal velocity i.e., the across-track velocity 
    '''
    
    l, VT, VN = map(np.asarray, (l, VT, VN))
    m = np.isfinite(l) & np.isfinite(VT) & np.isfinite(VN)
    l, VT, VN = l[m], VT[m], VN[m]

    if l.size < 4:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    def vn_model(l, x0, C, D):
        dl = l - x0
        return C * dl + D * dl**3

    def vt_model(l, x0, A, B):
        dl = l - x0
        return A + B * dl**2

    def fit_x0_from_vn(l, VN):
        x0_guess = l[np.argmin(np.abs(VN))]
        C_guess = 0.0
        D_guess = 0.0

        try:
            popt, _ = curve_fit(
                vn_model, l, VN,
                p0=[x0_guess, C_guess, D_guess],
                maxfev=10000
            )
            return popt
        except Exception:
            return np.array([np.nan, np.nan, np.nan])

    x0, C, D = fit_x0_from_vn(l, VN)
    if not np.isfinite(x0):
        return np.nan, np.nan, np.nan, np.nan, np.nan

    # focus on core
    mask = np.abs(l - x0) <= core_thresh
    l, VT, VN = l[mask], VT[mask], VN[mask]

    # refit
    x0, C, D = fit_x0_from_vn(l, VN)
    if not np.isfinite(x0):
        return np.nan, np.nan, np.nan, np.nan, np.nan

    dl = l - x0

    Xu = np.vstack([np.ones_like(dl), dl**2]).T
    A, B = np.linalg.lstsq(Xu, VT, rcond=None)[0]

    l0 = x0 # the eddy center's along track position
    r0 = A / C # the eddy center's across track position
    Q = np.array([[1., 0.], [0., 1.]])
    Omega =  C 
    w = 2 * Omega

    if plot:
        lfit = np.linspace(l.min(), l.max(), 500)

        if ax is None:
            fig, axs = plt.subplots(1, 2, figsize=(10, 4))
        else:
            axs = ax

        axs = np.asarray(axs).ravel()

        axs[0].scatter(l, VN)
        axs[0].plot(lfit, vn_model(lfit, x0, C, D))
        axs[0].axvline(x0, linestyle='--')
        axs[0].axhline(0, linewidth=0.8)
        axs[0].set_xlabel('l')
        axs[0].set_ylabel('VN')
        axs[0].set_title('VN fit')

        axs[1].scatter(l, VT)
        axs[1].plot(lfit, vt_model(lfit, x0, A, B))
        axs[1].axvline(x0, linestyle='--')
        axs[1].set_xlabel('l')
        axs[1].set_ylabel('VT')
        axs[1].set_title('VT fit')

        if ax is None:
            plt.tight_layout()
            plt.show()

    return l0, r0, w, Q, Omega

# Project quasi-straight data to transect data i.e., length along transect, along-track and across-track velocities
def project_sadcp_to_transect(x, y, u, v):
    x, y, u, v = map(np.asarray, (x, y, u, v))
    msk = np.isfinite(x) & np.isfinite(y) & np.isfinite(u) & np.isfinite(v)
    x, y, u, v = x[msk], y[msk], u[msk], v[msk]
    n = x.size
    if n == 0:
        return pd.DataFrame(columns=["x","y","u","v","V_N","V_T","l"]), np.nan

    # best-fit line y = m x + c
    A = np.c_[x, np.ones(n)]
    m, c0 = np.linalg.lstsq(A, y, rcond=None)[0]

    # unit tangent along line and normal (right-hand)
    t = np.array([1.0, m]); t /= np.linalg.norm(t)
    nvec = np.array([-t[1], t[0]])

    p0 = np.array([x[0], y[0]])
    pts = np.c_[x, y]
    s = (pts - p0) @ t          # signed along-track coordinate

    idx = np.argsort(s)
    s, u, v = s[idx], u[idx], v[idx]

    s_new = np.linspace(s.min(), s.max(), n)
    u_new = np.interp(s_new, s, u)
    v_new = np.interp(s_new, s, v)

    pts_new = p0 + np.outer(s_new, t)
    x_new, y_new = pts_new[:, 0], pts_new[:, 1]

    VT = u_new*t[0] + v_new*t[1]
    VN = u_new*nvec[0] + v_new*nvec[1]

    df = pd.DataFrame({"x": x_new, "y": y_new, "u": u_new, "v": v_new, "V_N": VN, "V_T": VT})
    df["l"] = s_new - s_new.min()   # 0 at start, increasing along transect

    return df, m

# Function to translate solo eddy center approximation to cartesian coordinates
def translate_solo_results(x_l_start, y_l_start, m, l0, r0):
    denom = np.sqrt(1 + m**2)
    x0 = (l0 - r0*m)/denom + x_l_start
    y0 = (l0*m + r0)/denom + y_l_start
    return x0, y0

# DOPPIO

def doppio(x1, y1, u1, v1, x2, y2, u2, v2, plot=False):

    def nan_return():
        nan2 = np.full((2, 2), np.nan)
        coeffs = pd.DataFrame(
            np.nan,
            index=['0', '1', '2', '3'],
            columns=['A', 'B', 'C', 'D']
        )
        return np.nan, np.nan, np.nan, np.full((2,2), np.nan), np.nan

    def clean(*arrs):
        arrs = [np.asarray(a, dtype=float).ravel() for a in arrs]
        m = np.logical_and.reduce([np.isfinite(a) for a in arrs])
        return [a[m] for a in arrs]

    def poly3(p, z):
        return p[0] + p[1]*z + p[2]*z**2 + p[3]*z**3

    def find_root(x, y, degree=3):
        if x.size < degree + 1:
            return np.nan

        try:
            coeffs = np.polyfit(x, y, degree)
            roots = np.roots(coeffs)
        except Exception:
            return np.nan

        real_roots = roots[np.isreal(roots)].real
        if real_roots.size == 0:
            return np.nan

        mid_x = np.median(x)
        return real_roots[np.argmin(np.abs(real_roots - mid_x))]

    x1, y1, u1, v1 = clean(x1, y1, u1, v1)
    x2, y2, u2, v2 = clean(x2, y2, u2, v2)

    if len(x1) < 4 or len(y2) < 4:
        return nan_return()

    pts1 = np.column_stack((x1, y1))
    pts2 = np.column_stack((x2, y2))

    common = np.array([p for p in pts1 if np.any(np.all(pts2 == p, axis=1))])

    if len(common) != 1:
        return nan_return()

    center_x, center_y = common[0]

    # Set origin at transect center
    x = x1 - center_x
    y = y2 - center_y

    # Initial guesses
    # Unconstrained cubic fits used only to estimate x0 and y0
    x0_guess = find_root(x, v1)
    y0_guess = find_root(y, u2)

    pA = np.polyfit(x, u1, 3)[::-1]
    pB = np.polyfit(x, v1, 3)[::-1]
    pC = np.polyfit(y, u2, 3)[::-1]
    pD = np.polyfit(y, v2, 3)[::-1]

    # Free parameters:
    # x0, y0,
    # A0,A1,A2,A3,
    # B1,B2,B3,   with B0 set so V1(x0)=0 exactly
    # C1,C2,C3,   with C0 set so U2(y0)=0 exactly
    # D0,D2,D3,   with D1 set so d1=-a1 exactly
    p0 = np.array([
        x0_guess, y0_guess,
        pA[0], pA[1], pA[2], pA[3],
        pB[1], pB[2], pB[3],
        pC[1], pC[2], pC[3],
        pD[0], pD[2], pD[3]
    ], dtype=float)

    def unpack(p):
        x0, y0 = p[0], p[1]

        A0, A1, A2, A3 = p[2:6]
        B1, B2, B3 = p[6:9]
        C1, C2, C3 = p[9:12]
        D0, D2, D3 = p[12:15]

        a1 = A1 + 2*A2*x0 + 3*A3*x0**2
        b1 = B1 + 2*B2*x0 + 3*B3*x0**2
        c1 = C1 + 2*C2*y0 + 3*C3*y0**2

        B0 = -x0 * b1 # V1(x0)=0 
        C0 = -y0 * c1 # U2(y0)=0 
        D1 = -a1 - 2*D2*y0 - 3*D3*y0**2 # d1=-a1

        A = np.array([A0, A1, A2, A3], dtype=float)
        B = np.array([B0, B1, B2, B3], dtype=float)
        C = np.array([C0, C1, C2, C3], dtype=float)
        D = np.array([D0, D1, D2, D3], dtype=float)

        return x0, y0, A, B, C, D

    def residuals(p):
        x0, y0, A, B, C, D = unpack(p)

        U1 = poly3(A, x)
        V1 = poly3(B, x)
        U2 = poly3(C, y)
        V2 = poly3(D, y)

        return np.r_[U1 - u1, V1 - v1, U2 - u2, V2 - v2]

    res = least_squares(residuals, p0, method='trf')

    if not res.success:
        return nan_return()

    x0, y0, A, B, C, D = unpack(res.x)

    A0, A1, A2, A3 = A
    B0, B1, B2, B3 = B
    C0, C1, C2, C3 = C
    D0, D1, D2, D3 = D

    # Tangent coefficients
    a1 = A1 + 2*A2*x0 + 3*A3*x0**2
    a0 = (A0 + A1*x0 + A2*x0**2 + A3*x0**3) - x0*a1

    b1 = B1 + 2*B2*x0 + 3*B3*x0**2
    b0 = -x0*b1

    c1 = C1 + 2*C2*y0 + 3*C3*y0**2
    c0 = -y0*c1

    d1 = D1 + 2*D2*y0 + 3*D3*y0**2
    d0 = (D0 + D1*y0 + D2*y0**2 + D3*y0**3) - y0*d1

    # Enforce d1 = -a1 numerically if tiny mismatch remains
    d1 = -a1

    radicand = -b1*c1 - a1**2
    if radicand <= 0 or not np.isfinite(radicand):
        return nan_return()

    Omega = np.sign(b1) * np.sqrt(radicand)
    if Omega == 0 or not np.isfinite(Omega):
        return nan_return()

    Q = (1 / Omega) * np.array([
        [b1,  -a1],
        [-a1, -c1]
    ], dtype=float)

    xc = (a0*a1 + c1*d0) / Omega**2
    yc = (a0*b1 - a1*d0) / Omega**2

    w = Omega * (Q[0,0] + Q[1,1])

    coeffs = pd.DataFrame(
        {
            'A': A,
            'B': B,
            'C': C,
            'D': D
        },
        index=['0', '1', '2', '3']
    )

    yu1 = a0 + a1*x
    yv1 = b0 + b1*x
    yu2 = c0 + c1*y
    yv2 = d0 + d1*y

    if plot:
        fig, axs = plt.subplots(2, 2, figsize=(8, 6))

        xs = np.linspace(x.min(), x.max(), 400)
        ys = np.linspace(y.min(), y.max(), 400)

        U1s = poly3(A, xs)
        V1s = poly3(B, xs)
        U2s = poly3(C, ys)
        V2s = poly3(D, ys)

        axs[0, 0].scatter(x, u1, s=10)
        axs[0, 0].plot(xs, U1s)
        axs[0, 0].plot(x, yu1)
        axs[0, 0].axvline(x0, ls='--')
        axs[0, 0].set_title("U1(x)")
        axs[0, 0].set_xlabel("x (m)")
        axs[0, 0].set_ylabel(r"(ms$^{-1}$)")

        axs[0, 1].scatter(x, v1, s=10)
        axs[0, 1].plot(xs, V1s)
        axs[0, 1].plot(x, yv1)
        axs[0, 1].axvline(x0, ls='--')
        axs[0, 1].axhline(0, ls='--')
        axs[0, 1].set_title("V1(x)")
        axs[0, 1].set_xlabel("x (m)")
        axs[0, 1].set_ylabel(r"(ms$^{-1}$)")

        axs[1, 0].scatter(y, u2, s=10)
        axs[1, 0].plot(ys, U2s)
        axs[1, 0].plot(y, yu2)
        axs[1, 0].axvline(y0, ls='--')
        axs[1, 0].axhline(0, ls='--')
        axs[1, 0].set_title("U2(y)")
        axs[1, 0].set_xlabel("y (m)")
        axs[1, 0].set_ylabel(r"(ms$^{-1}$)")

        axs[1, 1].scatter(y, v2, s=10)
        axs[1, 1].plot(ys, V2s)
        axs[1, 1].plot(y, yv2)
        axs[1, 1].axvline(y0, ls='--')
        axs[1, 1].set_title("V2(y)")
        axs[1, 1].set_xlabel("y (m)")
        axs[1, 1].set_ylabel(r"(ms$^{-1}$)")

        plt.tight_layout()
        plt.show()

    return xc + center_x, yc + center_y, w, Q, Omega

# LATTE

def latte(xi, yi, ui, vi):
    from scipy.optimize import least_squares
    xi, yi, ui, vi = map(lambda a: np.asarray(a, float), (xi, yi, ui, vi))
    m = np.isfinite(xi) & np.isfinite(yi) & np.isfinite(ui) & np.isfinite(vi)
    x, y, u_i, v_i = xi[m], yi[m], ui[m], vi[m]
    n = x.size
    if n < 3:
        return np.nan, np.nan, np.nan, np.full((2,2), np.nan), np.nan, np.nan

    p0 = np.array([x.mean(), y.mean(), 1., 0., 1.])  # xc,yc,Oq11,Oq12,Oq22

    def fun(p):
        xc,yc,Oq11,Oq12,Oq22 = p
        dx, dy = x-xc, y-yc
        u = -Oq22*dy - Oq12*dx
        v =  Oq11*dx + Oq12*dy
        r = np.empty(2*n); r[:n]=u-u_i; r[n:]=v-v_i
        return r

    def jac(p):
        xc,yc,Oq11,Oq12,Oq22 = p
        dx, dy = x-xc, y-yc
        J = np.zeros((2*n,5))
        # u = -Omega*q22*(y-yc) - Omega*q12*(x-xc)
        J[:n,0], J[:n,1], J[:n,3], J[:n,4] = Oq12, Oq22, -dx, -dy
        # v = Omega*q11*(x-xc) + Omega*q12*(y-yc)
        J[n:,0], J[n:,1], J[n:,2], J[n:,3] = -Oq11, -Oq12, dx, dy
        return J

    xc,yc,Oq11,Oq12,Oq22 = least_squares(fun, p0, jac=jac).x

    dx, dy = x-xc, y-yc
    u = -Oq22*dy - Oq12*dx
    v =  Oq11*dx + Oq12*dy
    err2 = ((u-u_i)**2 + (v-v_i)**2).sum()
    tot2 = ((u_i-u_i.mean())**2 + (v_i-v_i.mean())**2).sum()
    r2 = 1 - err2/tot2 if tot2 > 0 else np.nan

    w = Oq11+Oq22
    OQ = np.array([[Oq11,Oq12],[Oq12,Oq22]])
    det = OQ[0,0]*OQ[1,1] - OQ[0,1]*OQ[1,0]
    Omega = np.sign(Oq11)*np.sqrt(abs(det))
    Q = OQ/Omega if Omega != 0 else np.full((2,2), np.nan)

    return xc, yc, w, Q, Omega, r2


# Outer-core ESP paramter finder
    
def out_core_param_fit(
    rho2, Qr, vt,
    Omega0=None, Rc0=None,
    plot=False, ax=None,
    maxfev=10000, Rc_max=1e5,
    r2_flag=False,
    rho_plot_max=None, n_curve=400,
    km_flag=False,
    ci_flag=False,
    pred_flag=False,
):

    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt

    rho2 = np.asarray(rho2, float)
    Qr   = np.asarray(Qr, float)
    vt   = np.asarray(vt, float)

    m = np.isfinite(rho2) & np.isfinite(Qr) & np.isfinite(vt) & (rho2 >= 0) & (Qr != 0)
    if not np.any(m):
        return (np.nan, np.nan, np.nan, np.nan) if r2_flag else (np.nan, np.nan, np.nan)

    rho2 = rho2[m]
    Qr   = Qr[m]
    vt   = vt[m]

    rho = np.sqrt(rho2)
    vt = vt * (rho / Qr)

    def vt_model(r2, Omega, Rc):
        return Omega * np.sqrt(r2) * np.exp(-r2 / (Rc**2))

    i = np.nanargmax(np.abs(vt))
    rho_max = rho[i]

    if Rc0 is None:
        Rc0 = max(rho_max * np.sqrt(2), 1e-6)

    if Omega0 is None:
        denom = rho * np.exp(-rho2 / (Rc0**2))
        ok = np.abs(denom) > 0
        Omega0 = np.nanmedian(vt[ok] / denom[ok]) if np.any(ok) else 0

    if not np.isfinite(Omega0):
        Omega0 = 0

    pcov = None
    try:
        popt, pcov = curve_fit(
            vt_model, rho2, vt,
            p0=[Omega0, Rc0],
            bounds=([-np.inf, 1e-8], [np.inf, np.inf]),
            maxfev=maxfev
        )
        Omega_opt, Rc_opt = popt
    except:
        Omega_opt, Rc_opt = Omega0, Rc0

    if (not np.isfinite(Rc_opt)) or (Rc_opt > Rc_max):
        Omega_opt, Rc_opt = Omega0, Rc0
        pcov = None

    psi0_opt = -0.5 * Omega_opt * Rc_opt**2

    vt_fit = vt_model(rho2, Omega_opt, Rc_opt)

    ss_res = np.sum((vt - vt_fit)**2)
    ss_tot = np.sum((vt - vt.mean())**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    dof = max(len(vt) - 2, 1)
    sigma2 = ss_res / dof

    if plot:

        if ax is None:
            fig, ax = plt.subplots()

        if rho_plot_max is None:
            rho_plot_max = np.nanmax(rho)

        r_grid = np.linspace(0, rho_plot_max, n_curve)
        r2_grid = r_grid**2

        vt_grid = vt_model(r2_grid, np.abs(Omega_opt), Rc_opt)

        if km_flag:
            core_mask = rho <= 30
        else:
            core_mask = rho <= 30_000

        ax.scatter(rho[core_mask], np.abs(vt[core_mask]), s=10, color='m', label='Inner-core \nobserved')
        ax.scatter(rho[~core_mask], np.abs(vt[~core_mask]), s=10, color='g', label='Outer-core \nobserved')

        ax.plot(r_grid, np.abs(vt_grid), lw=2, color='b', label='')
        ax.axvline(Rc_opt / np.sqrt(2), ls='--', color='r', label='', lw=2)

        if pcov is not None:

            exp_term = np.exp(-r2_grid / (Rc_opt**2))

            dOmega = np.sqrt(r2_grid) * exp_term
            dRc = 2 * Omega_opt * (r2_grid**1.5) * exp_term / (Rc_opt**3)

            J = np.vstack([dOmega, dRc]).T
            var_model = np.einsum("ij,jk,ik->i", J, pcov, J)

            if ci_flag:
                se_model = np.sqrt(np.maximum(var_model, 0))
                lo = vt_grid - 1.96 * se_model
                hi = vt_grid + 1.96 * se_model

                ax.fill_between(
                    r_grid, np.abs(lo), np.abs(hi),
                    color='orange', alpha=.2, label='95% CI'
                )

            if pred_flag:
                se_pred = np.sqrt(np.maximum(var_model + sigma2, 0))
                lo = vt_grid - 1.96 * se_pred
                hi = vt_grid + 1.96 * se_pred

                ax.fill_between(
                    r_grid, lo, hi,
                    color='b', alpha=.15, label=''
                )

        ax.set_xlabel(r'$\rho$')
        ax.set_ylabel(r'$|v_t^\star|$')

        ax.set_title(
            f"Omega={Omega_opt:.3g}, Rc={Rc_opt:.3g}, psi0={psi0_opt:.3g}, R²={R2:.2f}"
        )
        ax.set_ylim(0, None); ax.set_xlim(0, None)

        ax.legend(loc='upper left')

    return (Rc_opt, psi0_opt, Omega_opt, R2) if r2_flag else (Rc_opt, psi0_opt, Omega_opt)


# Helper functions

def doppio_pipeliner(nxc, nyc, ut, vt, X_new, Y_new, r=30000.0):
    nan = (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    x = np.asarray(X_new[:, 0], float)
    y = np.asarray(Y_new[0, :], float)

    dx = np.nanmedian(np.abs(np.diff(x)))
    dy = np.nanmedian(np.abs(np.diff(y)))
    cell = np.nanmax([dx, dy])
    if not np.isfinite(cell) or cell == 0:
        return nan

    margin = int(np.ceil(r / cell))

    ic = int(np.clip(np.searchsorted(x, nxc), 1, x.size-1))
    ic -= (nxc - x[ic-1] < x[ic] - nxc)
    jc = int(np.clip(np.searchsorted(y, nyc), 1, y.size-1))
    jc -= (nyc - y[jc-1] < y[jc] - nyc)

    if ic < margin or ic >= x.size - margin or jc < margin or jc >= y.size - margin:
        return nan

    # x-transect (y = y[jc])
    i0 = np.searchsorted(x, nxc - r, side="left")
    i1 = np.searchsorted(x, nxc + r, side="right")
    x1 = x[i0:i1]
    y1 = np.full(x1.size, y[jc])
    u1 = ut[i0:i1, jc]
    v1 = vt[i0:i1, jc]

    # y-transect (x = x[ic])
    j0 = np.searchsorted(y, nyc - r, side="left")
    j1 = np.searchsorted(y, nyc + r, side="right")
    y2 = y[j0:j1]
    x2 = np.full(y2.size, x[ic])
    u2 = ut[ic, j0:j1]
    v2 = vt[ic, j0:j1]

    return x1, y1, u1, v1, x2, y2, u2, v2



def latte_source_selector(ds_sadcp, ds_sat, source="multi", z_target=37.0, plot=False):
    xc_pre, yc_pre = 0.0, 0.0
    q11, q12, q22 = 1.0, 0.0, 1.0
    rho_core, rho_outer = 35_000.0, 100_000.0

    nanQ = np.full((2, 2), np.nan)

    def empty_return():
        df_row = pd.DataFrame([{
            "source": source, "xc": np.nan, "yc": np.nan, "w": np.nan, "Q": nanQ,
            "Omega": np.nan, "Rc": np.nan, "psi0": np.nan,
            "time": time, "alpha": np.nan
        }])
        df_xyuv = pd.DataFrame(columns=["source", "core", "xi", "yi", "ui", "vi"])
        df_fit = pd.DataFrame(columns=["rho2", "Qr", "vt"])
        return df_row, df_xyuv, df_fit

    t_mean = ds_sadcp.time.mean().values
    time = t_mean

    sat_interp = ds_sat.interp(time=t_mean)
    X, Y = ds_sat.x.values, ds_sat.y.values
    xg, yg = X.ravel(), Y.ravel()
    ut, vt = sat_interp.u.values.ravel(), sat_interp.v.values.ravel()

    sat_mask = np.isfinite(xg) & np.isfinite(yg) & np.isfinite(ut) & np.isfinite(vt)
    xg, yg, ut, vt = xg[sat_mask], yg[sat_mask], ut[sat_mask], vt[sat_mask]

    z0 = float(ds_sadcp.z.sel(z=z_target, method="nearest").values)
    xd, yd = ds_sadcp.x.values, ds_sadcp.y.values
    ud = ds_sadcp.u.sel(z=z0).values
    vd = ds_sadcp.v.sel(z=z0).values

    sadcp_mask = np.isfinite(xd) & np.isfinite(yd) & np.isfinite(ud) & np.isfinite(vd)
    xd, yd, ud, vd = xd[sadcp_mask], yd[sadcp_mask], ud[sadcp_mask], vd[sadcp_mask]

    if source == "multi":
        xi = np.concatenate([xd, xg])
        yi = np.concatenate([yd, yg])
        ui = np.concatenate([ud, ut])
        vi = np.concatenate([vd, vt])
    elif source == "sadcp":
        xi, yi, ui, vi = xd, yd, ud, vd
    elif source == "sat":
        xi, yi, ui, vi = xg, yg, ut, vt
    else:
        return empty_return()

    dx, dy = xi - xc_pre, yi - yc_pre
    rho2 = q11 * dx**2 + 2 * q12 * dx * dy + q22 * dy**2
    core1 = rho2 <= rho_core**2

    if not np.any(core1):
        return empty_return()

    xc1, yc1, w1, Q1, Omega1, _ = latte(xi[core1], yi[core1], ui[core1], vi[core1])

    q11, q12, q22 = Q1[0, 0], Q1[0, 1], Q1[1, 1]
    dx, dy = xi - xc1, yi - yc1
    rho2 = q11 * dx**2 + 2 * q12 * dx * dy + q22 * dy**2
    core2 = rho2 <= rho_core**2

    if not np.any(core2):
        return empty_return()

    xc, yc, w, Q, Omega, _ = latte(xi[core2], yi[core2], ui[core2], vi[core2])
    alpha = axis_ratio_from_Q(Q)

    df_inner = pd.DataFrame({
        "source": source,
        "core": "inner",
        "xi": xi[core2],
        "yi": yi[core2],
        "ui": ui[core2],
        "vi": vi[core2]
    })

    q11, q12, q22 = Q[0, 0], Q[0, 1], Q[1, 1]
    dx, dy = xi - xc, yi - yc
    rho2 = q11 * dx**2 + 2 * q12 * dx * dy + q22 * dy**2
    outer = rho2 <= rho_outer**2
    outer_only = outer & ~core2

    df_outer = pd.DataFrame({
        "source": source,
        "core": "outer",
        "xi": xi[outer_only],
        "yi": yi[outer_only],
        "ui": ui[outer_only],
        "vi": vi[outer_only]
    })

    xi_o, yi_o = xi[outer], yi[outer]
    ui_o, vi_o = ui[outer], vi[outer]
    dx_o, dy_o = dx[outer], dy[outer]
    rho2_o = rho2[outer]

    vt_o = tangential_velocity(xi_o, yi_o, ui_o, vi_o, xc, yc, Q)
    Qr_o = np.sqrt((q11 * dx_o + q12 * dy_o)**2 + (q12 * dx_o + q22 * dy_o)**2)

    sign_mask = (vt_o <= 0) if (Omega < 0) else (vt_o >= 0)
    rho2_f = rho2_o[sign_mask]
    Qr_f = Qr_o[sign_mask]
    vt_f = vt_o[sign_mask]

    if len(rho2_f) == 0:
        return empty_return()

    Rc, psi0, Omega_opt = out_core_param_fit(
        rho2_f, Qr_f, vt_f, Omega0=Omega, plot=plot, pred_flag=True
    )
    w = Omega_opt * (q11 + q22)

    df_row = pd.DataFrame([{
        "source": source,
        "xc": xc,
        "yc": yc,
        "w": w,
        "Q": Q,
        "Omega": Omega_opt,
        "Rc": Rc,
        "psi0": psi0,
        "time": time,
        "alpha": alpha
    }])

    df_fit = pd.DataFrame({
        "rho2": rho2_f,
        "Qr": Qr_f,
        "vt": vt_f
    })

    return df_row, pd.concat([df_inner, df_outer], ignore_index=True), df_fit

def axis_ratio_from_Q(Q):
    lam = np.abs(np.linalg.eigvalsh(Q))
    return np.sqrt(lam.max() / lam.min())

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

def model_uv_at_xy(xi, yi, xc, yc, Q, Omega, Rc):
    dx = xi - xc
    dy = yi - yc

    q11, q12, q22 = Q[0,0], Q[0,1], Q[1,1]

    rho2 = q11*dx*dx + 2*q12*dx*dy + q22*dy*dy
    fac = Omega * np.exp(-rho2 / (Rc*Rc))

    uhat = -fac * (q12*dx + q22*dy)
    vhat =  fac * (q11*dx + q12*dy)

    return uhat, vhat
    
def vector_R2(u, v, uhat, vhat):
    m = np.isfinite(u) & np.isfinite(v) & np.isfinite(uhat) & np.isfinite(vhat)
    if not np.any(m):
        return np.nan
    u, v, uhat, vhat = u[m], v[m], uhat[m], vhat[m]
    err2 = (uhat - u)**2 + (vhat - v)**2
    u0, v0 = np.mean(u), np.mean(v)
    tot2 = (u - u0)**2 + (v - v0)**2
    return 1 - np.sum(err2)/np.sum(tot2) if np.sum(tot2) > 0 else np.nan









    