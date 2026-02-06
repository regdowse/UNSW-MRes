import numpy as np
import pandas as pd

# import sys
# sys.path.append("/home/z5297792/UNSW-MRes/MRes/modules") 
# from utils import 

def calculate_eddy(width=500e3, num_depth_layers=21, eta0=1, L=1e5, H=1000,
                   a=1, b=1, rho0=1025, f0=None,
                   alpha_1=0.01, alpha_2=0.01,
                   T0=20, dTdz=0.005, dSdz=0.01, taper_depth=None,
                   q11=1.0, q12=0.0, q22=1.0):
    """
    Compute 3D geostrophic U,V and fields for an eddy with prescribed
    quadratic-form shape coefficients q11,q12,q22 (constant in z).
    """
    if f0 is None:
        f0 = 2 * 7.29e-5 * np.sin(np.radians(-34))
    g = 9.81

    x = np.linspace(-width/2, width/2, 51)
    y = np.linspace(-width/2, width/2, 51)
    z = np.linspace(-H, 0, num_depth_layers)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    taper = np.exp(Z / taper_depth) if taper_depth is not None else 1.0

    # center shift with depth
    Xc = alpha_1 * Z
    Yc = alpha_2 * Z

    # elliptical radius-squared
    dX = X - Xc
    dY = Y - Yc
    r2Q = q11 * dX**2 + 2*q12 * dX * dY + q22 * dY**2

    phi = np.exp(-r2Q / L**2) * taper

    T = -2 * phi
    S = 1.5 * phi
    P = -rho0 * g * Z

    alpha = 2e-4
    beta  = 8e-4
    sigma = rho0 * (1 - alpha*(T - T0) + beta*(S - 35))
    sigma -= sigma.mean()

    dx, dy = x[1]-x[0], y[1]-y[0]
    dsdx = np.gradient(sigma, dx, axis=0, edge_order=2)
    dsdy = np.gradient(sigma, dy, axis=1, edge_order=2)

    U =  g/f0 * dsdy * taper
    V = -g/f0 * dsdx * taper

    # flip so z=0 is first index
    U, V, sigma, T, S, P = [np.flip(arr, axis=2) for arr in (U,V,sigma,T,S,P)]
    z = z[::-1]

    # transpose horizontal slices
    for k in range(U.shape[2]):
        U[:,:,k] = U[:,:,k].T
        V[:,:,k] = V[:,:,k].T

    return U, V, sigma, T, S, P, x/1000, y/1000, z/1000

def calculate_eddy_2D(width=500000, L=1e5, f0=None, rho0=1025, q11=1., q22=1., q12=0.):
    if f0 is None:
        f0 = 2 * 7.29E-5 * np.sin(np.radians(-34))
    g = 9.81
    x = np.linspace(-width // 2, width // 2, 101)
    y = np.linspace(-width // 2, width // 2, 101)
    x_2d, y_2d = np.meshgrid(x, y, indexing='ij')
    
    x_c = 0.0
    y_c = 0.0

    X = np.stack([x_2d - x_c, y_2d - y_c], axis=0) 
    Q = np.array([[q11, q12], [q12, q22]]) 
    
    r_c = np.sqrt(np.einsum('i...,ij,j...->...', X, Q, X))

    sigma = -2 * np.exp(-r_c**2 / L**2)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    d_sigma_dx = np.gradient(sigma, axis=0) / dx
    d_sigma_dy = np.gradient(sigma, axis=1) / dy
    
    U = -g / f0 * d_sigma_dy
    V = g / f0 * d_sigma_dx
    
    return U, V, x/1000, y/1000

def plot_ellipse(Q, center=(0, 0), scale=1):
    def normalize_matrix(A, norm_type='fro'):
        norm = np.linalg.norm(A, 'fro') if norm_type == 'fro' else np.max(np.abs(A))
        return A / norm if norm else A
    Q = normalize_matrix(Q)

    def swap_principal_axes(Q):
        eigvals, eigvecs = np.linalg.eigh(Q)
        return eigvecs @ np.diag(eigvals[::-1]) @ eigvecs.T

    Q = swap_principal_axes(Q)
    
    eigenvalues, eigenvectors = np.linalg.eigh(Q)
    if np.any(eigenvalues < 0):
        Q = np.array([[np.abs(Q[0,0]), Q[0,1]], [Q[1,0], np.abs(Q[1,1])]])

        def flip_Q_y(Q):
            F_y = np.diag([-1, 1])
            return F_y.T @ Q @ F_y

        Q = flip_Q_y(Q)

        eigenvalues, eigenvectors = np.linalg.eigh(Q)
        if np.any(eigenvalues < 0):
            return np.nan, np.nan
            
    a, b = np.sqrt(eigenvalues) * scale
    theta = np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])
    t = np.linspace(0, 2 * np.pi, 100)
    x, y = a * np.cos(t), b * np.sin(t)
    R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    x_ellipse, y_ellipse = R @ np.array([x, y]) + np.array(center).reshape(2, 1)
    return x_ellipse, y_ellipse













    

################################################ ESP Methods ################################################
def moca(l, VT, VN, Rc_max=1e5, plot_flag=False, df_flag=False):

    l, VT, VN = [np.asarray(a) for a in (l, VT, VN)]
    mask = ~np.isnan(l) & ~np.isnan(VT) & ~np.isnan(VT)
    l, VT, VN = l[mask], VT[mask], VN[mask]

    if len(l) == 0:
        nan2 = np.array([[np.nan, np.nan], [np.nan, np.nan]])
        return (
            np.nan,
            np.nan,  
            np.nan,  
            nan2,  
            np.nan,  
            np.nan,
            nan2
        )

    
    def find_root(x, y):
        coeffs = np.polyfit(x, y, 3)
        roots = np.roots(np.poly1d(coeffs))
        real_roots = roots[np.isreal(roots)].real
        mid = x[len(x)//2]
        return real_roots[np.argmin(np.abs(mid - real_roots))]
    
    def tang_at_root(x, y, rx):
        coeffs = np.polyfit(x, y, 3)
        deriv = np.polyder(coeffs)
        slope = np.polyval(deriv, rx)
        intercept = np.polyval(coeffs, rx) - slope * rx
        return slope, intercept
    
    def cubic_interpolate(x, y, root):
        coeffs = np.polyfit(x, y, 3)
        return np.polyval(coeffs, root)
    
    root = find_root(l, VN)
    c, b = tang_at_root(l, VN, root)  # c: slope, b: intercept
    a = cubic_interpolate(l, VT, root)
    
    l0 = -b / c
    r0 = a / c
    w = 2 * c
    
    xc, yc = l0, r0

    Aq11, Aq12, Aq22 = w/4, 0.0, w/4

    AQ = np.array([[Aq11, Aq12], [Aq12, Aq22]])
    detAQ = np.linalg.det(AQ)
    A = (abs(detAQ))**(0.5)
    A = np.sign(Aq11)*A
    Q = AQ / A
    q11, q12, q22 = Q[0,0], 0.0, Q[1,1]

    xi, yi, ui, vi = l, [0]*len(l), VT, VN
    
    # fit Rc, psi0, A
    df = psi_params(xc, yc, Q, xi, yi, ui, vi) # input data into m
    if A < 0:
        mask = df.vt <= 0
    else:
        mask = df.vt >= 0
    rho2, Qr, vt = df.rho2[mask], df.Qr[mask], df.vt[mask]

    Rc_opt, psi0_opt, A_opt = fit_psi_params(rho2, Qr, vt, A0=A,
                                             plot=plot_flag, Rc_max=Rc_max)

    if df_flag:
        return xc, yc, w, Q, Rc_opt, psi0_opt, A_opt, df
    return xc, yc, w, Q, Rc_opt, psi0_opt, A_opt

def obs_moca(l, VT, VN, Rc_max=1e5, plot_flag=False, df_flag=False):

    l, VT, VN = [np.asarray(a) for a in (l, VT, VN)]
    mask = ~np.isnan(l) & ~np.isnan(VT) & ~np.isnan(VT)
    l_all, VT_all, VN_all = l[mask], VT[mask], VN[mask]

    if len(l) == 0:
        nan2 = np.array([[np.nan, np.nan], [np.nan, np.nan]])
        return (
            np.nan,
            np.nan,  
            np.nan,  
            nan2,  
            np.nan,  
            np.nan,
            nan2
        )
    def find_root(x, y):
        coeffs = np.polyfit(x, y, 3)
        roots = np.roots(np.poly1d(coeffs))
        real_roots = roots[np.isreal(roots)].real
        mid = x[len(x)//2]
        return real_roots[np.argmin(np.abs(mid - real_roots))]
    def tang_at_root(x, y, rx):
        coeffs = np.polyfit(x, y, 3)
        deriv = np.polyder(coeffs)
        slope = np.polyval(deriv, rx)
        intercept = np.polyval(coeffs, rx) - slope * rx
        return slope, intercept
    def cubic_interpolate(x, y, root):
        coeffs = np.polyfit(x, y, 3)
        return np.polyval(coeffs, root)

    # Find point of closest approach
    root = find_root(l, VN)
    c, b = tang_at_root(l, VN, root)  # c: slope, b: intercept
    a = cubic_interpolate(l, VT, root)

    # Use core data to find xc, yc
    mask = np.abs(l-root) < 30e3 #m
    l, VT, VN = l_all[mask], VT_all[mask], VN_all[mask]
    
    root = find_root(l, VN)
    c, b = tang_at_root(l, VN, root)  # c: slope, b: intercept
    a = cubic_interpolate(l, VT, root)
    
    l0 = -b / c
    r0 = a / c
    w = 2 * c
    
    xc, yc = l0, r0

    Aq11, Aq12, Aq22 = w/4, 0.0, w/4

    AQ = np.array([[Aq11, Aq12], [Aq12, Aq22]])
    detAQ = np.linalg.det(AQ)
    A = (abs(detAQ))**(0.5)
    A = np.sign(Aq11)*A
    Q = AQ / A
    q11, q12, q22 = Q[0,0], 0.0, Q[1,1]

    xi, yi, ui, vi = l_all, [0]*len(l_all), VT_all, VN_all
    
    # fit Rc, psi0, A
    df = psi_params(xc, yc, Q, xi, yi, ui, vi) # input data into m
    if A < 0:
        mask = df.vt <= 0
    else:
        mask = df.vt >= 0
    rho2, Qr, vt = df.rho2[mask], df.Qr[mask], df.vt[mask]

    Rc_opt, psi0_opt, A_opt = fit_psi_params(rho2, Qr, vt, A0=A,
                                             plot=plot_flag, Rc_max=Rc_max)

    if df_flag:
        return xc, yc, w, Q, Rc_opt, psi0_opt, A_opt, df
    return xc, yc, w, Q, Rc_opt, psi0_opt, A_opt

def dopioe(x1, y1, u1, v1, x2, y2, u2, v2, Rc_max=1e6, plot_flag=False):

    x1, y1, u1, v1 = [np.asarray(a) for a in (x1, y1, u1, v1)]
    mask = ~np.isnan(x1) & ~np.isnan(y1) & ~np.isnan(u1) & ~np.isnan(v1)
    x1, y1, u1, v1 = x1[mask], y1[mask], u1[mask], v1[mask]

    x2, y2, u2, v2 = [np.asarray(a) for a in (x2, y2, u2, v2)]
    mask = ~np.isnan(x2) & ~np.isnan(y2) & ~np.isnan(u2) & ~np.isnan(v2)
    x2, y2, u2, v2 = x2[mask], y2[mask], u2[mask], v2[mask]
    
    if (len(x1) == 0) &  (len(x2) == 0):
        nan2 = np.array([[np.nan, np.nan], [np.nan, np.nan]])
        return (
            np.nan,
            np.nan,  
            np.nan,  
            nan2,  
            np.nan,  
            np.nan,
            nan2
        )

    def find_root(x, y, degree=3):
        """Fit a degree-3 polynomial to (x, y) and return the real root closest to x's midpoint."""
        coeffs = np.polyfit(x, y, degree)
        roots = np.roots(np.poly1d(coeffs))
        real_roots = roots[np.isreal(roots)].real
        mid_x = x[len(x) // 2]
        return real_roots[np.argmin(np.abs(real_roots - mid_x))]
    
    def tangent_at_root(x, y, root_x, degree=3):
        coeffs = np.polyfit(x, y, degree)
        slope = np.polyval(np.polyder(coeffs), root_x)
        intercept = np.polyval(coeffs, root_x) - slope * root_x
        y_fit = np.polyval(coeffs, x)
        r2 = 1 - np.sum((y - y_fit)**2) / np.sum((y - np.mean(y))**2) if np.sum((y - np.mean(y))**2) != 0 else 1.0
        return slope, intercept, r2

    points1 = set(zip(x1, y1))
    points2 = set(zip(x2, y2))
    common_points = points1 & points2
    
    if len(common_points) != 1:
        print(f"Warning: Expected 1 common point, found {len(common_points)}.")
        nan2 = np.array([[np.nan, np.nan], [np.nan, np.nan]])
        return (
            np.nan,
            np.nan,  
            np.nan,  
            nan2,  
            np.nan,  
            np.nan,
            nan2
        )

    center_x, center_y = next(iter(common_points))

    xx = x1 - center_x
    yy = y2 - center_y

    # Compute roots and tangents
    root_x = find_root(xx, v1)
    A1, A0, r2A = tangent_at_root(xx, u1, root_x)
    B1, B0, r2B = tangent_at_root(xx, v1, root_x)
    root_y = find_root(yy, u2)
    C1, C0, r2C = tangent_at_root(yy, u2, root_y)
    D1, D0, r2D = tangent_at_root(yy, v2, root_y)
    
    alpha = A0 if r2A > r2C else C0
    beta  = B0 if r2B > r2D else D0
    gamma = A1 if r2A > r2D else -D1
    
    Aq11 = B1 / 2
    Aq22 = -C1 / 2
    Aq12 = -gamma / 2
    denom = C1 * B1 + gamma**2
    if denom == 0:
        raise ZeroDivisionError("Denominator is zero.")
        
    xc = - (alpha * gamma + beta * C1) / denom  + center_x
    yc = (beta * gamma - alpha * B1) / denom + center_y
    w = 2 * (Aq11 + Aq22)

    AQ = np.array([[Aq11, Aq12], [Aq12, Aq22]])
    detAQ = np.linalg.det(AQ)
    A = (abs(detAQ))**(0.5)
    A = np.sign(Aq11)*A
    Q = AQ / A
    q11, q12, q22 = Q[0,0], Q[1,0], Q[1,1]

    # Remove duplicates (i.e., the center) from (x1, y1, u1, v1)
    mask = ~np.array([(x, y) in common_points for x, y in zip(x1, y1)])
    x1f = x1[mask]
    y1f = y1[mask]
    u1f = u1[mask]
    v1f = v1[mask]
    # Concatenate with unaltered second set
    xi = np.concatenate([x1f, x2])
    yi = np.concatenate([y1f, y2])
    ui = np.concatenate([u1f, u2])
    vi = np.concatenate([v1f, v2])

    # fit Rc, psi0, A
    df = psi_params(xc, yc, Q, xi, yi, ui, vi) # input data into m
    if A < 0:
        mask = df.vt <= 0
    else:
        mask = df.vt >= 0
    rho2, Qr, vt = df.rho2[mask], df.Qr[mask], df.vt[mask]

    Rc_opt, psi0_opt, A_opt = fit_psi_params(rho2, Qr, vt, A0=A,
                                             plot=plot_flag, Rc_max=Rc_max)
    return xc, yc, w, Q, Rc_opt, psi0_opt, A_opt


def espra(xi, yi, ui, vi, Rc_max=1e5, plot_flag=False, ax=None, r2_flag=False):
    xi, yi, ui, vi = [np.asarray(a) for a in (xi, yi, ui, vi)]
    mask = ~np.isnan(xi) & ~np.isnan(yi) & ~np.isnan(ui) & ~np.isnan(vi)
    xi, yi, ui, vi = xi[mask], yi[mask], ui[mask], vi[mask]
    
    from scipy.optimize import least_squares

    def residuals(params, x, y, u_i, v_i):
        xc, yc, Aq11, Aq12, Aq22 = params
        u = -2 * Aq22 * (y - yc) - 2 * Aq12 * (x - xc)
        v =  2 * Aq11 * (x - xc) + 2 * Aq12 * (y - yc)
        return np.concatenate([(u - u_i), (v - v_i)])

    def fit_params(x, y, u_i, v_i):
        xc_init, yc_init = np.mean(x), np.mean(y)
        Aq11_init, Aq12_init, Aq22_init = 1.0, 0.0, 1.0 
        params_init = [xc_init, yc_init, Aq11_init, Aq12_init, Aq22_init]
        result = least_squares(residuals, params_init, args=(x, y, u_i, v_i))
        return result.x 

    xc, yc, Aq11, Aq12, Aq22 = fit_params(xi, yi, ui, vi)

    # Predicted velocities
    u_pred = -2 * Aq22 * (yi - yc) - 2 * Aq12 * (xi - xc)
    v_pred =  2 * Aq11 * (xi - xc) + 2 * Aq12 * (yi - yc)
    # --- Vector R² computation ---
    err2 = (u_pred - ui)**2 + (v_pred - vi)**2
    mean_u, mean_v = np.mean(ui), np.mean(vi)
    tot2 = (ui - mean_u)**2 + (vi - mean_v)**2
    r2_core = 1 - np.sum(err2) / np.sum(tot2) if np.sum(tot2) != 0 else np.nan

    w = 2*(Aq11 + Aq22)

    AQ = np.array([[Aq11, Aq12], [Aq12, Aq22]])
    detAQ = np.linalg.det(AQ)
    A = (abs(detAQ))**(0.5)
    A = np.sign(Aq11)*A
    Q = AQ / A

    dx, dy = xi - xc, yi - yc
    rho2 = Q[0,0]*dx**2 + 2*Q[1,0]*dx*dy + Q[1,1]*dy**2
    Qr = np.sqrt((Q[0,0]*dx + Q[1,0]*dy)**2 + (Q[1,0]*dx + Q[1,1]*dy)**2)
    vt = tangential_velocity(xi, yi, ui, vi, xc, yc, Q)
    if A < 0:
        mask = vt <= 0
    else:
        mask = vt >= 0
    rho2, Qr, vt = rho2[mask], Qr[mask], vt[mask]

    Rc_opt, psi0_opt, A_opt, r2_outer_core = fit_psi_params(rho2, Qr, vt, A0=A,
                                                            plot=plot_flag, Rc_max=Rc_max, ax=ax, r2_flag=True)
    if r2_flag:
        return xc, yc, w, Q, Rc_opt, psi0_opt, A_opt, r2_core, r2_outer_core
    return xc, yc, w, Q, Rc_opt, psi0_opt, A_opt

def dopioe_pipeliner(nxc, nyc, cyc, ut, vt, X_new, Y_new, r=30000):

    R_grid = np.hypot(nxc - X_new, nyc - Y_new)
    ic, jc = map(int, np.unravel_index(np.argmin(R_grid), R_grid.shape))

    # DOPIOE wont work if too close to boundary
    x_new = X_new[:, 0]
    y_new = Y_new[0, :]
    dx = np.max(np.diff(x_new))  # spacing in x-direction
    dy = np.max(np.diff(y_new))  # spacing in y-direction
    cell_size = np.max([dx, dy])        # average cell size in Euclidean units
    margin = int(np.ceil(r / cell_size)) 

    if (ic < margin or ic >= X_new.shape[0] - margin or
            jc < margin or jc >= X_new.shape[1] - margin):
        return np.nan, np.nan, np.nan, np.array([[np.nan, np.nan],
                                                [np.nan, np.nan]]), np.nan, np.nan, np.nan, np.nan, np.nan

    # horizontal transect (constant y = y[jc])
    x_mask = np.abs(x_new - nxc) <= r
    x1 = x_new[x_mask]
    y1 = np.full_like(x1, y_new[jc])
    u1 = ut[x_mask, jc]
    v1 = vt[x_mask, jc]
    
    # vertical transect (constant x = x[ic])
    y_mask = np.abs(y_new - nyc) <= r
    y2 = y_new[y_mask]
    x2 = np.full_like(y2, x_new[ic])
    u2 = ut[ic, y_mask]
    v2 = vt[ic, y_mask]
    
    xc, yc, w, Q, _, _, A0 = dopioe(x1, y1, u1, v1, x2, y2, u2, v2)

    cyc_DOPIOE = 'CE' if w < 0 else 'AE'
    
    if (cyc_DOPIOE != cyc) or (np.hypot(nxc - xc, nyc - yc) > 50000):
        return np.nan, np.nan, np.nan, np.array([[np.nan, np.nan],
                                                [np.nan, np.nan]]), np.nan, np.nan, np.nan, np.nan, np.nan
    else:
        w *= 1e-3 # to s^-1

        radii = find_directional_radii(ut, vt, X_new, Y_new, xc, yc, Q)
        R = np.mean([radii['up'], radii['right'], radii['down'], radii['left']])

        q11, q12, q22 = Q[0,0], Q[0,1], Q[1,1]
        dx, dy = X_new - xc, Y_new - yc
        rho2 = q11*dx**2 + 2*q12*dx*dy + q22*dy**2
        rho_search = np.sqrt(np.where(rho2 < 0, np.nan, rho2))
        
        mask_outer = rho_search < max(min(R*1.75, 200000), 30000) 
        axi, ayi, aui, avi = X_new[mask_outer], Y_new[mask_outer], ut[mask_outer], vt[mask_outer]
        
        df = psi_params(xc, yc, Q, axi, ayi, aui, avi)
        Rc, psi0, A = fit_psi_params(df.rho2, df.Qr, df.vt, A0=A0, Rc_max=200000)

        if np.sign(A) != np.sign(w):
            Rc, psi0, A = np.nan, np.nan, A0
            return xc, yc, w, Q, Rc, psi0, A, R, df
        
    return xc, yc, w, Q, Rc, psi0, A, R, df








# Finding Rc
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

def fit_psi_params(rho2, Qr, vt, A0=None, Rc0=None, plot=False, ax=None,
                   maxfev=10000, Rc_max=1e5, r2_flag=False,
                   rho_plot_max=None, n_curve=400):
    import numpy as np
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    import pandas as pd

    d = pd.DataFrame({'rho2': rho2, 'Qr': Qr, 'vt': vt})
    m = (np.isfinite(d.rho2) & np.isfinite(d.vt) & np.isfinite(d.Qr) & (d.rho2 >= 0) & (d.Qr != 0))
    if not np.any(m):
        # "No valid rows after masking."
        return (np.nan, np.nan, np.nan, np.nan) if r2_flag else (np.nan, np.nan, np.nan)
        # raise ValueError("No valid rows after masking.")
    rho2 = d.rho2.values[m]
    vt   = d.vt.values[m]
    Qr   = d.Qr.values[m]

    vt = vt * (np.sqrt(rho2) / Qr)

    def vt_model(rho2_, A, Rc):
        return 2.0 * A * np.sqrt(rho2_) * np.exp(-rho2_ / (Rc**2))

    i = np.nanargmax(np.abs(vt))
    rho_max = np.sqrt(rho2[i])
    if Rc0 is None:
        Rc0 = max(rho_max * np.sqrt(2.0), 1e-6)

    denom = 2.0 * np.sqrt(rho2) * np.exp(-rho2 / (Rc0**2))
    ok = np.abs(denom) > 0
    if A0 is None:
        A0 = np.nanmedian(vt[ok] / denom[ok]) if np.any(ok) else 0.0
    if not np.isfinite(A0):
        A0 = 0.0

    popt, _ = curve_fit(vt_model, rho2, vt, p0=[A0, Rc0],
                        bounds=([-np.inf, 1e-8], [np.inf, np.inf]),
                        maxfev=maxfev)
    A_opt, Rc_opt = popt
    if Rc_opt > Rc_max:
        A_opt, Rc_opt = A0, Rc0

    psi0_opt = -A_opt * Rc_opt**2

    vt_fit = vt_model(rho2, *popt)
    ss_res = np.sum((vt - vt_fit)**2)
    ss_tot = np.sum((vt - np.mean(vt))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot != 0 else np.nan

    if plot:
        r_data = np.sqrt(rho2)
        if rho_plot_max is None:
            rho_plot_max = float(np.nanmax(r_data)) if r_data.size else Rc_opt
        r_grid = np.linspace(0.0, rho_plot_max, n_curve)
        vt_grid = vt_model(r_grid**2, A_opt, Rc_opt)

        if ax is None:
            _, ax = plt.subplots()
        # ax.scatter(r_data, np.abs(vt), s=20, label='Observed', marker='x')
        ax.scatter(r_data, np.abs(vt), s=20, label='Observed', marker='.')
        # ax.plot(r_grid, np.abs(vt_grid), label='Fit', lw=2, color='#ff7f0e')
        ax.plot(r_grid, np.abs(vt_grid), label='Fit', lw=3, color='#ff7f0e')
        # ax.axvline(x=Rc_opt/np.sqrt(2), ls='--', label=r'$\rho_{\max}$', lw=2, color='#ff7f0e')
        ax.axvline(x=Rc_opt/np.sqrt(2), ls='--', label=r'$\rho_{\max}$', lw=3, color='#ff7f0e')
        ax.set_xlabel(r'$\rho$')
        ax.set_ylabel(r'$|v_t^\star|$')
        # ax.legend()
        ax.set_title(f'Best Fit: A={A_opt:.4g}, Rc={Rc_opt:.4g}, psi0={psi0_opt:.4g}, R²={r2:.2f}')

    return (Rc_opt, psi0_opt, A_opt, r2) if r2_flag else (Rc_opt, psi0_opt, A_opt)

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

def psi_params(xc, yc, Q, xi, yi, ui, vi):
    dx, dy = xi - xc, yi - yc
    rho2 = Q[0,0]*dx**2 + 2*Q[1,0]*dx*dy + Q[1,1]*dy**2
    Qr = np.sqrt((Q[0,0]*dx + Q[1,0]*dy)**2 + (Q[1,0]*dx + Q[1,1]*dy)**2)
    vt = tangential_velocity(xi, yi, ui, vi, xc, yc, Q)
    df = pd.DataFrame({'rho2': rho2, 'Qr': Qr, 'vt': vt})
    return df

def plot_max_vt(X, Y, xc, yc, q11, q12, q22, Rc):
    dx = X - xc
    dy = Y - yc
    rho2 = q11*dx**2 + 2*q12*dx*dy + q22*dy**2

    mask = np.isclose(rho2, Rc**2/2, rtol=1e-3, atol=1e-3)
    x_ell = X[mask]
    y_ell = Y[mask]
    return x_ell, y_ell














################################################ TILT ################################################
def bearing(a, b):
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    angle_rad = np.arctan2(dx, dy)  # note the order: dx, dy
    angle_deg = np.degrees(angle_rad)
    bearing = (angle_deg + 360) % 360
    return bearing

def compute_tilt_data(dic, eddy, num=6, depth_int=10, max_depth=1000, var='Depth'):
    
    df_tilt_data = pd.DataFrame(columns=['Eddy', 'Day', 'TiltDis', 'TiltDir'])
    
    diffs_x = {}
    diffs_y = {}
        
    for d, day in enumerate(dic.keys()):
    
        df = dic[day].copy()
        if var == 'Depth':
            df[var] = -df[var]
        df = df[df[var] <= max_depth]
        # don’t drop rows — keep all depths, even if x or y are NaN
        df = df.set_index(var).sort_index()

        if len(df):
            depths = df.index.values
            # interpolate at every 10 m from 0 to max_depth
            target_depths = np.arange(0, max_depth+1, depth_int)
            valid = target_depths[
                (target_depths >= depths.min()) &
                (target_depths <= depths.max())
            ]
            if len(valid) < 2:
                continue
        
            x_i = np.interp(valid, depths, df['x'].values, left=np.nan, right=np.nan)
            y_i = np.interp(valid, depths, df['y'].values, left=np.nan, right=np.nan)
        
            dx = np.diff(x_i)
            dy = np.diff(y_i)
        
            # use the actual depth levels (valid[:-1]) as the Series index
            idx = valid[:-1]
            diffs_x[f'$t_{{{d}}}$'] = pd.Series(dx, index=idx)
            diffs_y[f'$t_{{{d}}}$'] = pd.Series(dy, index=idx)

        else:
            idx = [depth_int]
            diffs_x[f'$t_{{{d}}}$'] = pd.Series(np.array([np.nan]*len(idx)), index=idx)
            diffs_y[f'$t_{{{d}}}$'] = pd.Series(np.array([np.nan]*len(idx)), index=idx)
            
    
    # now construct your DataFrames simply by passing the dict-of-series:
    df_X_all = pd.DataFrame(diffs_x)
    df_Y_all = pd.DataFrame(diffs_y)
    
    for ref_day in range(num //2, len(dic) - num //2):
    
        df_X = df_X_all.iloc[:, ref_day - num // 2:ref_day + num // 2 + 1]
        df_Y = df_Y_all.iloc[:, ref_day - num // 2:ref_day + num // 2 + 1]

        if var == 'rho': # rhos are not lined up nicely like depths.
            df_X = df_X.dropna()
            df_Y = df_Y.dropna()
        
        # Calculation of variability at each depth
        df_data = pd.DataFrame()
        df_data[r'$\Delta x$'] = df_X.mean(axis=1)
        df_data[r'$\Delta y$'] = df_Y.mean(axis=1)
        df_data[r'$\sum{\Delta x}$'] = df_data[r'$\Delta x$'].cumsum()
        df_data[r'$\sum{\Delta y}$'] = df_data[r'$\Delta y$'].cumsum()
        df_data[r'$\sigma^2_{\Delta x}$'] = df_X.var(axis=1)
        df_data[r'$\sigma^2_{\Delta y}$'] = df_Y.var(axis=1)
        df_data[r'Total $\sigma^2$'] = df_data[r'$\sigma^2_{\Delta x}$'] + df_data[r'$\sigma^2_{\Delta y}$']
        df_data['weight'] = 1 / df_data[r'Total $\sigma^2$']
        df_data[var] = df_data.index 
        df_data
        
        # Line of Best Fit
        
        # your data arrays of shape (N,)
        x = df_data[r'$\sum{\Delta x}$'].values
        y = df_data[r'$\sum{\Delta y}$'].values
        z = df_data[var].values
        w = df_data['weight'].values
        
        # 1. compute weighted mean
        W = np.sum(w)
        mean = np.array([np.dot(w, x),
                         np.dot(w, y),
                         np.dot(w, z)]) / W
        
        # 2. center and weight the data
        X = np.vstack((x, y, z)).T
        Xc = X - mean
        Xw = Xc * np.sqrt(w)[:, None]
        
        # 3. SVD on weighted, centered data
        try:
            flag = 0
            _, _, Vt = np.linalg.svd(Xw, full_matrices=False)
        except Exception:
            flag = 1
            
        if flag:
            
            df_tilt_data.loc[len(df_tilt_data)] = {'Eddy': eddy, 'Day': list(dic.keys())[ref_day][3:], 'TiltDis': np.nan, 'TiltDir': np.nan}
            
        else:
            
            direction = Vt[0]   # principal axis
            
            # The best-fit line is:  p(t) = mean + t * direction
            t = np.linspace((np.max(z)-mean[2])/direction[2], (np.min(z)-mean[2])/direction[2], 2)    
            p = mean[None, :] + t[:, None] * direction  
            # or equivalently
            p = mean + np.outer(t, direction)          
            
            # then split back out if you need x,y,z separately:
            x_line, y_line, z_line = p.T
            
            tilt_dist = np.hypot(x_line[0]-x_line[1], y_line[0]-y_line[1])
            
            top_idx = np.where(np.abs(z_line)==np.min(np.abs(z_line)))[0][0]
            if top_idx == 1:
                btm_idx = 0
            else:
                btm_idx = 1
            top = [x_line[top_idx], y_line[top_idx], z_line[top_idx]]
            btm = [x_line[btm_idx], y_line[btm_idx], z_line[btm_idx]]
            tilt_direc = ( bearing(btm, top) + 20 ) % 360
        
            df_tilt_data.loc[len(df_tilt_data)] = {'Eddy': eddy, 'Day': list(dic.keys())[ref_day][3:], 'TiltDis': tilt_dist, 'TiltDir': tilt_direc}
        
    df_tilt_data['Day'] = df_tilt_data['Day'].astype(int)
        
    return df_tilt_data






    





    



















def unit_det(Q, symmetrize=True):
    Q = np.asarray(Q, dtype=float)
    d = np.linalg.det(Q)
    if not np.isfinite(d) or d <= 0: raise ValueError("det(Q) must be positive and finite")
    s = 1/np.sqrt(d)
    return s*Q, s

def gaussian_vel_reconstruction(xc, yc, q11, q12, q22, Rc, psi0, X=None, Y=None):

    if X is None:
        width = 200
        x = np.linspace(xc-width, xc+width, 51)
        y = np.linspace(yc-width, yc+width, 51)
        X, Y = np.meshgrid(x, y)

    dx, dy = X-xc, Y-yc
    rho      = q11*dx**2 + 2*q12*dx*dy + q22*dy**2
    rho_x    = 2*q11*dx   + 2*q12*dy
    rho_y    = 2*q12*dx   + 2*q22*dy
    exp_t    = np.exp(-rho/Rc**2)
    u   =  psi0/Rc**2 * rho_y * exp_t
    v   = -psi0/Rc**2 * rho_x * exp_t

    return u, v, X, Y

from scipy.sparse import diags, eye, kron, csr_matrix
from scipy.sparse.linalg import spsolve
def solve_w(U, V, x, y, z, f=-7.7e-5, N2=5e-3):

    # tic = time.time()

    x = x * 1000
    y = y * 1000
    # z = z * 1000

    dx, dy, dz = x[1] - x[0], y[1] - y[0], z[1] - z[0]
    nx, ny, nz = U.shape

    # build Q
    dudz = np.gradient(U, dz, axis=2)
    dvdz = np.gradient(V, dz, axis=2)
    dvdx = np.gradient(V, dx, axis=0)
    dvdy = np.gradient(V, dy, axis=1)
    dudx = np.gradient(U, dx, axis=0)
    dudy = np.gradient(U, dy, axis=1)

    Qx =  f * (dudz * dvdx + dvdz * dvdy)
    Qy = -f * (dudz * dudx + dvdz * dudy)

    S = 2 * (np.gradient(Qx, dx, axis=0) + np.gradient(Qy, dy, axis=1))
    b = S.ravel(order='F')

    # finite difference operators
    ex, ey = np.ones(nx), np.ones(ny)
    Lx = diags([ex, -2*ex, ex], [-1, 0, 1], shape=(nx, nx)) / dx**2
    Ly = diags([ey, -2*ey, ey], [-1, 0, 1], shape=(ny, ny)) / dy**2
    Ix, Iy = eye(nx), eye(ny)

    # nonuniform Lz
    Lz = np.zeros((nz, nz))
    for i in range(1, nz - 1):
        dzm = z[i] - z[i - 1]
        dzp = z[i + 1] - z[i]
        Lz[i, i - 1] =  2 / (dzm * (dzm + dzp))
        Lz[i, i]     = -2 / (dzm * dzp)
        Lz[i, i + 1] =  2 / (dzp * (dzm + dzp))
    Lz = csr_matrix(Lz)
    Iz = eye(nz)

    A = N2 * (kron(kron(Iz, Iy), Lx) + kron(kron(Iz, Ly), Ix)) \
      + f**2 * kron(kron(Lz, Iy), Ix)

    w = spsolve(A.tocsr(), b).reshape((nx, ny, nz), order='F')

    # toc = time.time()
    # print(f"Elapsed time: {toc - tic:.4f} seconds")

    return w

def extract_transect_center(u, v, X, Y, x0, y0, r=30):

    x, y = X[:,0], Y[0,:]
    
    # find grid point closest to the eddy centre
    dis = np.hypot(X - x0, Y - y0)
    ic, jc = np.unravel_index(np.argmin(dis), dis.shape)

    # horizontal transect (constant y = y[jc])
    x_mask = np.abs(x - x0) < r
    x1 = x[x_mask]
    y1 = np.full_like(x1, y[jc])
    u1 = u[x_mask, jc]
    v1 = v[x_mask, jc]

    # vertical transect (constant x = x[ic])
    y_mask = np.abs(y - y0) < r
    y2 = y[y_mask]
    x2 = np.full_like(y2, x[ic])
    u2 = u[ic, y_mask]
    v2 = v[ic, y_mask]

    # find intersection point
    points1 = set(zip(x1, y1))
    points2 = set(zip(x2, y2))
    common = points1 & points2
    if not common:
        raise ValueError("No common points found.")
    center_x, center_y = common.pop()

    # compute offsets
    xx = x1 - center_x
    yy = y2 - center_y

    return {
        'x1': x1, 'y1': y1, 'u1': u1, 'v1': v1,
        'x2': x2, 'y2': y2, 'u2': u2, 'v2': v2,
        'center_x': center_x, 'center_y': center_y,
        'xx': xx, 'yy': yy
    }

def calc_ow(u, v, dx, dy):
    dudy, dudx = np.gradient(u, dy, dx)
    dvdy, dvdx = np.gradient(v, dy, dx)
    
    # normal strain, shear strain, vorticity
    S_n   = dudx - dvdy
    S_s   = dvdx + dudy
    omega = dvdx - dudy
    
    # Okubo–Weiss parameter
    OW = S_n**2 + S_s**2 - omega**2
    return OW

def plot_isosurface(ax, Xn, Yn, zn, Un, Vn, level=-0.2, elev=13, azim=135, flag=False):
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from skimage.measure import marching_cubes
    from scipy.ndimage import map_coordinates
    dx, dy = Xn[0,1]-Xn[0,0], Yn[1,0]-Yn[0,0]
    ow = calc_ow(Un, Vn, dx*1000, dx*1000, flag=flag)
    sigma_ow = normalize_matrix(ow)

    Xg, Yg, Zg = np.meshgrid(Xn[:, 0], Yn[0, :], zn/1000, indexing='xy' if flag else 'ij')

    verts, faces, normals, values = marching_cubes(sigma_ow, level=level)
    pts = verts.T
    real_x = map_coordinates(Xg, pts, order=1)
    real_y = map_coordinates(Yg, pts, order=1)
    real_z = map_coordinates(Zg, pts, order=1)
    real_verts = np.vstack((real_x, real_y, real_z)).T

    mesh = Poly3DCollection(real_verts[faces], alpha=0.3)
    mesh.set_edgecolor('k')
    ax.add_collection3d(mesh)
    ax.set_xlim(Xg.min(), Xg.max())
    ax.set_ylim(Yg.min(), Yg.max())
    ax.set_zlim(Zg.min(), Zg.max())
    ax.invert_zaxis()
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_zlabel('Depth (km)')
    ax.view_init(elev=elev, azim=azim)

def smooth(x, y, num=1000, window=100):

    from scipy.interpolate import interp1d
    """
    Smooth x vs y by:
      1) interpolating onto a uniform y-grid
      2) applying a running nan‑mean of width `window`
      3) interpolating back to the original y

    Any output point whose window contains only NaNs will be NaN.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)

    # 1) uniform y grid & interpolation
    y_uniform = np.linspace(y.min(), y.max(), num)
    f_interp = interp1d(y, x, kind='linear', fill_value='extrapolate')
    x_uniform = f_interp(y_uniform)

    # 2) running nan-mean on the uniform grid
    def moving_nanmean(a, win):
        """Return an array the same size as a, where each point is the mean of
        the surrounding win points (including itself), ignoring NaNs."""
        valid = ~np.isnan(a)
        # fill NaNs with zero so they don't contribute to sum
        a_fill = np.where(valid, a, 0.0)
        # convolution kernels
        kernel = np.ones(win, dtype=float)
        # sums of valid data
        sums   = np.convolve(a_fill, kernel, mode='same')
        # counts of valid entries
        counts = np.convolve(valid.astype(float), kernel, mode='same')
        # compute mean; outside domain or where counts==0 gives NaN
        with np.errstate(divide='ignore', invalid='ignore'):
            mu = sums / counts
        mu[counts == 0] = np.nan
        return mu

    x_smooth_uniform = moving_nanmean(x_uniform, window)

    # 3) back to original y
    f_smooth = interp1d(y_uniform, x_smooth_uniform,
                        kind='linear', fill_value='extrapolate')
    x_smooth = f_smooth(y)

    return x_smooth

def robust_smooth(y, x, win=5, poly=3, k=3.5, s=None):
    from scipy.signal import savgol_filter
    from scipy.interpolate import UnivariateSpline
    i = np.argsort(y); y0, x0 = y[i], x[i]
    x1 = savgol_filter(x0, win|1, poly, mode='interp')
    r = x0 - x1
    mad = np.median(np.abs(r - np.median(r))) + 1e-12
    mask = np.abs(r) <= k*1.4826*mad
    spl = UnivariateSpline(y0[mask], x0[mask], s=s)
    return spl(y0), y0, x0

def normalize_matrix(matrix, mask_value=np.nan):
    valid_mask = np.where(matrix == mask_value, 0, 1)
    valid_mean = np.nansum(matrix) / np.sum(valid_mask)
    valid_std = np.sqrt(np.nansum(valid_mask * (matrix - valid_mean) ** 2) / np.sum(valid_mask))
    return (matrix - valid_mean) / valid_std

def rossby_number(vort, lat_deg):
    omega = 7.2921e-5  # Earth's rotation rate (rad/s)
    phi = np.radians(lat_deg)
    f = np.abs(2 * omega * np.sin(phi))
    return vort / f

def deg_to_m(lat):
    R = 6371000
    rad = np.radians(lat)
    return (np.pi/180)*R*np.sqrt((np.cos(rad))**2 + 1)

def find_icjc(x, y, X_grid, Y_grid):
    from scipy.spatial import cKDTree
    tree = cKDTree(np.column_stack((X_grid.ravel(), Y_grid.ravel())))
    xcs = np.asarray(x, dtype=float)
    ycs = np.asarray(y, dtype=float)
    valid = np.isfinite(xcs) & np.isfinite(ycs)
    ics = np.full(xcs.shape, -1, dtype=int)
    jcs = np.full(ycs.shape, -1, dtype=int)
    if valid.any():
        _, ind = tree.query(np.column_stack((xcs[valid], ycs[valid])))
        ii, jj = np.unravel_index(ind, X_grid.shape)
        ics[valid] = ii
        jcs[valid] = jj
    return ics, jcs

def nencioli(u, v, lon, lat, a, b):
    """
    Identify the points in the domain which satisfy the four velocity constraints for eddy detection.

    Parameters:
    - u, v: 2D velocity fields for u and v components
    - lon, lat: Longitude and Latitude matrices
    - mask: Matrix defining sea (1) and land points (0)
    - a, b: Parameters used for constraints

    Returns:
    - eddy_uv: Positions that satisfy the first two constraints (for debugging)
    - eddy_c: Positions satisfying the first three constraints (for debugging)
    - eddy: Positions of the eddy centers with their type (cyclonic=1, anticyclonic=-1) (opposite for me)
    """

    borders = max(a, b) + 1

    # Compute velocity magnitude
    vel = np.sqrt(u**2 + v**2)

    # Initialize arrays for storing eddy centers
    eddy_uv = np.zeros((0, 2))
    eddy_c = np.zeros((0, 2))
    eddy = np.zeros((0, 3))

    # Get domain dimensions
    bound = vel.shape

    # Loop through each latitudinal section
    for i in range(borders, len(v) - borders + 1):
        wrk = v[i, :]  # Latitudinal section of v

        # First constraint: zero crossing in v component
        s = np.sign(wrk)
        indx = np.where(np.diff(s) != 0)[0]
        indx = indx[(indx >= borders) & (indx < len(wrk) - borders)]

        for ii in indx:
            var = 0  # Eddy type (0 = no eddy, 1 = cyclonic, -1 = anticyclonic)
            if wrk[ii] >= 0:  # Anticyclonic
                if wrk[ii - a] > wrk[ii] and wrk[ii + 1 + a] < wrk[ii + 1]:
                    var = -1
            elif wrk[ii] < 0:  # Cyclonic
                if wrk[ii - a] < wrk[ii] and wrk[ii + 1 + a] > wrk[ii + 1]:
                    var = 1

            # Second constraint: u component reversal
            if var != 0:
                if var == -1:
                    if (u[i - a, ii] <= 0 and u[i - a, ii] <= u[i - 1, ii] and
                        u[i + a, ii] >= 0 and u[i + a, ii] >= u[i + 1, ii]):
                        eddy_uv = np.vstack([eddy_uv, [lat[i, ii], lon[i, ii]], [lat[i, ii + 1], lon[i, ii + 1]]])
                    else:
                        var = 0
                elif var == 1:
                    if (u[i - a, ii] >= 0 and u[i - a, ii] >= u[i - 1, ii] and
                        u[i + a, ii] <= 0 and u[i + a, ii] <= u[i + 1, ii]):
                        eddy_uv = np.vstack([eddy_uv, [lat[i, ii], lon[i, ii]], [lat[i, ii + 1], lon[i, ii + 1]]])
                    else:
                        var = 0

                # Third constraint: velocity minimum
                if var != 0:
                    srch = vel[i - b:i + b, ii - b:ii + b + 1]
                    slat = lat[i - b:i + b, ii - b:ii + b + 1]
                    slon = lon[i - b:i + b, ii - b:ii + b + 1]
                    X, Y = np.unravel_index(np.argmin(srch), srch.shape)
                    srch2 = vel[max(i - b + X - 1 - b, 0):min(i - b + X - 1 + b, bound[0]),
                                max(ii - b + Y - 1 - b, 0):min(ii - b + Y - 1 + b, bound[1])]

                    if np.min(srch2) == np.min(srch):
                        eddy_c = np.vstack([eddy_c, [slat[X, Y], slon[X, Y]]])
                    else:
                        var = 0

                # Fourth constraint: vector rotation (simplified version)
                d = a - 1
                if var != 0:
                    # Find indices of the estimated center in the large domain
                    i1, i2 = np.where((lat == slat[X, Y]) & (lon == slon[X, Y]))

                    i1, i2 = int(i1[0]), int(i2[0])
                    
                    # Extract velocities within "a-1" points from the estimated center
                    u_small = u[max(i1 - d, 0):min(i1 + d, bound[0]), max(i2 - d, 0):min(i2 + d, bound[1])]
                    v_small = v[max(i1 - d, 0):min(i1 + d, bound[0]), max(i2 - d, 0):min(i2 + d, bound[1])]
                    
                    # Apply constraint only if there are no NaNs in u_small
                    if not np.isnan(u_small).any():
                        # Boundary velocities
                        u_bound = np.concatenate([u_small[0, :], u_small[1:, -1], u_small[-1, -2::-1], u_small[-2::-1, 0]])
                        v_bound = np.concatenate([v_small[0, :], v_small[1:, -1], v_small[-1, -2::-1], v_small[-2::-1, 0]])

                        # Vector defining which quadrant each boundary vector belongs to
                        quadrants = np.zeros_like(u_bound)
                        quadrants[(u_bound >= 0) & (v_bound >= 0)] = 1
                        quadrants[(u_bound < 0) & (v_bound >= 0)] = 2
                        quadrants[(u_bound < 0) & (v_bound < 0)] = 3
                        quadrants[(u_bound >= 0) & (v_bound < 0)] = 4
                        
                        # Identify the first fourth quadrant vector
                        spin = np.where(quadrants == 4)[0]
                        
                        # Apply the constraint only if the rotation is complete and not all vectors are in the fourth quadrant
                        if spin.size > 0 and spin.size != quadrants.size:
                            # If vectors start in the 4th quadrant, add 4 to all quadrant positions after the first occurrence
                            if spin[0] == 0:
                                spin = np.where(quadrants != 4)[0]
                                spin = spin[0] - 1
                                
                            if not isinstance(spin, np.ndarray):
                                spin = np.array([int(spin)])
                            quadrants[spin[-1] + 1:] += 4
                            
                            # Inspect vector rotation: no consecutive vectors should be more than one quadrant apart
                            # and there should be no backward rotation
                            if not np.any(np.diff(quadrants) > 1) and not np.any(np.diff(quadrants) < 0):
                                eddy = np.vstack([eddy, [slat[X, Y], slon[X, Y], var]])


    # Process eddy results (sorting and removing duplicates)
    eddy = np.unique(eddy, axis=0)
    eddy_uv = np.unique(eddy_uv, axis=0)
    eddy_c = np.unique(eddy_c, axis=0)
    # Adjust for the Southern Hemisphere (flip cyclonic/anticyclonic labels)
    # eddy[eddy[:, 0] < 0, 2] = -eddy[eddy[:, 0] < 0, 2]
    eddy[:, 2] = -eddy[:, 2]
    # Swap for personal preference 
    eddy[:, [0, 1]] = eddy[:, [1, 0]]

    return eddy_uv, eddy_c, eddy

def tilt_distance_LI(x, y, z, zmin=None, zmax=None):
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

def project_sadcp_to_transect(x, y, u, v):

    num_points = len(x)

    df = pd.DataFrame({'x': x, 'y': y, 'u': u, 'v': v})
    df = df.sort_values(by='x').reset_index(drop=True)
    x, y, u, v = np.array(df['x']), np.array(df['y']), np.array(df['u']), np.array(df['v'])

    # Fit a best-fit line (y = m*x + c) using linear regression.
    A = np.vstack([x, np.ones(len(x))]).T
    m, _ = np.linalg.lstsq(A, y, rcond=None)[0]
    # Create a unit direction vector along the best-fit line.
    direction = np.array([1, m])
    direction = direction / np.linalg.norm(direction)
    # Project each (x, y) point onto the line using the first point as reference.
    p0 = np.array([x[0], y[0]])
    points = np.column_stack((x, y))
    projections = np.dot(points - p0, direction)
    # Sort projections and corresponding velocities.
    sort_idx = np.argsort(projections)
    proj_sorted = projections[sort_idx]
    u_sorted = u[sort_idx]
    v_sorted = v[sort_idx]
    # Interpolate the velocity components on evenly spaced positions.
    new_distances = np.linspace(proj_sorted.min(), proj_sorted.max(), num_points)
    new_u = np.interp(new_distances, proj_sorted, u_sorted)
    new_v = np.interp(new_distances, proj_sorted, v_sorted)
    # Calculate new (x, y) points along the line.
    new_points = p0 + np.outer(new_distances, direction)
    # Create and return the output dataframe.
    x, y = new_points[:, 0], new_points[:, 1]

    cos_theta = 1 / np.sqrt(1+m**2)
    sin_theta = m / np.sqrt(1+m**2)
    V_N = -new_u * sin_theta + new_v * cos_theta
    V_T = new_v * sin_theta + new_u * cos_theta
        
    df_projected = pd.DataFrame({
        'x': x,
        'y': y,
        'u': new_u,
        'v': new_v,
        'V_N': V_N,
        'V_T': V_T,
        # 'l': new_distances
    })

    df_projected = df_projected.sort_values(by='x').reset_index(drop=True)

    df_projected['l'] = np.hypot(df_projected['x']-df_projected['x'].iloc[0], df_projected['y']-df_projected['y'].iloc[0])

    return df_projected, m

def translate_moca_results(x_l_start, y_l_start, m, l0, r0):

    x0 = (l0-r0*m)/np.sqrt(1+m**2) + x_l_start
    y0 = (l0*m+r0)/np.sqrt(1+m**2) + y_l_start
    
    return x0, y0

def eccentricity_from_Q(Q):
    eigvals = np.linalg.eigvals(Q)
    lam1, lam2 = np.sort(eigvals)
    axis_ratio = lam1 / lam2
    eccentricity = np.sqrt(1 - (lam2/lam1)**2)
    return axis_ratio, eccentricity

def eccentricity(df):
    lam1 = 0.5*(df.sq11 + df.sq22) + np.sqrt(((df.sq11 - df.sq22)/2)**2 + df.sq12**2)  # major
    lam2 = 0.5*(df.sq11 + df.sq22) - np.sqrt(((df.sq11 - df.sq22)/2)**2 + df.sq12**2)  # minor
    return np.sqrt(1 - (lam2/lam1)**2)

def axis_ratio(df):
    lam1 = 0.5*(df.sq11 + df.sq22) + np.sqrt(((df.sq11 - df.sq22)/2)**2 + df.sq12**2)  # major
    lam2 = 0.5*(df.sq11 + df.sq22) - np.sqrt(((df.sq11 - df.sq22)/2)**2 + df.sq12**2)  # minor
    return np.sqrt(lam1/lam2)

def ellipse_aspect_ratio(q11, q12, q22, eps=1e-12):
    # q11 = np.asarray(np.abs(q11), float)
    q11 = np.asarray(q11, float)
    q12 = np.asarray(q12, float)
    # q22 = np.asarray(np.abs(q22), float)
    q22 = np.asarray(q22, float)
    tr  = q11 + q22
    rad = np.sqrt((q11 - q22)**2 + 4*q12**2)
    lam_min = 0.5*(tr - rad)
    lam_max = 0.5*(tr + rad)
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.sqrt(lam_max / np.maximum(lam_min, eps))
    # mark non-SPD cases explicitly as nan
    r = np.where(lam_min > 0, r, np.nan)
    return r.item()

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
            rho2 = (
                row.q11 * dx**2
                + 2 * row.q12 * dx * dy
                + row.q22 * dy**2
            )
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



























