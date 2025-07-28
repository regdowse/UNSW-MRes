import numpy as np

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
    dsdx = np.gradient(sigma, dx, axis=0)
    dsdy = np.gradient(sigma, dy, axis=1)

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
    
    x_c = 0
    y_c = 0

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













    


def moca(l, VT, VN, Rc_init=4.0, psi0_init=200.0):

    if np.any(np.isnan(VT)):
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

    Q11, Q12, Q22 = w/4, 0, w/4
    Q = np.array([[Q11, Q12], [Q12, Q22]])
    xi, yi, ui, vi = l, [0]*len(l), VT, VN
    xc, yc = l0, r0

    # psi0 = estimate_psi0(xi, yi, ui, vi, xc, yc, Q[0,0], Q[1,0], Q[1,1], psi0_init=psi0_init)
    # q11_init = -Rc_init**2/psi0 * Q11
    # q12_init = -Rc_init**2/psi0 * Q12
    # q22_init = -Rc_init**2/psi0 * Q22
    # Rc, q = estimate_Rc(xi, yi, ui, vi, xc, yc, psi0, 
    #                  q11_init=q11_init, q12_init=q12_init, q22_init=q22_init)

    s_init = - Rc_init**2 / psi0_init
    params = fit_eddy(xi, yi, ui, vi, xc, yc, Q[0,0], Q[1,0], Q[1,1], p0=(s_init, Rc_init, psi0_init))
    Rc, psi0, q11, q12, q22 = params['Rc'], params['psi0'], params['q11'], params['q12'], params['q22']
    q = np.array([[q11, q12], [q12, q22]])
    
    return l0, r0, w, Q, Rc, psi0, q

def dopioe(x1, y1, u1, v1, x2, y2, u2, v2, Rc_init=4.0, psi0_init=200.0):

    if np.any(np.isnan(u1)) or np.any(np.isnan(u2)):
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
    
    Q11 = B1 / 2
    Q22 = -C1 / 2
    Q12 = -gamma / 2
    denom = C1 * B1 + gamma**2
    if denom == 0:
        raise ZeroDivisionError("Denominator is zero.")
        
    xc = - (alpha * gamma + beta * C1) / denom  + center_x
    yc = (beta * gamma - alpha * B1) / denom + center_y
    w = 2 * (Q11 + Q22)

    Q = np.array([[Q11, Q12], [Q12, Q22]])

    # Remove duplicates from (x1, y1, u1, v1)
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

    # psi0 = estimate_psi0(xi, yi, ui, vi, xc, yc, Q[0,0], Q[1,0], Q[1,1], psi0_init=psi0_init)
    # q11_init = -Rc_init**2/psi0 * Q11
    # q12_init = -Rc_init**2/psi0 * Q12
    # q22_init = -Rc_init**2/psi0 * Q22
    # Rc, q = estimate_Rc(xi, yi, ui, vi, xc, yc, psi0, 
    #                  q11_init=q11_init, q12_init=q12_init, q22_init=q22_init)

    s_init = - Rc_init**2 / psi0_init
    params = fit_eddy(xi, yi, ui, vi, xc, yc, Q[0,0], Q[1,0], Q[1,1], p0=(s_init, Rc_init, psi0_init))
    Rc, psi0, q11, q12, q22 = params['Rc'], params['psi0'], params['q11'], params['q12'], params['q22']
    q = np.array([[q11, q12], [q12, q22]])

    return xc, yc, w, Q, Rc, psi0, q

def espra(xi, yi, ui, vi, Rc_init=4.0, psi0_init=200.0):

    if np.any(np.isnan(ui)):
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
    
    from scipy.optimize import least_squares

    def residuals(params, x, y, u_i, v_i):
        xc, yc, Q11, Q12, Q22 = params
        u = -2 * Q22 * (y - yc) - 2 * Q12 * (x - xc)
        v =  2 * Q11 * (x - xc) + 2 * Q12 * (y - yc)
        return np.concatenate([(u - u_i), (v - v_i)])

    def fit_params(x, y, u_i, v_i):
        xc_init, yc_init = np.mean(x), np.mean(y)
        Q11_init, Q12_init, Q22_init = 1.0, 0.0, 1.0  # Initial guesses
        params_init = [xc_init, yc_init, Q11_init, Q12_init, Q22_init]
        result = least_squares(residuals, params_init, args=(x, y, u_i, v_i))
        return result.x 

    xc, yc, Q11, Q12, Q22 = fit_params(xi, yi, ui, vi)

    w = 2*(Q11 + Q22)

    Q = np.array([[Q11, Q12], [Q12, Q22]])

    # psi0 = estimate_psi0(xi, yi, ui, vi, xc, yc, Q[0,0], Q[1,0], Q[1,1], psi0_init=psi0_init)
    # q11_init = -Rc_init**2/psi0 * Q11
    # q12_init = -Rc_init**2/psi0 * Q12
    # q22_init = -Rc_init**2/psi0 * Q22
    # Rc, q = estimate_Rc(xi, yi, ui, vi, xc, yc, psi0, 
    #                  q11_init=q11_init, q12_init=q12_init, q22_init=q22_init)

    s_init = - Rc_init**2 / psi0_init
    params = fit_eddy(xi, yi, ui, vi, xc, yc, Q[0,0], Q[1,0], Q[1,1], p0=(s_init, Rc_init, psi0_init))
    Rc, psi0, q11, q12, q22 = params['Rc'], params['psi0'], params['q11'], params['q12'], params['q22']
    q = np.array([[q11, q12], [q12, q22]])

    return xc, yc, w, Q, Rc, psi0, q

def fit_eddy(xi, yi, ui, vi, xc, yc, Q11, Q12, Q22, p0=(1,1,1)):
    from scipy.optimize import least_squares
    def residuals(p):
        s, Rc, psi0 = p
        q11, q12, q22 = s*Q11, s*Q12, s*Q22
        dx, dy = xi-xc, yi-yc
        rho      = q11*dx**2 + 2*q12*dx*dy + q22*dy**2
        rho_x    = 2*q11*dx   + 2*q12*dy
        rho_y    = 2*q12*dx   + 2*q22*dy
        exp_t    = np.exp(-rho/Rc**2)
        u_pred   =  psi0/Rc**2 * rho_y * exp_t
        v_pred   = -psi0/Rc**2 * rho_x * exp_t
        return np.concatenate((u_pred-ui, v_pred-vi))
    res = least_squares(residuals, p0)
    s_opt, Rc_opt, psi0_opt = res.x
    return {
      'q11': s_opt*Q11,
      'q12': s_opt*Q12,
      'q22': s_opt*Q22,
      'Rc':  Rc_opt,
      'psi0': psi0_opt,
      's': s_opt
    }


# def estimate_psi0(x, y, u, v, xc, yc, Q11, Q12, Q22, psi0_init=200):
#     from scipy.optimize import minimize_scalar
#     """Least‑squares estimate of ψ0 given core Q_ij and observed (u,v)."""
#     x = np.array(x) if isinstance(x, list) else x
#     y = np.array(y) if isinstance(y, list) else y
#     u = np.array(u) if isinstance(u, list) else u
#     v = np.array(v) if isinstance(v, list) else v

#     sign_guess = -np.sign(Q11)        # + AE, - CE
#     guess     = sign_guess * psi0_init      # magnitude ≈ 100 in your units
#     bounds    = (1e-3, 1e6) if sign_guess > 0 else (-1e6, -1e-3)

#     # shifts and handy algebra
#     dx, dy   = x - xc, y - yc
#     Phi      = Q11*dx**2 + 2*Q12*dx*dy + Q22*dy**2
#     dPhi_dx  = 2*Q11*dx + 2*Q12*dy
#     dPhi_dy  = 2*Q12*dx + 2*Q22*dy

#     mask     = ~np.isnan(u) & ~np.isnan(v)

#     def residual(psi0):
#         e  = np.exp(Phi[mask] / psi0)
#         um = e * dPhi_dy[mask]
#         vm = -e * dPhi_dx[mask]        
#         return np.mean((u[mask]-um)**2 + (v[mask]-vm)**2)

#     res = minimize_scalar(residual, method='bounded',
#                           bounds=bounds, options={'xatol':1e-8})

#     psi0 = res.x
#     if np.abs(psi0) > .9 * np.max(np.abs(bounds)):
#         psi0 = np.nan

#     return psi0 

# def estimate_Rc(xi, yi, ui, vi, xc, yc, psi0, q11_init=1.0, q12_init=0.0, q22_init=1.0, Rc_init=2):

#     from scipy.optimize import least_squares
#     if np.any(np.isnan(ui)):
#         return np.nan

#     def residuals(params, x, y, u_i, v_i, xc, yc, psi0):
#         q11, q12, q22, Rc = params
#         dx, dy = x - xc, y - yc
#         rho   = q11*dx**2 + 2*q12*dx*dy + q22*dy**2
#         rho_x = 2*q11*dx   + 2*q12*dy
#         rho_y = 2*q12*dx   + 2*q22*dy
#         exp_t = np.exp(-rho/Rc**2)
#         u_pred =  psi0/Rc**2 * rho_y * exp_t
#         v_pred = -psi0/Rc**2 * rho_x * exp_t
#         return np.concatenate([u_pred - u_i, v_pred - v_i])

#     # initial guesses + non‐negative bounds on q11,q22,Rc
#     x0 = [q11_init, q12_init, q22_init, Rc_init]
#     lb = [-np.inf,  -np.inf,   -np.inf,     1e-6]
#     ub = [ np.inf,   np.inf,    np.inf,   np.inf]

#     try:
#         res = least_squares(
#             residuals,
#             x0,
#             bounds=(lb, ub),
#             args=(xi, yi, ui, vi, xc, yc, psi0)
#         )
#     except ValueError as e:
#         # if "`x0` is infeasible" in str(e):
#         return np.nan, np.array([[np.nan, np.nan],
#                                  [np.nan, np.nan]])
#         # else:
#         #     raise  # re-raise other ValueErrors
    
#     q11, q12, q22, Rc = res.x
#     q = np.array([[q11, q12], [q12, q22]])
#     return Rc, q

























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

# def gaussian_vel_reconstruction(x0, y0, Q11, Q12, Q22, Rc, psi0, X=None, Y=None, flag=True): 

#     if flag:
#         s = - Rc**2 / psi0
#     else:
#         s = 1
    
#     q11 = s * Q11
#     q12 = s * Q12
#     q22 = s * Q22

#     if X is None:
#         width = 200
#         x = np.linspace(x0-width, x0+width, 51)
#         y = np.linspace(y0-width, y0+width, 51)
#         X, Y = np.meshgrid(x, y)
    
#     dx, dy = X - x0, Y - y0
    
#     phi   = q11*dx**2 + 2*q12*dx*dy + q22*dy**2
#     phi_x = 2*q11*dx  + 2*q12*dy
#     phi_y = 2*q22*dy  + 2*q12*dx
    
#     # 5) build Gaussian streamfunction with that R
#     exp_term = np.exp(-phi / Rc**2)
#     psi_x = -phi_x / Rc**2 * exp_term
#     psi_y = -phi_y / Rc**2 * exp_term
    
#     u =  psi_y * psi0
#     v = -psi_x * psi0

#     return u, v, X, Y

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

def normalize_matrix(matrix, mask_value=np.nan):
    valid_mask = np.where(matrix == mask_value, 0, 1)
    valid_mean = np.nansum(matrix) / np.sum(valid_mask)
    valid_std = np.sqrt(np.nansum(valid_mask * (matrix - valid_mean) ** 2) / np.sum(valid_mask))
    return (matrix - valid_mean) / valid_std






































