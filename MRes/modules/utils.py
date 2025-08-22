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













    


def moca(l, VT, VN):

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

    psi0 = find_optimal_psi0(xi, yi, ui, vi,
                      xc, yc, Q11, Q12, Q22)

    Rc = calc_tang_vel_max_r(xc, yc, xi, yi, ui, vi)

    s = -Rc**2/psi0
    q = s*Q

    return l0, r0, w, Q, Rc, psi0, q

def dopioe(x1, y1, u1, v1, x2, y2, u2, v2):

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

    psi0 = find_optimal_psi0(xi, yi, ui, vi,
                      xc, yc, Q11, Q12, Q22)

    Rc = calc_tang_vel_max_r(xc, yc, xi, yi, ui, vi)

    s = -Rc**2/psi0
    q = s*Q

    return xc, yc, w, Q, Rc, psi0, q

def espra(xi, yi, ui, vi):

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

    psi0 = find_optimal_psi0(xi, yi, ui, vi,
                      xc, yc, Q11, Q12, Q22)

    Rc = calc_tang_vel_max_r(xc, yc, xi, yi, ui, vi)

    s = -Rc**2/psi0
    q = s*Q

    return xc, yc, w, Q, Rc, psi0, q










# Finding psi0

def find_optimal_psi0(xi, yi, ui, vi,
                      xc, yc, Q11, Q12, Q22,
                      bounds=(1e-6, 1e6),
                      method='bounded',
                      exp_clip=700):
    from scipy.optimize import minimize_scalar
    """
    Estimate the ψ₀ that minimises
      R1(ψ₀) = Σ_i [ (u_i + β_i e^{γ_i/ψ₀})^2 + (v_i - α_i e^{γ_i/ψ₀})^2 ],
    where
      γ_i = Q11·dx² + 2·Q12·dx·dy + Q22·dy²,
      α_i = 2·Q11·dx + 2·Q12·dy,
      β_i = 2·Q22·dy + 2·Q12·dx.
    """
    if np.isnan(Q11):
        return np.nan
    sign = 'positive' if Q11 < 0 else 'negative'
    
    xi = np.asarray(xi, dtype=np.float64)
    yi = np.asarray(yi, dtype=np.float64)
    ui = np.asarray(ui, dtype=np.float64)
    vi = np.asarray(vi, dtype=np.float64)

    dx = xi - xc
    dy = yi - yc
    gamma = Q11*dx**2 + 2*Q12*dx*dy + Q22*dy**2
    alpha = 2*Q11*dx + 2*Q12*dy
    beta  = 2*Q22*dy + 2*Q12*dx

    def obj(psi0):
        # Stable exponent
        z = np.clip(gamma / psi0, -exp_clip, exp_clip)
        E = np.exp(z)
        u_pred = -beta * E
        v_pred =  alpha * E
        s = np.sum((ui - u_pred)**2 + (vi - v_pred)**2)
        # Penalise any non-finite objective (keeps optimiser sane)
        return s if np.isfinite(s) else 1e300

    def solve_on_interval(lo, hi):
        return minimize_scalar(obj, bounds=(lo, hi), method=method)

    # Choose intervals
    if sign == 'positive':
        candidates = [solve_on_interval(bounds[0], bounds[1])]
    elif sign == 'negative':
        candidates = [solve_on_interval(-bounds[1], -bounds[0])]
    else:  # auto: try both and pick the best
        r_pos = solve_on_interval(bounds[0],  bounds[1])
        r_neg = solve_on_interval(-bounds[1], -bounds[0])
        candidates = [r_pos, r_neg]

    best = min(candidates, key=lambda r: r.fun if np.isfinite(r.fun) else np.inf)
    return best.x

    

# Finding Rc

def calc_tang_vel(xc, yc, xp, yp, up, vp):
    xp, yp = np.asarray(xp), np.asarray(yp)
    up, vp = np.asarray(up), np.asarray(vp)
    dx = xp - xc
    dy = yp - yc
    r = np.hypot(dx, dy)
    v_theta = (-up * dy + vp * dx) / r
    v_theta = np.where(r>0, v_theta, 0.0)
    return v_theta

def calc_tang_vel_max_r(xc, yc, xp, yp, up, vp, cyc=None):
    xp, yp = np.asarray(xp), np.asarray(yp)
    up, vp = np.asarray(up), np.asarray(vp)
    dx = xp - xc
    dy = yp - yc
    r = np.hypot(dx, dy)
    v_theta = (-up * dy + vp * dx) / r

    if cyc is None:
        v_theta = np.abs(v_theta)
    elif cyc == 'AE':
        mask = v_theta >= 0
        if np.sum(mask) == 0:
            return np.nan
        v_theta = np.abs(v_theta[mask])
        r = r[mask]
    elif cyc == 'CE':
        mask = v_theta <= 0
        if np.sum(mask) == 0:
            return np.nan
        v_theta = np.abs(v_theta[mask])
        r = r[mask]
    else:
        return np.nan

    v_theta = np.where(r>0, v_theta, 0.0)

    i_peak = np.nanargmax(v_theta)
    r_peak =  r[i_peak]
    
    return r_peak

def find_directional_radii(u, v, x, y, xc, yc, calc_tang_vel, return_index=False):
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
            vt = abs(calc_tang_vel(xc, yc, x[i, j], y[i, j], u[i, j], v[i, j]))
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

def eddy_core_radius(r, v_theta):
    r, v = np.asarray(r), np.asarray(v_theta)
    m = ~np.isnan(r)&~np.isnan(v)
    r, v = r[m], v[m]
    order = np.argsort(r)
    r, v = r[order], v[order]

    R2_list, slopes = [], []
    for n in range(2, len(r)+1):
        # slope of vt=r·Ω through the origin
        Ω = np.dot(r[:n], v[:n]) / np.dot(r[:n], r[:n])
        v_fit = Ω * r[:n]
        ss_res = np.sum((v[:n] - v_fit)**2)
        ss_tot = np.sum((v[:n] - v[:n].mean())**2)
        R2_list.append(1 - ss_res/ss_tot)
        slopes.append(Ω)

    R2_list = smooth(R2_list, np.arange(len(R2_list)), num=len(R2_list), window=round(len(R2_list)*.1))

    i_best = int(np.argmax(R2_list))
    r_core   = r[i_best+1]     # +1 because R2_list[0] used r[:2], so index→r[1]
    Ω_uniform = slopes[i_best]
    return r_core, Ω_uniform, R2_list, slopes





















    






    





    





















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

def rossby_number(vort, lat_deg):
    omega = 7.2921e-5  # Earth's rotation rate (rad/s)
    phi = np.radians(lat_deg)
    f = np.abs(2 * omega * np.sin(phi))
    return vort / f

def deg_to_m(lat):
    R = 6371000
    rad = np.radians(lat)
    return (np.pi/180)*R*np.sqrt((np.cos(rad))**2 + 1)

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


































