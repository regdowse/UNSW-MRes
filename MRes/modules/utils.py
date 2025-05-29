import numpy as np

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
    
    return l0, r0, w

def moca_grid(u, v, x, y, nic, njc, r):
    
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

    # Get center coordinates
    center_x = x[nic]
    center_y = y[njc]
    
    # Distance along column (x-direction)
    dx = np.abs(x - center_x)
    ix = np.where(dx <= r)[0]
    v1 = v[ix, njc]
    u1 = u[ix, njc]
    if np.any(np.isnan(u1)):
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    xx = x[ix]
    base = xx[0]
    l = xx - base
    yy = np.full_like(xx, center_y)
    
    root = find_root(l, v1)
    c, b = tang_at_root(l, v1, root)  # c: slope, b: intercept
    a = cubic_interpolate(l, u1, root)
    
    x0 = -b / c 
    y0 = a / c 
    w = 2 * c

    q11, q12, q22 = w/4, 0, w/4
    psi0_opt, Rc_opt = fit_Rc_VN1(xx, yy, v1, q11, q12, q22, x0, y0)
    
    return x0 + base, y0+ yy[0], w, l, xx, yy, u1, v1, Rc_opt, psi0_opt

def dopioe(x1, y1, u1, v1, x2, y2, u2, v2):
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

    if np.any(np.isnan(u1)) or np.any(np.isnan(u2)):
        return np.nan, np.nan, np.nan, np.nan

    points1 = set(zip(x1, y1))
    points2 = set(zip(x2, y2))
    common_points = points1 & points2
    
    if len(common_points) != 1:
        print(f"Warning: Expected 1 common point, found {len(common_points)}.")
        return np.nan, np.nan, np.nan, np.nan
    
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
    
    q11 = B1 / 2
    q22 = -C1 / 2
    q12 = -gamma / 2
    denom = C1 * B1 + gamma**2
    if denom == 0:
        raise ZeroDivisionError("Denominator is zero.")
        
    x0 = - (alpha * gamma + beta * C1) / denom
    y0 = (beta * gamma - alpha * B1) / denom
    w = 2 * (q11 + q22)

    Q = np.array([[q11, q12], [q12, q22]])

    xi = np.concatenate([x1, x2])
    yi = np.concatenate([y1, y2])
    ui = np.concatenate([u1, u2])
    vi = np.concatenate([v1, v2])
    
    Rc_opt, psi0_opt = espra_Rc(xi, yi, ui, vi, x0, y0, q11, q12, q22)

    return x0 + center_x, y0 + center_y, w, Q, Rc_opt, psi0_opt

def dopioe_grid(nic, njc, r, u, v, X, Y):

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
    
    # Get center coordinates
    center_x = X.T[nic, njc]
    center_y = Y.T[nic, njc]
    
    # Distance along column (x-direction)
    x_col = X.T[:, njc]
    dx = np.abs(x_col - center_x)
    ix = np.where(dx <= r)[0]
    u1 = u[ix, njc]
    v1 = v[ix, njc]
    x1 = x_col[ix]
    y1 = Y.T[ix, njc]
    xx = x1 - center_x
    
    # Distance along row (y-direction)
    y_row = Y.T[nic, :]
    dy = np.abs(y_row - center_y)
    iy = np.where(dy <= r)[0]
    u2 = u[nic, iy]
    v2 = v[nic, iy]
    x2 = X.T[nic, iy]
    y2 = y_row[iy]
    yy = y2 - center_y

    if np.any(np.isnan(u1)) or np.any(np.isnan(u2)):
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

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
    
    q11 = B1 / 2
    q22 = -C1 / 2
    q12 = -gamma / 2
    denom = C1 * B1 + gamma**2
    if denom == 0:
        raise ZeroDivisionError("Denominator is zero.")
        
    x0 = - (alpha * gamma + beta * C1) / denom + center_x
    y0 = (beta * gamma - alpha * B1) / denom + center_y
    w = 2 * (q11 + q22)

    Q = np.array([[q11, q12], [q12, q22]])

    xi = np.concatenate([x1, x2])
    yi = np.concatenate([y1, y2])
    ui = np.concatenate([u1, u2])
    vi = np.concatenate([v1, v2])
    
    Rc_opt, psi0_opt = espra_Rc(xi, yi, ui, vi, x0, y0, q11, q12, q22)

    return x0, y0, w, Q, Rc_opt, psi0_opt, x1, y1, u1, v1, x2, y2, u2, v2, xx, yy

def espra(xi, yi, ui, vi):

    if np.any(np.isnan(ui)):
        return np.nan, np.nan, np.array([[np.nan, np.nan], [np.nan, np.nan]]), np.nan, np.nan, np.nan
    
    from scipy.optimize import least_squares

    def residuals(params, x, y, u_i, v_i):
        x0, y0, q11, q12, q22 = params
        u = -2 * q22 * (y - y0) - 2 * q12 * (x - x0)
        v =  2 * q11 * (x - x0) + 2 * q12 * (y - y0)
        return np.concatenate([(u - u_i), (v - v_i)])

    def fit_params(x, y, u_i, v_i):
        x0_init, y0_init = np.mean(x), np.mean(y)
        q11_init, q12_init, q22_init = 1.0, 0.0, 1.0  # Initial guesses
        params_init = [x0_init, y0_init, q11_init, q12_init, q22_init]
        result = least_squares(residuals, params_init, args=(x, y, u_i, v_i))
        return result.x 

    x0, y0, q11, q12, q22 = fit_params(xi, yi, ui, vi)

    w = 2*(q11 + q22)

    Q = np.array([[q11, q12], [q12, q22]])

    Rc_opt, psi0_opt = espra_Rc(xi, yi, ui, vi, x0, y0, q11, q12, q22)
    
    return x0, y0, w, Q, Rc_opt, psi0_opt

def fit_Rc_VN1(x1, y1, v1, Q11, Q12, Q22, x0, y0,
               psi0_0=30.0, Rc0=10.0,
               bounds_scale=(-np.inf, np.inf), bounds_Rc=(1e-6, np.inf),
               max_expo=700.0, penalty=1e20):
    from scipy.optimize import minimize
    y_T = y1[0]

    def obj(p):
        psi0, Rc = p
        # enforce Rc>0
        if Rc <= 0:
            return penalty

        # q11 = - Rc**2 / psi0 * Q11
        # q12 = - Rc**2 / psi0 * Q12
        # q22 = - Rc**2 / psi0 * Q22

        q11 = - psi0 / Rc**2 * Q11
        q12 = - psi0 / Rc**2 * Q12
        q22 = - psi0 / Rc**2 * Q22

        # un-scaled quad for exponent
        C = q11
        B = -2*q11*x0 + 2*q12*(y_T - y0)
        E = q11*x0**2 - 2*q12*x0*(y_T - y0) + q22*(y_T - y0)**2

        expo = -(C*x1**2 + B*x1 + E) / Rc**2

        # **reject any candidate that would overflow exp()**
        if np.any(expo > max_expo):
            return penalty

        VN = -(2*C*x1 + B) / Rc**2 * np.exp(expo) * psi0

        # catch any remaining NaN/Inf
        if not np.isfinite(VN).all():
            return penalty

        return np.sum((v1 - VN)**2)

    p0 = np.array([psi0_0, Rc0])
    res = minimize(obj, p0,
                   bounds=[bounds_scale, bounds_Rc],
                   method='L-BFGS-B')
    if not res.success:
        return np.nan, np.nan
        # raise RuntimeError("Optimization failed: " + res.message)
        
    return res.x  # (psi0_opt, Rc_opt)



def fit_Rc_VN2(x2, y2, u2, Q11, Q12, Q22, x0, y0,
           psi0_0=30.0, Rc0=10.0,
           bounds_psi0=(-np.inf, np.inf), bounds_Rc=(1e-6, np.inf)):
    from scipy.optimize import minimize

    # q11, q22 = np.abs(q11), np.abs(q22)

    # relative transect offset
    x_T = x2[0]

    def obj(p):
        psi0, Rc = p
        if Rc <= 0:
            return np.inf

        # q11 = - Rc**2 / psi0 * Q11
        # q12 = - Rc**2 / psi0 * Q12
        # q22 = - Rc**2 / psi0 * Q22

        q11 = - psi0 / Rc**2 * Q11
        q12 = - psi0 / Rc**2 * Q12
        q22 = - psi0 / Rc**2 * Q22

        # build un‐scaled quadratic
        C = Q22
        B = 2*y0*Q22 - 2*Q12*(x_T-x0)
        E = Q11*(x_T-x0)**2 - 2*Q12*(x_T-x0)*y0 + Q22*y0**2

        # modelled V_N
        expo = -(C*y2**2 - B*y2 + E) / Rc**2
        VN = -(-2*C*y2 + B) / Rc**2 * np.exp(expo) * psi0

        return np.sum((u2 - VN)**2)

    p0 = np.array([psi0_0, Rc0])

    res = minimize(obj, p0,
                   bounds=[bounds_psi0, bounds_Rc],
                   method='L-BFGS-B')
    if not res.success:
        raise RuntimeError("Optimization failed: " + res.message)

    return res.x  # psi0_opt, Rc_opt

def espra_Rc(xi, yi, ui, vi, x0, y0, Q11, Q12, Q22):
    from scipy.optimize import least_squares
    if np.any(np.isnan(ui)):
        return np.nan, np.nan

    # q11, q22 = np.abs(q11), np.abs(q22)

    def residuals(params, x, y, u_i, v_i):
        Rc, psi0 = params

        # q11 = - Rc**2 / psi0 * Q11
        # q12 = - Rc**2 / psi0 * Q12
        # q22 = - Rc**2 / psi0 * Q22

        q11 = - psi0 / Rc**2 * Q11
        q12 = - psi0 / Rc**2 * Q12
        q22 = - psi0 / Rc**2 * Q22

        phi = Q11*(x - x0)**2 + 2*Q12*(x - x0)*(y - y0) + Q22*(y - y0)**2
        phi_x = 2*Q11*(x - x0) + 2*Q12*(y - y0)
        phi_y = 2*Q22*(y - y0) + 2*Q12*(x - x0)

        factor = - 1 / Rc**2
        exp_term = np.exp(factor * phi)

        u = -factor * phi_y * exp_term * psi0
        v = factor * phi_x * exp_term * psi0

        return np.concatenate([(u - u_i), (v - v_i)])

    # Initial guesses: Rc=10, psi0=1
    params_init = [10.0, 30.0]
    bounds_lower = [1e-6, -1000]   # Rc ≥ 1, scale ≥ 0.01
    bounds_upper = [np.inf, 1000] # Rc ≤ 20, scale ≤ 100

    result = least_squares(residuals, params_init, bounds=(bounds_lower, bounds_upper), args=(xi, yi, ui, vi))
    Rc_opt, psi0_opt = result.x
    return Rc_opt, psi0_opt

def gaussian_vel_reconstruction(x0, y0, Q11, Q12, Q22, Rc, psi0):

    # q11, q22 = np.abs(q11), np.abs(q22)

    # q11 = - Rc**2 / psi0 * Q11
    # q12 = - Rc**2 / psi0 * Q12
    # q22 = - Rc**2 / psi0 * Q22

    q11 = - psi0 / Rc**2 * Q11
    q12 = - psi0 / Rc**2 * Q12
    q22 = - psi0 / Rc**2 * Q22
    
    width = 200
    x = np.linspace(x0-width, x0+width, 51)
    y = np.linspace(y0-width, y0+width, 51)
    X, Y = np.meshgrid(x, y)
    
    dx, dy = X - x0, Y - y0
    
    phi   = q11*dx**2 + 2*q12*dx*dy + q22*dy**2
    phi_x = 2*q11*dx  + 2*q12*dy
    phi_y = 2*q22*dy  + 2*q12*dx
    
    # 5) build Gaussian streamfunction with that R
    exp_term = np.exp(-phi / Rc**2)
    psi_x = -phi_x / Rc**2 * exp_term
    psi_y = -phi_y / Rc**2 * exp_term
    
    u = -psi_y * psi0
    v =  psi_x * psi0

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














