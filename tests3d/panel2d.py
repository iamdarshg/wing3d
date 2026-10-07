"""Independent 2D Hess-Smith panel method (constant source+vortex) for ground truth."""
import numpy as np


def panel_2d(cl_pts, alpha_deg, n_wake=0):
    """cl_pts: closed airfoil loop (M,2) CCW starting at TE. Returns cl."""
    p = np.asarray(cl_pts, dtype=float)
    m = len(p) - 1  # assume last == first
    if not np.allclose(p[0], p[-1]):
        p = np.vstack([p, p[:1]])
        m = len(p) - 1
    a = p[:-1]
    b = p[1:]
    mid = 0.5 * (a + b)
    dx, dy = b[:, 0] - a[:, 0], b[:, 1] - a[:, 1]
    leng = np.hypot(dx, dy)
    tx, ty = dx / leng, dy / leng
    nx, ny = ty, -tx  # outward for CCW? ensure below
    # orient: signed area >0 CCW; outward normal = (ty, -tx) for CCW
    area2 = 0.5 * np.sum(a[:, 0] * b[:, 1] - b[:, 0] * a[:, 1])
    if area2 < 0:
        nx, ny = -nx, -ny
    al = np.radians(alpha_deg)
    V = np.array([np.cos(al), np.sin(al)])
    A = np.zeros((m + 1, m + 1))
    rhs = np.zeros(m + 1)
    for i in range(m):
        for j in range(m):
            # constant-source + constant-vortex influence (Hess & Smith)
            ex, ey = mid[i, 0] - a[j, 0], mid[i, 1] - a[j, 1]
            # local frame of panel j
            e1 = ex * tx[j] + ey * ty[j]
            e2 = -ex * ty[j] + ey * tx[j]
            L = leng[j]
            if i == j:
                u_s, w_s = 0.0, 0.5  # source: normal 1/2
                u_v, w_v = -0.5, 0.0  # vortex: tangential -1/2 (local)
            else:
                import math
                r1 = math.hypot(e1, e2)
                r2 = math.hypot(e1 - L, e2)
                th1 = math.atan2(e2, e1)
                th2 = math.atan2(e2, e1 - L)
                u_s = -0.5 / math.pi * math.log(r2 / max(r1, 1e-300))
                w_s = -(th2 - th1) / (2 * math.pi)
                u_v = -w_s
                w_v = u_s
            # global: velocity = u*t + w*n for source; vortex rotated
            # source velocity in global:
            usx = u_s * tx[j] + w_s * nx[j]
            usy = u_s * ty[j] + w_s * ny[j]
            # vortex velocity (streamfunction conjugate):
            uvx = u_v * tx[j] + w_v * nx[j]
            uvy = u_v * ty[j] + w_v * ny[j]
            A[i, j] = usx * nx[i] + usy * ny[i]  # source col (normal)
            A[i, m] += 0  # vortex strength single unknown accumulated below
        # vortex column: sum over panels of vortex influence
        A[i, m] = 0.0
        for j in range(m):
            ex, ey = mid[i, 0] - a[j, 0], mid[i, 1] - a[j, 1]
            e1 = ex * tx[j] + ey * ty[j]
            e2 = -ex * ty[j] + ey * tx[j]
            L = leng[j]
            if i == j:
                u_v, w_v = -0.5, 0.0
            else:
                import math
                r1 = math.hypot(e1, e2)
                r2 = math.hypot(e1 - L, e2)
                th1 = math.atan2(e2, e1)
                th2 = math.atan2(e2, e1 - L)
                u_s = -0.5 / math.pi * math.log(r2 / max(r1, 1e-300))
                w_s = -(th2 - th1) / (2 * math.pi)
                u_v = -w_s
                w_v = u_s
            uvx = u_v * tx[j] + w_v * nx[j]
            uvy = u_v * ty[j] + w_v * ny[j]
            A[i, m] += uvx * nx[i] + uvy * ny[i]
        rhs[i] = -(V[0] * nx[i] + V[1] * ny[i])
    # Kutta: tangential velocity equal at first/last panel (TE)
    for j in range(m):
        for (i, sgn) in ((0, 1.0), (m - 1, -1.0)):
            ex, ey = mid[i, 0] - a[j, 0], mid[i, 1] - a[j, 1]
            e1 = ex * tx[j] + ey * ty[j]
            e2 = -ex * ty[j] + ey * tx[j]
            L = leng[j]
            if i == j:
                u_s, w_s = 0.0, 0.5
                u_v, w_v = -0.5, 0.0
            else:
                import math
                r1 = math.hypot(e1, e2)
                r2 = math.hypot(e1 - L, e2)
                th1 = math.atan2(e2, e1)
                th2 = math.atan2(e2, e1 - L)
                u_s = -0.5 / math.pi * math.log(r2 / max(r1, 1e-300))
                w_s = -(th2 - th1) / (2 * math.pi)
                u_v = -w_s
                w_v = u_s
            usx = u_s * tx[j] + w_s * nx[j]
            usy = u_s * ty[j] + w_s * ny[j]
            uvx = u_v * tx[j] + w_v * nx[j]
            uvy = u_v * ty[j] + w_v * ny[j]
            A[m, j] += sgn * (usx * tx[i] + usy * ty[i])
        # vortex col
        A[m, m] += 0
    # vortex column for kutta row
    for j in range(m):
        for (i, sgn) in ((0, 1.0), (m - 1, -1.0)):
            ex, ey = mid[i, 0] - a[j, 0], mid[i, 1] - a[j, 1]
            e1 = ex * tx[j] + ey * ty[j]
            e2 = -ex * ty[j] + ey * tx[j]
            L = leng[j]
            if i == j:
                u_v, w_v = -0.5, 0.0
            else:
                import math
                r1 = math.hypot(e1, e2)
                r2 = math.hypot(e1 - L, e2)
                th1 = math.atan2(e2, e1)
                th2 = math.atan2(e2, e1 - L)
                u_s = -0.5 / math.pi * math.log(r2 / max(r1, 1e-300))
                w_s = -(th2 - th1) / (2 * math.pi)
                u_v = -w_s
                w_v = u_s
            uvx = u_v * tx[j] + w_v * nx[j]
            uvy = u_v * ty[j] + w_v * ny[j]
            A[m, m] += sgn * (uvx * tx[i] + uvy * ty[i])
    rhs[m] = -((V[0] * tx[0] + V[1] * ty[0]) - (V[0] * tx[m - 1] + V[1] * ty[m - 1]))
    x = np.linalg.solve(A, rhs)
    gamma = x[m]
    return 2 * gamma  # cl = 2*Gamma/(V c), c=1, V=1
