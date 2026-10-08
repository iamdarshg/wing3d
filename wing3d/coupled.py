"""3D viscous-inviscid coupling: streamline IBL + transpiration (blowing).

Outer loop (Carter/Veldman-style quasi-simultaneous spirit, relaxed):
  inviscid solve -> trace surface streamlines -> IBL (Thwaites/Michel/Head)
  -> displacement blowing velocity -> transpiration source in next solve.
Converges CL/CD including friction drag + displacement decambering.
"""
import numpy as np


def trace_streamlines(mesh, vel, ds_factor=0.6, max_steps=400):
    """Follow surface velocity from STAGNATION seeds (full LE->TE tracks).

    Seeds are low-speed (stagnation/attachment) panels; each streamline
    marches downstream with the local velocity until it leaves the body,
    stagnates, or repeats. This gives proper boundary-layer histories
    (Thwaites must integrate from stagnation, not mid-chord fragments).
    Returns list of streamlines + coverage mask.
    """
    from scipy.spatial import cKDTree
    C = mesh.centroid
    tree = cKDTree(C)
    spd = np.linalg.norm(vel, axis=1)
    n = mesh.npanels
    covered = np.zeros(n, dtype=bool)
    lines = []
    sizes = np.sqrt(mesh.area)
    # stagnation threshold: slowest 5% (attachment line/stagnation region)
    thresh = np.quantile(spd, 0.05)
    seeds = np.where(spd <= thresh)[0]
    # thin seeds spatially (avoid redundant neighbors)
    kept = []
    for s in seeds:
        if all(np.linalg.norm(C[s] - C[k]) > 2.5 * sizes[s] for k in kept):
            kept.append(int(s))
    for seed in kept:
        panels, ss, ues = [], [], []
        i = int(seed)
        s = 0.0
        seen = set()
        for _ in range(max_steps):
            if i in seen:
                break
            seen.add(i)
            v = vel[i]
            sp = np.linalg.norm(v)
            if sp < 1e-9:
                break
            panels.append(i)
            ss.append(s)
            ues.append(sp)
            covered[i] = True
            d = v / sp
            ds = ds_factor * sizes[i]
            pnext = C[i] + d * ds
            dd, j = tree.query(pnext)
            if dd > 3.0 * sizes[i] or int(j) in seen:
                if dd <= 3.0 * sizes[i]:
                    panels.append(int(j))
                    ss.append(s + ds)
                    ues.append(max(np.linalg.norm(vel[int(j)]), 1e-9))
                    covered[int(j)] = True
                break
            s += ds
            i = int(j)
        if len(panels) >= 5:
            lines.append({'panels': np.array(panels),
                          's': np.array(ss), 'Ue': np.array(ues)})
    # fill gaps: seed uncovered panels (short tracks for coverage only)
    for seed in np.argsort(spd):
        if covered[int(seed)]:
            continue
        panels, ss, ues = [], [], []
        i = int(seed)
        s = 0.0
        seen = set()
        for _ in range(60):
            if i in seen:
                break
            seen.add(i)
            v = vel[i]
            sp = np.linalg.norm(v)
            if sp < 1e-9:
                break
            panels.append(i)
            ss.append(s)
            ues.append(sp)
            covered[i] = True
            d = v / sp
            ds = ds_factor * sizes[i]
            dd, j = tree.query(C[i] + d * ds)
            if dd > 3.0 * sizes[i] or int(j) in seen:
                break
            s += ds
            i = int(j)
        if len(panels) >= 5:
            lines.append({'panels': np.array(panels),
                          's': np.array(ss), 'Ue': np.array(ues),
                          'fragment': True})
        if covered.all():
            break
    return lines, covered


def smooth_surface_field(mesh, F, passes=1):
    """Area-weighted neighbor-average smoothing of a surface field.

    Removes panel-discretization singularities (LE/corner velocity
    spikes) while preserving the smooth physical field. Standard
    de-singularization before boundary-layer tracing: the IBL
    integrates dUe/ds from stagnation, where a single spiked panel
    corrupts the entire downstream history (Thwaites lambda blowup).
    """
    F = np.asarray(F, dtype=float).copy()
    n = mesh.npanels
    # panel adjacency via shared verts
    vert_to_pan = {}
    for i, f in enumerate(mesh.faces):
        for v in f:
            vert_to_pan.setdefault(int(v), []).append(i)
    nbr = [set() for _ in range(n)]
    for lst in vert_to_pan.values():
        for i in lst:
            nbr[i].update(lst)
    nbr = [sorted(s - {i}) for i, s in enumerate(nbr)]
    w = np.asarray(mesh.area, dtype=float)
    for _ in range(passes):
        G = F.copy() if F.ndim == 1 else F.copy()
        if F.ndim == 1:
            for i in range(n):
                js = nbr[i]
                if js:
                    tot = w[i] + w[js].sum()
                    G[i] = (w[i] * F[i] + (w[js] * F[js]).sum()) / tot
        else:
            for i in range(n):
                js = nbr[i]
                if js:
                    tot = w[i] + w[js].sum()
                    G[i] = (w[i] * F[i] + (w[js, None] * F[js]).sum(
                        axis=0)) / tot
        F = G
    return F


def panel_ibl_map(mesh, lines, nu):
    """Run IBL per streamline; map Cf/dstar/d(Ue d*)/ds back to panels."""
    from .ibl import solve_streamline
    n = mesh.npanels
    cf = np.zeros(n)
    dstar = np.zeros(n)
    dfds = np.zeros(n)
    sep = np.zeros(n, dtype=bool)
    count = np.zeros(n)
    nsep_lines = 0
    for ln in lines:
        r = solve_streamline(ln['s'], ln['Ue'], nu)
        if r['separated']:
            nsep_lines += 1
        F = r['Ue'] * r['dstar']
        # smooth + differentiate along the streamline (robust to noise)
        w = min(5, len(F) // 2 * 2 + 1)
        if w >= 3 and len(F) >= 3:
            ker = np.ones(w) / w
            Fs = np.convolve(np.pad(F, w // 2, mode='edge'), ker,
                             mode='valid')
        else:
            Fs = F
        dF = np.gradient(Fs, ln['s'])
        dF = np.clip(dF, -0.05, 0.05)
        for k, p in enumerate(ln['panels']):
            if k >= len(r['cf']):
                break
            cf[p] += r['cf'][k]
            dstar[p] += r['dstar'][k]
            dfds[p] += dF[k]
            sep[p] = sep[p] or (r['i_sep'] >= 0 and k >= r['i_sep'])
            count[p] += 1
    m = count > 0
    cf[m] /= count[m]
    dstar[m] /= count[m]
    dfds[m] /= count[m]
    # fill uncovered from nearest covered
    if not m.all():
        from scipy.spatial import cKDTree
        tree = cKDTree(mesh.centroid[m])
        _, jj = tree.query(mesh.centroid[~m])
        idx = np.where(m)[0][jj]
        cf[~m] = cf[idx]
        dstar[~m] = dstar[idx]
        dfds[~m] = dfds[idx]
    return {'cf': cf, 'dstar': dstar, 'dfds': dfds, 'sep': sep,
            'sep_fraction': float(sep.mean()),
            'sep_lines': nsep_lines, 'n_lines': len(lines)}


def blowing_rhs(mesh, Brow, vel, bl, relax=1.0):
    """Transpiration source from streamline displacement growth.

    Uses Vn = d(Ue*d*)/ds computed along smooth streamlines (bl['dfds']),
    mapped to panels. Returns extra_rhs contribution.
    """
    vn = np.clip(np.asarray(bl['dfds'], dtype=float), -0.05, 0.05)
    sigma_t = 2.0 * vn  # blowing velocity -> transpiration source
    return -(Brow @ sigma_t), vn


def solve_coupled(mesh, wakes, vinf, Re, Lref=1.0, itmax=12, relax=0.3,
                  kutta_mode='doublet', hist=None):
    """Full viscous-coupled solve. Returns (res, forces, info)."""
    from .solver import solve, assemble
    from .forces import pressure_forces
    V = np.asarray(vinf, dtype=float)
    vmag = np.linalg.norm(V)
    nu = vmag * Lref / Re
    A, Brow = assemble(mesh, wakes, kutta_mode=kutta_mode)
    # sref: planform-ish (caller overrides via info)
    extra = np.zeros(mesh.npanels + len(wakes))
    if hist is None:
        hist = []
    bl = None
    n = mesh.npanels
    for it in range(itmax):
        res = solve(mesh, wakes, V, A=A, Brow=Brow, extra_rhs=extra[:n],
                    kutta_mode=kutta_mode)
        f = pressure_forces(mesh, res['cp'], V, sref=1.0)
        # friction from current BL (previous iter) for history
        cdf = 0.0 if bl is None else float(
            (bl['cf'] * mesh.area).sum() / 1.0)
        hist.append({'CL': f['CL'], 'CD': f['CDp'] + cdf, 'it': it})
        vmag_f = np.linalg.norm(res['vel'], axis=1)
        vmag_s = smooth_surface_field(mesh, vmag_f, passes=1)
        # rebuild vectors with smoothed magnitude, original direction
        dirn = res['vel'] / np.maximum(vmag_f, 1e-12)[:, None]
        lines, _ = trace_streamlines(mesh, dirn * vmag_s[:, None])
        bl = panel_ibl_map(mesh, lines, nu)
        rhs_new, _ = blowing_rhs(mesh, Brow, res['vel'], bl)
        extra_new = np.zeros_like(extra)
        extra_new[:mesh.npanels] = rhs_new
        extra = (1 - relax) * extra + relax * extra_new
        if it >= 2:
            dcl = abs(hist[-1]['CL'] - hist[-2]['CL'])
            dcd = abs(hist[-1]['CD'] - hist[-2]['CD'])
            if dcl < 1e-4 and dcd < 1e-4:
                break
    # final friction drag along drag dir
    from .forces import pressure_forces as pf
    f = pf(mesh, res['cp'], V, sref=1.0)
    vhat = V / vmag
    lift = np.array([0., 0., 1.])
    lift = lift - (lift @ vhat) * vhat
    lift /= max(np.linalg.norm(lift), 1e-300)
    side = np.cross(vhat, lift)
    drag = np.cross(lift, side)
    cdf_vec = (bl['cf'][:, None] * (-res['vel'] / np.maximum(
        np.linalg.norm(res['vel'], axis=1, keepdims=True), 1e-9))
        * mesh.area[:, None]).sum(axis=0)
    info = {'history': hist, 'bl': bl, 'cdf_vec': cdf_vec,
            'drag_dir': drag, 'iterations': it + 1, 'Re': Re}
    return res, f, info
