"""3D viscous-inviscid coupling: streamline IBL + transpiration (blowing).

Outer loop (Carter/Veldman-style quasi-simultaneous spirit, relaxed):
  inviscid solve -> trace surface streamlines -> IBL (Thwaites/Michel/Head)
  -> displacement blowing velocity -> transpiration source in next solve.
Converges CL/CD including friction drag + displacement decambering.
"""
import numpy as np


def trace_streamlines(mesh, vel, ds_factor=0.6, max_steps=400):
    """Follow surface velocity from stagnation seeds. Returns list of
    streamlines; each is a dict with panels, s, Ue arrays."""
    from scipy.spatial import cKDTree
    C = mesh.centroid
    tree = cKDTree(C)
    spd = np.linalg.norm(vel, axis=1)
    n = mesh.npanels
    covered = np.zeros(n, dtype=bool)
    lines = []
    sizes = np.sqrt(mesh.area)
    order = np.argsort(spd)  # start at stagnation, expand outward
    for seed in order:
        if covered[seed] and len(lines) > 4:
            # still allow a few seeds for coverage, then rely on fill
            pass
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
                # end of reachable surface (edge/wake/TE)
                panels.append(int(j)) if dd <= 3.0 * sizes[i] else None
                if dd <= 3.0 * sizes[i]:
                    ss.append(s + ds)
                    ues.append(max(np.linalg.norm(vel[int(j)]), 1e-9))
                    covered[int(j)] = True
                break
            s += ds
            i = int(j)
        if len(panels) >= 5:
            lines.append({'panels': np.array(panels),
                          's': np.array(ss), 'Ue': np.array(ues)})
        if covered.all():
            break
    return lines, covered


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
        lines, _ = trace_streamlines(mesh, res['vel'])
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
