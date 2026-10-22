"""CLI: run the TSD full-potential solver on an STL geometry.

Usage:
    python -m wing3d.tsdcli <file.stl> [--mach 0.9] [--alpha 1.25]
                             [--span 1.6] [--epochs 8] [--json out.json]

Reads the STL (binary or ASCII), auto-detects the extrusion axis, slices
the mid-span airfoil section, builds the section profile (thickness +
camber line) and runs the full-potential TSD with the supersonic BCs.
Prints CL/CD, section cl, camber magnitude and thickness ratio.
"""
import argparse
import json
import os
import sys

import numpy as np

# campaign-validated solver path: tests3d/ is repo-adjacent, not installed
_TESTS3D = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'tests3d')
if _TESTS3D not in sys.path:
    sys.path.insert(0, _TESTS3D)


def extract_section(mesh, station=None, n=400):
    """Slice an airfoil section out of an STL mesh.

    Auto-detects the thin (extrusion/span) axis and slices perpendicular
    to it at `station` (default: the geometric middle of the mesh extent).
    Returns dict(x, yu, yl, camber_max, thick_max, span_axis, station).
    """
    V = np.asarray(mesh.verts, dtype=float)
    ext = V.max(axis=0) - V.min(axis=0)
    # section plane: the section repeats along the LONGEST axis
    # (span for a wing, extrusion length for an airfoil slab)
    thin = int(np.argmax(ext))
    other = [a for a in range(3) if a != thin]
    ch = other[int(np.argmax([ext[a] for a in other]))]   # chord axis
    th = other[int(np.argmin([ext[a] for a in other]))]   # thickness axis
    lo, hi = V[:, thin].min(), V[:, thin].max()
    if station is None:
        station = 0.5 * (lo + hi)
    tol = max(1e-6, 0.02 * (hi - lo))   # slab: 2% of span; wing: wider
    m = np.abs(V[:, thin] - station) <= tol
    if m.sum() < 8:
        # fall back: take the nearest slice that has points
        order = np.argsort(np.abs(V[:, thin] - station))
        m = np.zeros(len(V), dtype=bool)
        m[order[:max(8, len(V) // 50)]] = True
    P = V[m]
    x1, x2 = P[:, ch].min(), P[:, ch].max()
    chord = max(x2 - x1, 1e-12)
    z = P[:, th]
    xe = P[:, ch]
    # pass 1: coarse uniform bins give a robust mean-line reference (this
    # is the binning verified against the campaign tables during geo
    # reconciliation). Pass 2 reuses it to classify upper/lower instead of
    # a per-bin median split: thin cosine bins can catch points of only
    # one surface, and a median split would then misassign them.
    n0 = 200
    xq0 = np.linspace(x1, x2, n0)
    bw0 = chord / (n0 - 1)
    idx0 = np.clip(np.round((xe - x1) / bw0).astype(int), 0, n0 - 1)
    mu0 = np.full(n0, np.nan)
    ml0 = np.full(n0, np.nan)
    for k in range(n0):
        zz = z[idx0 == k]
        if len(zz) == 0:
            continue
        zm = np.median(zz)
        u = zz[zz > zm]
        l = zz[zz < zm]
        if len(u) == 0 and len(l) == 0:
            mu0[k] = ml0[k] = zm          # tip bin: surfaces meet at zm
            continue
        if len(u):
            mu0[k] = u.max()
        if len(l):
            ml0[k] = l.min()
    ok0 = np.isfinite(mu0) & np.isfinite(ml0)
    if ok0.sum() < 3:
        raise ValueError('could not slice a closed section from the STL')
    ref0 = 0.5 * (mu0[ok0] + ml0[ok0])
    xr0 = xq0[ok0]
    # fill gaps, then a short median filter rejects isolated bad bins
    xg0 = np.linspace(x1, x2, n0)
    ref0 = np.interp(xg0, xr0, ref0)
    from scipy.signal import medfilt
    ref0 = medfilt(ref0, 5)
    ref = np.interp(xe, xg0, ref0)
    # pass 2: cosine-clustered output stations (dense at LE/TE like the
    # mesher), split at the local mean-line reference: the mean line
    # separates the two surfaces exactly by construction.
    beta = np.linspace(0.0, np.pi, n)
    xq = x1 + 0.5 * (1.0 - np.cos(beta)) * chord
    edges = np.empty(n + 1)
    edges[1:-1] = 0.5 * (xq[1:] + xq[:-1])
    edges[0] = x1 - 1e-9
    edges[-1] = x2 + 1e-9
    yu = np.full(n, np.nan)
    yl = np.full(n, np.nan)
    for k in range(n):
        sel = (xe >= edges[k]) & (xe < edges[k + 1])
        if not np.any(sel):
            continue
        zz = z[sel]
        zr = ref[sel]
        up = zz[zz >= zr]
        lo = zz[zz <= zr]
        if len(up):
            yu[k] = up.max()
        if len(lo):
            yl[k] = lo.min()
    have_u = np.isfinite(yu)
    have_l = np.isfinite(yl)
    if not have_u.any() or not have_l.any():
        raise ValueError('could not slice a closed section from the STL')
    yu = np.interp(xq, xq[have_u], yu[have_u])
    yl = np.interp(xq, xq[have_l], yl[have_l])
    mean = 0.5 * (yu + yl)
    half = 0.5 * (yu - yl)
    # camber metric: peak-to-peak of the mean line (datum-independent).
    # NOTE: do NOT subtract the mean-line average here -- that is a
    # constant z-shift (a translation, not a de-rotation); it both
    # understates the camber of aft-loaded sections (0.0127 -> 0.0074 for
    # the RAE variant) and changes the section passed to the solver.
    return {
        'x': (xq - x1) / chord,
        'yu': yu / chord,
        'yl': yl / chord,
        'camber_max': float(np.max(mean) - np.min(mean)) / chord,
        'thick_max': float(np.max(half) * 2 / chord),
        'span_axis': thin,
        'station': float(station),
        'chord': float(chord),
        'mesh_extent': [float(e) for e in ext],
    }


def run_tsd_on_stl(path, mach=0.9, alpha_deg=1.25, span=1.6, epochs=8,
                   station=None):
    """Full TSD run on an STL. Returns a result dict."""
    from wing3d.geometry import read_stl, section_loop
    from wing3d.tsd import section_profile_xy
    from tsd_m_sup import run_sup
    from wing3d.transonic import m_crit, predict
    mesh = read_stl(path)
    sec = extract_section(mesh, station=station)
    # bin-level slice noise (~1e-4 in z) must not leak into the wall BC,
    # whose camber/thickness SLOPES feed the solver directly (unsmoothed
    # slopes shifted the transonic branch: CL 0.25 -> 0.16 at M0.85).
    # Ablation knob (audit 2026-10-22, savgol flagged SHAKY/tuned):
    # WING3D_RAW_SECTION=1 -> raw slopes, no savgol.
    if os.environ.get('WING3D_RAW_SECTION', '') == '1':
        yu_s = np.array(sec['yu'])
        yl_s = np.array(sec['yl'])
    else:
        from scipy.signal import savgol_filter
        w = min(41, len(sec['x']) // 4 * 2 + 1)
        yu_s = savgol_filter(sec['yu'], w, 3)
        yl_s = savgol_filter(sec['yl'], w, 3)
        # savgol edge handling extrapolates a polynomial (nose overshoot);
        # pin the endpoints to the raw slice values instead
        yu_s[0], yu_s[-1] = sec['yu'][0], sec['yu'][-1]
        yl_s[0], yl_s[-1] = sec['yl'][0], sec['yl'][-1]
    prof = section_profile_xy(sec['x'], yu_s, yl_s)
    sec_loop = section_loop(sec['x'], yu_s, sec['x'], yl_s, 24)
    # campaign-validated path: fixed k=1 camber imprint + supersonic BCs.
    # (the free-k closure was diagnosed non-viable -- k collapse / spurious
    # TE unloading; see tests3d kscan/kroot diagnostics)
    free_k = False
    L, xc, ix, j0 = run_sup(mach, span=span, epochs=epochs,
                             free_k=free_k, profile=prof, section=sec_loop)
    xg = np.linspace(0.02, 0.95, 40)
    cp0min = float(np.min(prof['slope'](xg))) * 3.0   # panel-suction proxy
    r = predict(mach, alpha_deg)
    out = {
        'stl': os.path.abspath(path),
        'mach': float(mach),
        'alpha_deg': float(alpha_deg),
        'npanels': int(mesh.npanels),
        'nverts': int(len(mesh.verts)),
        'thick_max': sec['thick_max'],
        'camber_max': sec['camber_max'],
        'mcrit': float(m_crit(min(cp0min, -0.1))),
        'CL': float(L['CL']),
        'CD': float(L['CD']),
        'section_cl': float(np.trapezoid(
            (L['cpl'][ix, j0] - L['cpu'][ix, j0]), xc)),
        'surrogate_CL': (r['CL'] if isinstance(r.get('CL'), float)
                         and r['CL'] == r['CL'] else None),
    }
    return out


def main(argv=None):
    p = argparse.ArgumentParser(
        prog='python -m wing3d.tsdcli',
        description='Run the TSD full-potential solver on an STL geometry.')
    p.add_argument('stl', help='path to the STL file (binary or ASCII)')
    p.add_argument('--mach', type=float, default=0.9,
                   help='freestream Mach (default 0.9)')
    p.add_argument('--alpha', type=float, default=1.25,
                   help='angle of attack in degrees (default 1.25)')
    p.add_argument('--span', type=float, default=1.6,
                   help='TSD wing span in chords (default 1.6)')
    p.add_argument('--epochs', type=int, default=8,
                   help='solver iterations (default 8)')
    p.add_argument('--station', type=float, default=None,
                   help='slice station along the thin axis (default: mid)')
    p.add_argument('--json', type=str, default=None,
                   help='optional output JSON path')
    a = p.parse_args(argv)
    if not os.path.exists(a.stl):
        print('error: STL not found: %s' % a.stl, file=sys.stderr)
        return 2
    r = run_tsd_on_stl(a.stl, mach=a.mach, alpha_deg=a.alpha,
                       span=a.span, epochs=a.epochs, station=a.station)
    print('TSD full-potential result for %s' % r['stl'])
    print('  M=%.3f  alpha=%.2fdeg  span=%.1f chords' % (
        r['mach'], r['alpha_deg'], a.span))
    print('  mesh: %d panels / %d vertices' % (r['npanels'], r['nverts']))
    print('  section: t/c=%.3f  camber=%.4f  Mcrit=%.2f' % (
        r['thick_max'], r['camber_max'], r['mcrit']))
    print('  CL=%.4f  section_cl=%.4f  CD=%.4f' % (
        r['CL'], r['section_cl'], r['CD']))
    if a.json:
        with open(a.json, 'w') as f:
            json.dump(r, f, indent=1)
        print('  saved %s' % a.json)
    return 0


if __name__ == '__main__':
    sys.exit(main())
