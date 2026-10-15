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
    xq = np.linspace(x1, x2, n)
    # group by chord position, upper/lower by thickness coordinate
    binw = chord / (n - 1)
    idx = np.clip(np.round((P[:, ch] - x1) / binw).astype(int), 0, n - 1)
    # classify per bin by median split: upper/lower x-coordinates interleave
    # (surfaces have different x-distributions), so a bin can catch only one
    # side -- a naive min/max over the bin would then assign the upper value
    # to yl. Median-split keeps the two surfaces separate; missing sides are
    # filled by interpolation below.
    yu = np.full(n, np.nan)
    yl = np.full(n, np.nan)
    z = P[:, th]
    for k in range(n):
        sel = idx == k
        if not np.any(sel):
            continue
        zz = z[sel]
        zmed = np.median(zz)
        up = zz[zz > zmed]
        lo = zz[zz < zmed]
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
    # camber can be non-zero even for a ROTATED section: the mid-chord
    # mean height is the rotation offset; report the shape part only
    mean0 = mean - mean.mean()
    return {
        'x': (xq - x1) / chord,
        'yu': (yu - mean.mean()) / chord,
        'yl': (yl - mean.mean()) / chord,
        'camber_max': float(np.max(np.abs(mean0)) / chord),
        'thick_max': float(np.max(half) * 2 / chord),
        'span_axis': thin,
        'station': float(station),
        'chord': float(chord),
        'mesh_extent': [float(e) for e in ext],
    }


def run_tsd_on_stl(path, mach=0.9, alpha_deg=1.25, span=1.6, epochs=8,
                   station=None):
    """Full TSD run on an STL. Returns a result dict."""
    from wing3d.geometry import read_stl
    from wing3d.tsd import section_profile_xy
    from tsd_m_sup import run_sup
    from wing3d.transonic import m_crit, predict
    mesh = read_stl(path)
    sec = extract_section(mesh, station=station)
    prof = section_profile_xy(sec['x'], sec['yu'], sec['yl'])
    # use the campaign-validated path (panel circulation + supersonic BCs)
    free_k = mach < 0.95
    L, xc, ix, j0 = run_sup(mach, span=span, epochs=epochs,
                             free_k=free_k, profile=prof)
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
        'of_anchor_CL': (r['CL'] if isinstance(r.get('CL'), float)
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
