"""Complex 3D test shapes: lofted bodies, Ahmed-like car, F-16-like config.

All single watertight surfaces (lofted) except wing/fuselage junctions on
the multi-surface configs, which use overlapping closed solids
(documented approximation). Waked builders return (mesh, wakes, sref)
with per-surface TE wakes for lifting panels.
"""
import numpy as np
from .geometry import Mesh
from .primitives import merge, transform
from .geometry import build_wing


def rounded_rect(hw, hh, y0=0.0, z0=0.0, n=12, corner=0.3):
    """Rounded-rectangle loop (y,z) centered at (y0,z0)."""
    pts = []
    # param each corner arc + edges; simple superellipse instead:
    ts = np.linspace(0, 2 * np.pi, n, endpoint=False)
    p = 3.0
    for t in ts:
        c, s = np.cos(t), np.sin(t)
        y = hw * np.sign(c) * abs(c) ** (2 / p) + y0
        z = hh * np.sign(s) * abs(s) ** (2 / p) + z0
        pts.append([y, z])
    return np.array(pts)


def loft(stations):
    """Loft quads through sections. stations = [(x, loop(N,2) as (y,z))].

    Caps both ends with fans. All loops must have equal length.
    """
    n = len(stations[0][1])
    verts = []
    index = {}
    for j, (x, loop) in enumerate(stations):
        for i, (y, z) in enumerate(loop):
            index[(i, j)] = len(verts)
            verts.append([x, y, z])
    faces = []
    ns = len(stations)
    for j in range(ns - 1):
        for i in range(n):
            i2 = (i + 1) % n
            faces.append([index[(i, j)], index[(i2, j)],
                          index[(i2, j + 1)], index[(i, j + 1)]])
    for j in (0, ns - 1):
        pts = np.array([verts[index[(i, j)]] for i in range(n)])
        ci = len(verts)
        verts.append(pts.mean(axis=0))
        for i in range(n):
            i2 = (i + 1) % n
            faces.append([index[(i, j)], index[(i2, j)], ci])
    return Mesh(np.array(verts), faces)


def build_car(slant_deg=25.0, n_sec=16, n_loop=20):
    """Ahmed-body-like slantback (procedural, meters-ish, x aft).

    Length ~1.04, width 0.39, height 0.29, ground clearance 0.05.
    Nose rounded via shrinking sections; slantback via wedge profile.
    """
    L, W, H, gnd = 1.044, 0.389, 0.288, 0.05
    slant = np.radians(slant_deg)
    # height profile: nose ramp up, roof, slant down, base
    xs = np.linspace(0, L, n_sec)
    stations = []
    for x in xs:
        if x < 0.15:  # rounded nose: grow from small section
            s = 0.35 + 0.65 * (x / 0.15)
            hw, htop = W / 2 * s, gnd + (H - gnd) * s
            hbot = gnd + 0.02 * (1 - s)
        elif x < 0.55:  # forebody/midbox full section
            hw, htop, hbot = W / 2, H, gnd
        elif x < 0.55 + (H - 0.10) / np.tan(slant):  # slant panel
            t = (x - 0.55)
            hw = W / 2
            htop = H - t * np.tan(slant)
            hbot = gnd
        else:  # vertical base
            hw, htop, hbot = W / 2, 0.10 + gnd, gnd
        loop = rounded_rect(hw, (htop - hbot) / 2, 0.0, (htop + hbot) / 2,
                            n_loop)
        stations.append((x, loop))
    return loft(stations)


def build_f16(scale=1.0):
    """Simplified F-16-like config (qualitative geometry, x aft, meters).

    Fuselage loft + canopy bump + clipped-delta wing + h/v tails. Wing and
    tails overlap the fuselage (closed intersecting solids).
    """
    L = 15.0 * scale
    # fuselage stations: x, half-width, center height, half-height
    prof = [(0.00, 0.02, 0.0, 0.02), (0.06, 0.25, 0.0, 0.25),
            (0.15, 0.55, 0.05, 0.55), (0.30, 0.75, 0.1, 0.75),
            (0.45, 0.80, 0.1, 0.80), (0.60, 0.78, 0.1, 0.85),
            (0.75, 0.70, 0.1, 0.80), (0.90, 0.55, 0.1, 0.60),
            (0.97, 0.35, 0.1, 0.40), (1.00, 0.12, 0.1, 0.14)]
    stations = []
    for fx, hw, zc, hh in prof:
        loop = rounded_rect(hw * scale, hh * scale, 0.0, zc * scale, 16)
        stations.append((fx * L, loop))
    fuse = loft(stations)
    # canopy bump (scaled half-ellipsoid merged on top)
    from .primitives import ellipsoid
    canopy = ellipsoid(a=2.2 * scale, b=0.45 * scale, c=0.45 * scale,
                       nlat=10, nlon=16, center=(0.32 * L, 0, 0.75 * scale))
    # main wing: clipped delta-ish (taper + sweep); 0012 (thin sections
    # destabilize the Dirichlet system where solids overlap). Uniform
    # spacing avoids cosine tip-TE slivers on the tapered tip.
    wing = build_wing('0012', span=9.0 * scale, chord=3.2 * scale,
                      taper=0.4, sweep_deg=32.0, n_chord=16, n_span=8,
                      cosine=False)
    wing = transform(wing, translate=(0.45 * L, 0, -0.15 * scale))
    # horizontal tails
    htail = build_wing('0012', span=5.0 * scale, chord=1.8 * scale,
                       taper=0.5, sweep_deg=40.0, n_chord=10, n_span=6,
                       cosine=False)
    htail = transform(htail, translate=(0.88 * L, 0, 0.05 * scale))
    # vertical tail: rotate a wing -90 deg about x-axis (span y -> vertical z)
    vtail = build_wing('0012', span=2.6 * scale, chord=2.0 * scale,
                       taper=0.5, sweep_deg=45.0, n_chord=10, n_span=6,
                       cosine=False)
    v = np.asarray(vtail.verts)
    v3 = np.stack([v[:, 0], -v[:, 2], v[:, 1]], axis=1)
    from .geometry import Mesh as _M
    vtail2 = _M(v3 + np.array([0.86 * L, 0, 1.80 * scale]),
                [list(f) for f in vtail.faces])
    parts = merge([fuse, canopy, wing, htail, vtail2])
    parts.meta = {'span': 9.0 * scale, 'length': L, 'kind': 'f16-like'}
    return parts


def _place_wing(kw, translate):
    """build_wing + translate (verts and TE points together)."""
    from .solver import wing_strips_from_structured, build_wake, WakeStrip
    w = build_wing(**kw)
    t = np.asarray(translate, dtype=float)
    w.verts = np.asarray(w.verts) + t
    meta = w.meta
    meta['te_points'] = [np.asarray(p) + t for p in meta['te_points']]
    w._compute()
    return w


def _vtail_up(kw, translate):
    """Vertical tail: wing rotated span->z, then translated (with TE pts)."""
    from .geometry import Mesh as _M
    w = build_wing(**kw)
    v = np.asarray(w.verts, dtype=float)
    v3 = np.stack([v[:, 0], -v[:, 2], v[:, 1]], axis=1)
    t = np.asarray(translate, dtype=float)
    out = _M(v3 + t, [list(f) for f in w.faces])
    tep = [np.array([p[0], -p[2], p[1]]) + t
           for p in w.meta['te_points']]
    out.meta = dict(w.meta)
    out.meta['te_points'] = tep
    out._compute()
    return out


def _surface_wakes(part, offset, length=4.0, n_panels=3, cut=None):
    """WakeStrip list for one lifting part; ids shifted by offset.

    cut: None | ('y', ycut) | ('z', zcut) -- skip strips whose TE
    midpoint lies inside the body (|y|<ycut or z<zcut). Returns
    (wakes, skipped_global_pids) so forces can mask buried panels.
    """
    from .solver import (wing_strips_from_structured, build_wake)
    meta = part.meta
    strips, panel_id = wing_strips_from_structured(meta['n_loop'],
                                                  meta['n_span'])
    info = []
    skipped = []
    npl = meta['n_loop'] - 1
    for (up, lo, j) in strips:
        pa = np.asarray(meta['te_points'][j])
        pb = np.asarray(meta['te_points'][j + 1])
        mid = 0.5 * (pa + pb)
        drop = False
        if cut is not None:
            ax, lim = cut
            v = mid[1] if ax == 'y' else mid[2]
            drop = (abs(v) < lim) if ax == 'y' else (v < lim)
        if drop:
            for i in range(npl):
                skipped.append(offset + panel_id[(i, j)])
            continue
        info.append((up, lo, pa, pb))
    wakes = build_wake(part, info, length=length, n_panels=n_panels)
    for w in wakes:
        w.upper_te += offset
        w.lower_te += offset
    return wakes, skipped


def _merge_with_wakes(parts, wake_specs, wake_kw=None):
    """Merge bodies; shift WakeStrip ids by part panel offsets.

    wake_specs: list of (part, cut) with cut as in _surface_wakes.
    Records buried-panel mask in merged.meta['skipped'] for forces.
    """
    kw = dict(wake_kw or {})
    kw.pop('cut', None)
    offsets = []
    n = 0
    for p in parts:
        offsets.append(n)
        n += p.npanels
    wakes = []
    skipped = []
    spec = {id(p): c for (p, c) in wake_specs}
    for p, off in zip(parts, offsets):
        if id(p) in spec:
            w, s = _surface_wakes(p, off, cut=spec[id(p)], **kw)
            wakes.extend(w)
            skipped.extend(s)
    merged = merge(parts)
    if merged.npanels != n:
        print('WARNING _merge_with_wakes: %d panels dropped in merge '
              '(wake/skip offsets may shift)' % (n - merged.npanels))
    merged.meta = dict(getattr(parts[0], 'meta', {}))
    merged.meta['skipped'] = np.array(sorted(set(skipped)), dtype=int)
    return merged, wakes


def build_f16_waked(scale=0.2, wake_kw=None):
    """F-16-like with TE wakes on wing + h/v tails. Returns
    (mesh, wakes, sref, info). sref = main-wing planform."""
    L = 15.0 * scale
    prof = [(0.00, 0.02, 0.0, 0.02), (0.06, 0.25, 0.0, 0.25),
            (0.15, 0.55, 0.05, 0.55), (0.30, 0.75, 0.1, 0.75),
            (0.45, 0.80, 0.1, 0.80), (0.60, 0.78, 0.1, 0.85),
            (0.75, 0.70, 0.1, 0.80), (0.90, 0.55, 0.1, 0.60),
            (0.97, 0.35, 0.1, 0.40), (1.00, 0.12, 0.1, 0.14)]
    stations = []
    for fx, hw, zc, hh in prof:
        loop = rounded_rect(hw * scale, hh * scale, 0.0, zc * scale, 16)
        stations.append((fx * L, loop))
    fuse = loft(stations)
    from .primitives import ellipsoid
    canopy = ellipsoid(a=2.2 * scale, b=0.45 * scale, c=0.45 * scale,
                       nlat=10, nlon=16, center=(0.32 * L, 0, 0.75 * scale))
    wing = _place_wing(dict(code='0012', span=9.0 * scale, chord=3.2 * scale,
                            taper=0.4, sweep_deg=32.0, n_chord=16, n_span=8,
                            cosine=False),
                       (0.45 * L, 0, -0.15 * scale))
    htail = _place_wing(dict(code='0012', span=5.0 * scale, chord=1.8 * scale,
                             taper=0.5, sweep_deg=40.0, n_chord=10, n_span=6,
                             cosine=False),
                        (0.88 * L, 0, 0.05 * scale))
    vtail = _vtail_up(dict(code='0012', span=2.6 * scale, chord=2.0 * scale,
                           taper=0.5, sweep_deg=45.0, n_chord=10, n_span=6,
                           cosine=False),
                      (0.86 * L, 0, 1.80 * scale))
    parts = [fuse, canopy, wing, htail, vtail]
    mesh, wakes = _merge_with_wakes(
        parts, [(wing, ('y', 0.85 * scale)), (htail, ('y', 0.60 * scale)),
                (vtail, None)], wake_kw or {})
    sref = float((9.0 * scale) * (3.2 * scale) * (1 + 0.4) / 2)
    mesh.meta = {'span': 9.0 * scale, 'length': L, 'kind': 'f16-like-waked',
                 'skipped': mesh.meta.get('skipped', [])}
    return mesh, wakes, sref, {'span': 9.0 * scale, 'length': L}


def build_a330(scale=1.0, wake_kw=None):
    """Twin-aisle airliner (tube + swept wing + h/v tails).

    Meters at scale=1 (L~59, span~60). Returns (mesh, wakes, sref, info)
    with sref = wing planform. Qualitative loft, overlapping solids.
    """
    L = 59.0 * scale
    R = 2.8 * scale
    # fuselage: nose, constant, tail upsweep (y-z loops, z center varies)
    prof = [(0.00, 0.02, 0.0), (0.04, 0.45, 0.0), (0.10, 0.80, 0.0),
            (0.18, 0.97, 0.0), (0.30, 1.00, 0.0), (0.60, 1.00, 0.0),
            (0.78, 0.95, 0.02), (0.88, 0.75, 0.10), (0.95, 0.45, 0.22),
            (1.00, 0.10, 0.30)]
    stations = []
    for fx, rr, zc in prof:
        loop = rounded_rect(rr * R, rr * R, 0.0, zc * R, 16)
        stations.append((fx * L, loop))
    fuse = loft(stations)
    # wing: tapered, swept 30deg, mounted low at ~42% length
    croot = 11.0 * scale
    wing = _place_wing(dict(code='0012', span=60.0 * scale, chord=croot,
                            taper=0.28, sweep_deg=30.0, n_chord=14,
                            n_span=10, cosine=False),
                       (0.38 * L, 0, -0.9 * scale))
    htail = _place_wing(dict(code='0012', span=19.0 * scale, chord=5.5 * scale,
                             taper=0.35, sweep_deg=35.0, n_chord=10,
                             n_span=6, cosine=False),
                        (0.86 * L, 0, 0.6 * scale))
    vtail = _vtail_up(dict(code='0012', span=9.0 * scale, chord=6.5 * scale,
                           taper=0.4, sweep_deg=40.0, n_chord=10, n_span=6,
                           cosine=False),
                      (0.85 * L, 0, 2.2 * scale))
    parts = [fuse, wing, htail, vtail]
    mesh, wakes = _merge_with_wakes(
        parts, [(wing, ('y', 2.9 * scale)), (htail, ('y', 2.3 * scale)),
                (vtail, ('z', 2.6 * scale))], wake_kw or {})
    sref = float((60.0 * scale) * croot * (1 + 0.28) / 2)
    mesh.meta = {'span': 60.0 * scale, 'length': L, 'kind': 'a330-like',
                 'skipped': mesh.meta.get('skipped', [])}
    return mesh, wakes, sref, {'span': 60.0 * scale, 'length': L}


def build_shuttle(scale=1.0, wake_kw=None):
    """Orbiter-like delta + body + single vertical tail.

    Meters at scale=1 (L~37, span~24). Returns (mesh, wakes, sref,
    info). Double-delta approximated single taper; body loft blunt.
    """
    L = 37.0 * scale
    # body stations: blunt nose -> cylinder -> boat tail (up-swept)
    prof = [(0.00, 0.03, 0.0), (0.05, 0.35, 0.0), (0.12, 0.62, 0.0),
            (0.25, 0.85, 0.0), (0.55, 1.00, 0.0), (0.75, 0.98, 0.05),
            (0.88, 0.85, 0.18), (0.96, 0.60, 0.30), (1.00, 0.30, 0.35)]
    R = 3.2 * scale
    stations = []
    for fx, rr, zc in prof:
        loop = rounded_rect(rr * R, rr * R * 0.95, 0.0, zc * R, 16)
        stations.append((fx * L, loop))
    body = loft(stations)
    # delta wing: high sweep, strong taper, mounted aft/low
    wing = _place_wing(dict(code='0012', span=24.0 * scale, chord=13.0 * scale,
                            taper=0.12, sweep_deg=55.0, n_chord=16, n_span=8,
                            cosine=False),
                       (0.52 * L, 0, -1.2 * scale))
    vtail = _vtail_up(dict(code='0012', span=7.0 * scale, chord=5.0 * scale,
                           taper=0.5, sweep_deg=45.0, n_chord=10, n_span=6,
                           cosine=False),
                      (0.80 * L, 0, 2.0 * scale))
    parts = [body, wing, vtail]
    mesh, wakes = _merge_with_wakes(
        parts, [(wing, ('y', 3.4 * scale)), (vtail, ('z', 3.5 * scale))],
        wake_kw or {})
    sref = float((24.0 * scale) * (13.0 * scale) * (1 + 0.12) / 2)
    mesh.meta = {'span': 24.0 * scale, 'length': L, 'kind': 'shuttle-like',
                 'skipped': mesh.meta.get('skipped', [])}
    return mesh, wakes, sref, {'span': 24.0 * scale, 'length': L}
