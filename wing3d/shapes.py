"""Complex 3D test shapes: lofted bodies, Ahmed-like car, F-16-like config.

All single watertight surfaces (lofted) except wing/fuselage junctions on
the F-16, which use overlapping closed solids (documented approximation).
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
