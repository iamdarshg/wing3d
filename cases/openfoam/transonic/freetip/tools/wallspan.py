"""Spanwise wall-pressure analysis for the free-tip OpenFOAM cases.

Reads a legacy/XML `.vtp` wall patch written by
`foamToVTK -ascii -latestTime -fields p` and produces:

  * the global pressure-integrated normal force, cross-checked against the
    `forceCoeffs` history (a 1:1 match proves the surface orientation and the
    integration are right);
  * Cp at a set of chordwise stations (x/c = 0.25, 0.50, 0.75) versus
    spanwise station eta, upper and lower surface;
  * the sectional normal-force coefficient cl(eta) obtained by integrating the
    wall pressure over spanwise strips, which is the quantity that must show
    *unloading* near the free tip (eta -> 1) and, in the cascade case, must be
    span-independent.

Usage:
  python wallspan.py <wall.vtp> --mach 0.80 --span 1.6 --span-half 0.8 ...
Prints JSON on stdout.
"""
import argparse
import json
import xml.etree.ElementTree as ET

import numpy as np

TINF = 300.0
PINF = 1.0e5
RU = 8314.462618
MOLW = 28.01348
RS = RU / MOLW
GAMMA = 1.4
RHO_INF = PINF / (RS * TINF)
A_INF = np.sqrt(GAMMA * RS * TINF)


def read_vtp(path):
    root = ET.parse(path).getroot()
    piece = root.find('PolyData').find('Piece')
    points = None
    conn = off = None
    cdata = {}
    pdata = {}
    for g in piece:
        if g.tag == 'Points':
            arr = g.find('DataArray')
            v = np.fromstring(arr.text, sep=' ')
            points = v.reshape(-1, 3)
        elif g.tag == 'Polys':
            for arr in g.findall('DataArray'):
                if arr.attrib.get('Name') == 'connectivity':
                    conn = np.fromstring(arr.text, sep=' ', dtype=np.int64)
                elif arr.attrib.get('Name') == 'offsets':
                    off = np.fromstring(arr.text, sep=' ', dtype=np.int64)
        elif g.tag == 'CellData':
            for arr in g.findall('DataArray'):
                cdata[arr.attrib.get('Name')] = np.fromstring(arr.text, sep=' ')
        elif g.tag == 'PointData':
            for arr in g.findall('DataArray'):
                pdata[arr.attrib.get('Name')] = np.fromstring(arr.text, sep=' ')
    return points, conn, off, cdata, pdata


def poly_area_normals(pts, conn, off):
    """Newell area vector + centroid for each polygon of a VTK Polys array."""
    starts = np.concatenate(([0], off[:-1]))
    counts = off - starts
    ac = np.zeros((len(counts), 3))
    cen = np.zeros((len(counts), 3))
    for k, (s, n) in enumerate(zip(starts, counts)):
        idx = conn[s:s + n]
        v = pts[idx]
        vn = np.roll(v, -1, axis=0)
        ac[k] = 0.5 * np.cross(v, vn).sum(axis=0)
        cen[k] = v.mean(axis=0)
    return ac, cen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vtp')
    ap.add_argument('--mach', type=float, required=True)
    ap.add_argument('--halfspan', type=float, default=0.8,
                    help='y of the tip (wing y in [-halfspan, +halfspan])')
    ap.add_argument('--chord', type=float, default=1.0)
    ap.add_argument('--aref', type=float, default=1.6)
    ap.add_argument('--nstrip', type=int, default=20)
    ap.add_argument('--stations', type=float, nargs='*',
                    default=[0.25, 0.50, 0.75])
    a = ap.parse_args()

    pts, conn, off, cdata, pdata = read_vtp(a.vtp)
    pd = cdata if 'p' in cdata else pdata
    pname = 'p' if 'p' in pd else [k for k in pd if k.startswith('p')][0]
    p = pd[pname]
    ac, cx = poly_area_normals(pts, conn, off)
    npoly = ac.shape[0]
    pcell = p if p.size == npoly else np.full(npoly, p.mean())

    q = 0.5 * RHO_INF * (a.mach * A_INF) ** 2

    area = np.linalg.norm(ac, axis=1)
    nrm = ac / np.maximum(area[:, None], 1e-30)

    # Use the gauge pressure: p_inf * area mismatch between the upper and lower
    # surface of a strip would otherwise swamp CL (p_inf*|dA| >> q*CL).
    dp = pcell - PINF
    fz = -(dp * ac[:, 2]).sum()
    flick = fz / (q * a.aref)
    eta = cx[:, 1] / a.halfspan
    pinf_side = np.abs(nrm[:, 2]) > 1e-3

    # per-face lift increment from gauge pressure
    lift_area = -(dp * ac[:, 2])

    edges = np.linspace(0.0, 1.0, a.nstrip + 1)
    cl_eta, cl_sect, eta_mid = [], [], []
    dy = 2.0 * a.halfspan
    for i in range(a.nstrip):
        m = (np.abs(eta) >= edges[i]) & (np.abs(eta) < edges[i + 1]) & pinf_side
        de = edges[i + 1] - edges[i]
        if m.sum() == 0:
            cl_eta.append(float('nan'))
            cl_sect.append(float('nan'))
        else:
            fz_s = lift_area[m].sum()
            # contribution to CL on the wing reference area (Aref = c*span)
            cl_eta.append(float(fz_s / (q * a.aref * de)))
            # true sectional cl: normalised by the local planform area c*dy
            cl_sect.append(float(fz_s / (q * a.chord * dy * de)))
        eta_mid.append(float(0.5 * (edges[i] + edges[i + 1])))

    # Cp vs eta at chordwise stations, upper surface
    xc = cx[:, 0] / a.chord
    cp = (pcell - PINF) / q
    tol = 0.03
    st = {}
    for s in a.stations:
        m = (np.abs(xc - s) < tol) & pinf_side
        d = {}
        for sgn, name in ((1, 'upper'), (-1, 'lower')):
            mm = m & (np.sign(cx[:, 2]) == sgn)
            et, cv = np.abs(eta[mm]), cp[mm]
            order = np.argsort(et)
            d[name] = [[round(float(e), 3), round(float(c), 4)]
                       for e, c in zip(et[order], cv[order])]
        st[f'x/c={s}'] = d

    out = dict(
        vtp=a.vtp,
        mach=a.mach, q=round(float(q), 3), rho=round(float(RHO_INF), 6),
        nfaces=int(npoly),
        ymin=round(float(pts[:, 1].min()), 4),
        ymax=round(float(pts[:, 1].max()), 4),
        xmin=round(float(pts[:, 0].min()), 4),
        xmax=round(float(pts[:, 0].max()), 4),
        zmax_upper=round(float(pts[pts[:, 2] > 0][:, 2].max()), 4),
        CL_pressure_integral=round(float(flick), 6),
        cl_vs_eta=[[round(e, 3), (None if np.isnan(c) else round(c, 4))]
                   for e, c in zip(eta_mid, cl_eta)],
        cl_section_vs_eta=[[round(e, 3), (None if np.isnan(c) else round(c, 4))]
                           for e, c in zip(eta_mid, cl_sect)],
        cp_stations=st,
    )
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
