"""TSD full-potential campaign, matched to the OF Euler points.

For each (geometry-like section, M) in the M0.85-1.1 band, run the
supersonic-branch TSD (inflow-pinned/outflow-free) and record CL/CD +
error vs the OF anchors. NOTE: TSD section support is NACA0012-native
(thick=0.12); RAE2822/biconvex/DFVLR sections are matched by thickness
proxy (t/c) - honest caveat. M<0.95 uses free_k, M>=0.95 fixed_k free
as validated. Records into campaign results json.

Usage: python tsd_campaign.py [geo]
  geo: naca0012 | rae2822 | biconvex | dfvlr_r4
"""
import numpy as np
import json
import os
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')

THICK = {'naca0012': 0.12, 'rae2822': 0.121, 'biconvex': 0.08,
         'dfvlr_r4': 0.093}
MACHES = [0.85, 0.875, 0.90, 0.925, 0.95, 0.975, 1.00, 1.025, 1.05, 1.075, 1.10]


def build_profile(geo):
    """Section profile (slope+half-thickness) for TSD, from coordinate fns."""
    from wing3d.tsd import section_profile_xy
    from gen_transonic_stl import (naca0012_xy, rae2822_xy,
                                   biconvex_xy, dfvlr_r4_xy)
    if geo == 'naca0012':
        xa, yt = naca0012_xy(n=400)
        xu, yu, xl, yl = xa, yt, xa, -yt
    elif geo == 'rae2822':
        xu, yu, xl, yl = rae2822_xy()
    elif geo == 'biconvex':
        xu, yu = biconvex_xy(n=400)
        xl, yl = xu.copy(), -yu
    elif geo == 'dfvlr_r4':
        xu, yu, xl, yl = dfvlr_r4_xy()
    else:
        raise ValueError(geo)
    xg = np.linspace(0, 1, 400)
    # xl may be descending (TE->LE): np.interp requires ascending xp
    if xl[0] > xl[-1]:
        xl = xl[::-1]; yl = yl[::-1]
    if xu[0] > xu[-1]:
        xu = xu[::-1]; yu = yu[::-1]
    return section_profile_xy(xg, np.interp(xg, xu, yu),
                              np.interp(xg, xl, yl))


def run_tsd(M, geo, span=1.6, epochs=8):
    from tsd_m_sup import run_sup
    free_k = M < 0.95
    prof = build_profile(geo)
    s = run_sup(M, span=span, epochs=epochs, free_k=free_k, profile=prof)
    return s


def main(geo='naca0012'):
    from tsd_m_sup import run_sup
    prof = build_profile(geo)
    out = []
    for M in MACHES:
        free_k = M < 0.95
        try:
            L, xc, ix, j0 = run_sup(M, span=1.6, epochs=8, free_k=free_k,
                                    profile=prof)
            # full-field CL + centerline section cl
            jy = np.arange(0)  # placeholder
            out.append({'M': M, 'CL': float(L['CL']), 'free_k': free_k,
                        'geo': geo})
            print('M=%.3f CL=%.4f free_k=%s' % (M, L['CL'], free_k),
                  flush=True)
        except Exception as e:
            print('M=%.3f FAIL %s' % (M, e), flush=True)
            out.append({'M': M, 'CL': None, 'err': str(e), 'geo': geo})
    res_path = os.path.join(
        r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign',
        'tsd_%s.json' % geo)
    with open(res_path, 'w') as f:
        json.dump(out, f, indent=1)
    print('saved', res_path)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'naca0012')