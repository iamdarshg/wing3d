"""TSD solve convergence study: FD grid + farfield box + epoch count.

Quantifies how much of each campaign TSD CL depends on the TSD
discretization at fixed section/panel imprint (campaign grid = g80):
  g80   :  80x32x32  (campaign grid)
  g120  : 120x48x48
  g160  : 160x64x64
  box   :  80x32x32 with the farfield box enlarged ~45%
  ep16  :  80x32x32, 16 epochs (campaign uses 8)

Usage: python tsd_convergence.py [geo ...]   (default: rae, naca, dfvlr)
Writes tests3d/tsd_convergence.json incrementally.
"""
import json
import os
import sys
import time

sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')

from tsd_campaign import build_profile, build_section
from tsd_coupled import build_case
from tsd_m_sup import run_sup

CAMPDIR = os.path.join('cases', 'openfoam', 'transonic', 'campaign')
M = 0.85

CONFIGS = [
    ('g80',  80, 32, 32, None, 8),
    ('g120', 120, 48, 48, None, 8),
    ('g160', 160, 64, 64, None, 8),
    ('box',  80, 32, 32,
     {'xlim': (-12, 16), 'ylim': (-11, 11), 'zlim': (-11, 11)}, 8),
    ('ep16', 80, 32, 32, None, 16),
]


def campaign_cl(geo):
    rows = json.load(open(os.path.join(CAMPDIR, 'tsd_%s.json' % geo)))
    for r in rows:
        if abs(r['M'] - M) < 1e-9 and r.get('CL'):
            return float(r['CL'])
    return None


def main(geos):
    out = []
    outpath = os.path.join('tests3d', 'tsd_convergence.json')
    for geo in geos:
        prof = build_profile(geo)
        sec = build_section(geo)
        ref = campaign_cl(geo)
        print('== %s M=%.2f campaign CL=%s' % (geo, M, ref), flush=True)
        for tag, nx, ny, nz, box, epochs in CONFIGS:
            t0 = time.time()
            base, _r, _m = build_case(0.3, 1.25, False, True, box, nx, ny,
                                      nz, 1.6, section=sec)
            L, xc, ix, j0 = run_sup(M, span=1.6, epochs=epochs, free_k=False,
                                    profile=prof, section=sec, base=base)
            cl = float(L['CL'])
            d = (cl - ref) / ref * 100 if ref else float('nan')
            dt = time.time() - t0
            print('   %-5s %3dx%dx%d ep=%d box=%d CL=%.5f  d=%+.2f%%  (%.0fs)'
                  % (tag, nx, ny, nz, epochs, bool(box), cl, d, dt),
                  flush=True)
            out.append({'geo': geo, 'M': M, 'tag': tag, 'grid': [nx, ny, nz],
                        'epochs': epochs, 'box': box, 'CL': cl,
                        'campaign_CL': ref, 'delta_pct': d, 'secs': dt})
            json.dump(out, open(outpath, 'w'), indent=1)
    print('DONE ->', outpath, flush=True)


if __name__ == '__main__':
    main(sys.argv[1:] or ['rae2822', 'naca0012', 'dfvlr_r4'])
