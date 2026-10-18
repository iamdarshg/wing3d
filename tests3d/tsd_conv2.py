"""TSD follow-up: single (geo, M) run at the refined config g120 ep16.

Follow-up to tsd_convergence.py (which found g80/g120/g160/ep16 at M0.85):
this confirms g120+ep16 at the campaign M points and other (geo, M).

Usage: python tsd_conv2.py <geo> <M>      (one process per run)
Writes tests3d/tsd_conv2_<geo>_M<M>.json.
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
NX, NY, NZ, EP = 120, 48, 48, 16


def campaign_cl(geo, M):
    rows = json.load(open(os.path.join(CAMPDIR, 'tsd_%s.json' % geo)))
    for r in rows:
        if abs(r['M'] - M) < 1e-9 and r.get('CL'):
            return float(r['CL'])
    return None


def main():
    geo = sys.argv[1]
    M = float(sys.argv[2])
    prof = build_profile(geo)
    sec = build_section(geo)
    ref = campaign_cl(geo, M)
    print('== %s M=%.3f campaign CL=%s (g120 ep16)'
          % (geo, M, ref), flush=True)
    t0 = time.time()
    base, _r, _m = build_case(0.3, 1.25, False, True, None, NX, NY, NZ,
                              1.6, section=sec)
    L, xc, ix, j0 = run_sup(M, span=1.6, epochs=EP, free_k=False,
                            profile=prof, section=sec, base=base)
    cl = float(L['CL'])
    d = (cl - ref) / ref * 100 if ref else float('nan')
    rec = {'geo': geo, 'M': M, 'grid': [NX, NY, NZ], 'epochs': EP,
           'CL': cl, 'campaign_CL': ref, 'delta_pct': d,
           'secs': time.time() - t0}
    out = 'tests3d/tsd_conv2_%s_M%.3f.json' % (geo, M)
    json.dump(rec, open(out, 'w'), indent=1)
    print('RESULT g120ep16 %s M=%.3f CL=%.5f d=%+.2f%% (%.0fs) -> %s'
          % (geo, M, cl, d, rec['secs'], out), flush=True)


if __name__ == '__main__':
    main()
