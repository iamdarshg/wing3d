"""TSD follow-up: single (geo, M) run at the refined config g120 ep16.

Follow-up to tsd_convergence.py (which found g80/g120/g160/ep16 at M0.85):
this confirms g120+ep16 at the campaign M points and other (geo, M).

Usage: python tsd_conv2.py <geo> <M>                   (g120 ep16)
       python tsd_conv2.py <geo> <M> <nx> <ny> <nz> [ep]
                                                       (custom grid/epochs;
                                                        e.g. rae g160 ep16)
Writes tests3d/tsd_conv2_<geo>_M<M>.json (custom grids get a _g<grid> tag;
default-grid runs keep the plain name).
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
    nx, ny, nz, ep = NX, NY, NZ, EP
    if len(sys.argv) >= 6:
        nx, ny, nz = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    if len(sys.argv) >= 7:
        ep = int(sys.argv[6])
    label = '' if (nx, ny, nz) == (NX, NY, NZ) else '_g%dx%dx%d' % (nx, ny, nz)
    prof = build_profile(geo)
    sec = build_section(geo)
    ref = campaign_cl(geo, M)
    print('== %s M=%.3f campaign CL=%s (g%dx%dx%d ep%d)'
          % (geo, M, ref, nx, ny, nz, ep), flush=True)
    t0 = time.time()
    base, _r, _m = build_case(0.3, 1.25, False, True, None, nx, ny, nz,
                              1.6, section=sec)
    out = 'tests3d/tsd_conv2_%s_M%.3f%s.json' % (geo, M, label)

    def checkpoint(ep_done, cl_now, nsup):
        # per-epoch checkpoint: a memory kill must not wipe the run
        json.dump({'geo': geo, 'M': M, 'grid': [nx, ny, nz], 'epochs': ep,
                   'epochs_done': ep_done + 1, 'CL': cl_now,
                   'campaign_CL': ref,
                   'delta_pct': (cl_now - ref) / ref * 100 if ref else
                   float('nan'),
                   'secs': time.time() - t0, 'complete': False},
                  open(out, 'w'), indent=1)

    L, xc, ix, j0 = run_sup(M, span=1.6, epochs=ep, free_k=False,
                            profile=prof, section=sec, base=base,
                            callback=checkpoint)
    cl = float(L['CL'])
    d = (cl - ref) / ref * 100 if ref else float('nan')
    rec = {'geo': geo, 'M': M, 'grid': [nx, ny, nz], 'epochs': ep,
           'epochs_done': ep, 'CL': cl, 'campaign_CL': ref,
           'delta_pct': d, 'secs': time.time() - t0, 'complete': True}
    json.dump(rec, open(out, 'w'), indent=1)
    print('RESULT g%dx%dx%d ep%d %s M=%.3f CL=%.5f d=%+.2f%% (%.0fs) -> %s'
          % (nx, ny, nz, ep, geo, M, cl, d, rec['secs'], out), flush=True)


if __name__ == '__main__':
    main()
