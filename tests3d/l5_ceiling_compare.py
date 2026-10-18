"""L5 ceiling probe comparison: dfvlr M0.925 / M1.000 re-run on the 86k
mesh vs the corpus 45.6k values (direct mesh-dissipation test).

Usage: python tests3d/l5_ceiling_compare.py
Expects the pulled probe finals:
  tests3d/l5_probe_dfvlr_M0.925_86k.dat
  tests3d/l5_probe_dfvlr_M1.000_86k.dat
(12-column coefficient.dat lines; CL is column 5, CD column 7.)
"""
import json
import os

D = os.path.join('cases', 'openfoam', 'transonic', 'campaign')
PROBES = {
    0.925: os.path.join('tests3d', 'l5_probe_dfvlr_M0.925_86k.dat'),
    1.000: os.path.join('tests3d', 'l5_probe_dfvlr_M1.000_86k.dat'),
}


def main():
    rows = json.load(open(os.path.join(D, 'final_table.json')))
    print('%-7s %10s %10s %8s | %10s %10s %8s | %s'
          % ('M', 'CL 45.6k', 'CL 86k', 'dCL', 'CD 45.6k', 'CD 86k',
             'dCD', 'verdict'))
    for M, path in sorted(PROBES.items()):
        old = next(r for r in rows if r['geo'] == 'dfvlr_r4'
                   and abs(r['M'] - M) < 1e-9)
        cols = open(path).read().split()
        cl_new, cd_new = float(cols[4]), float(cols[6])
        dcl = (cl_new - old['of_cl']) / old['of_cl'] * 100
        dcd = (cd_new - old['of_cd']) / old['of_cd'] * 100
        if abs(dcl) >= 3:
            v = 'resolution matters (|dCL| >= 3%)'
        elif abs(dcl) >= 1:
            v = 'marginally resolved (1-3%)'
        else:
            v = 'flat (<1%) - OF saturation genuine'
        print('%-7.3f %10.6f %10.6f %+7.2f%% | %10.6f %10.6f %+7.2f%% | %s'
              % (M, old['of_cl'], cl_new, dcl, old['of_cd'], cd_new, dcd, v))


if __name__ == '__main__':
    main()
