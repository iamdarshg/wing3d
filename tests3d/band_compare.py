"""M0.85-1.1 band: OF Euler vs TSD fullpot comparison + error table.

Reads the campaign coefficient data (host copy) + tsd_*.json, produces
per-(geo,M) rows with % error vs OF and vs TSD. This is the flagship
comparison the user asked for.
"""
import os
import json
import glob
import numpy as np

D = r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign'
COEFF = os.path.join(D, 'container_coeffs')
GEOS = ['naca0012', 'rae2822', 'biconvex', 'dfvlr_r4']


def of_rows():
    out = {}
    for f in glob.glob(os.path.join(COEFF, '*.coeff')):
        name = os.path.basename(f)[:-6]
        geo = name.split('_M')[0]
        M = float(name.split('_M')[1].split('_')[0])
        rows = [l.split() for l in open(f).read().splitlines()
                if l.strip() and not l.strip().startswith('#')]
        if not rows:
            continue
        # cols: 1=Time 2=Cd 5=Cl (7 is Cl(r), NOT Cd)
        t, cl, cd = float(rows[-1][0]), float(rows[-1][4]), float(rows[-1][1])
        out.setdefault(geo, {})[M] = (t, cl, cd)
    return out


def tsd_rows():
    out = {}
    for g in GEOS:
        p = os.path.join(D, 'tsd_%s.json' % g)
        if not os.path.exists(p):
            continue
        for r in json.load(open(p)):
            out.setdefault(g, {})[r['M']] = r['CL']
    return out


def main():
    of = of_rows()
    tsd = tsd_rows()
    print('%-10s %-6s %-8s %-10s %-10s %-10s' % (
        'geo', 'M', 'status', 'OF_CL', 'TSD_CL', 'err%'))
    rows_out = []
    for g in GEOS:
        for M in sorted(of.get(g, {})):
            t, cl, cd = of[g][M]
            tv = tsd.get(g, {}).get(M)
            err = ((tv - cl) / cl * 100) if (tv is not None and abs(cl) > 1e-9) else None
            st = 'CONV' if t > 0.2 else 'part'
            print('%-10s %-6.3f %-8s %-10.5f %-10s %-10s' % (
                g, M, st, cl,
                '%.5f' % tv if tv is not None else '-',
                '%+.1f' % err if err is not None else '-'))
            rows_out.append({'geo': g, 'M': M, 'of_cl': cl, 'of_cd': cd,
                             'tsd_cl': tv, 'err_pct': err, 't': t})
    with open(os.path.join(D, 'band_compare.json'), 'w') as f:
        json.dump(rows_out, f, indent=1)
    print('saved band_compare.json (%d rows)' % len(rows_out))


if __name__ == '__main__':
    main()