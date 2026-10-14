"""M0.85-1.1 band analysis: OF campaign vs TSD vs legacy anchors.

Outputs a per-(geo,M) table with OF CL/CD, TSD CL, and % errors, plus
the interpolated band model for the UI.
"""
import os
import json
import glob
import numpy as np

CAMP = r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign'
OUT = os.path.join(CAMP, 'band_results.json')

# legacy OF anchors (AR1.6 slab, NACA0012, a=1.25) for continuity check
LEGACY = {0.8: (0.0686, 0.0485), 0.885: (0.0696, 0.0512),
          0.95: (0.0830, 0.0650), 1.0: (0.0790, 0.0770),
          1.1: (0.0533, 0.1030)}


def parse_case(name):
    """Read Cl/Cd from the log.run final converged values."""
    log = os.path.join(CAMP, name, 'log.run')
    if not os.path.exists(log):
        return None
    cl = cd = None
    try:
        with open(log, 'r', errors='ignore') as f:
            for line in f:
                s = line.strip()
                if s.startswith('Cl:'):
                    cl = float(s.split()[1])
                elif s.startswith('Cd:'):
                    cd = float(s.split()[1])
    except Exception:
        return None
    return (cl, cd)


def collect_of():
    rows = {}
    for d in sorted(glob.glob(os.path.join(CAMP, '*_M*.done'))):
        name = os.path.basename(d)[:-5]
        r = parse_case(name)
        if r is None:
            continue
        geo = name.split('_M')[0]
        M = float(name.split('_M')[1].split('_')[0])
        rows.setdefault(geo, {})[M] = r
    return rows


def collect_tsd():
    rows = {}
    for p in glob.glob(os.path.join(CAMP, 'tsd_*.json')):
        geo = os.path.basename(p)[4:-5]
        try:
            data = json.load(open(p))
        except Exception:
            continue
        rows[geo] = {r['M']: r.get('CL') for r in data}
    return rows


def main():
    of = collect_of()
    tsd = collect_tsd()
    out = {'legacy_naca0012': LEGACY, 'rows': []}
    for geo in sorted(set(list(of) + list(tsd))):
        for M in sorted(set(list(of.get(geo, {})) + list(tsd.get(geo, {})))):
            ofv = of.get(geo, {}).get(M)
            tsv = tsd.get(geo, {}).get(M)
            row = {'geo': geo, 'M': M, 'of_cl': ofv[0] if ofv else None,
                   'of_cd': ofv[1] if ofv else None,
                   'tsd_cl': tsv if tsv else None}
            if ofv and tsv and abs(ofv[0]) > 1e-9:
                row['err_pct'] = (tsv - ofv[0]) / ofv[0] * 100
            out['rows'].append(row)
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=1)
    # console table
    print('%-10s %-6s %10s %10s %10s %10s' % (
        'geo', 'M', 'OF_CL', 'OF_CD', 'TSD_CL', 'err%'))
    for r in sorted(out['rows'], key=lambda x: (x['geo'], x['M'])):
        if r['of_cl'] is None:
            continue
        print('%-10s %-6.3f %10.4f %10.4f %10s %10s' % (
            r['geo'], r['M'], r['of_cl'], r['of_cd'],
            '%.4f' % r['tsd_cl'] if r['tsd_cl'] is not None else '-',
            '%+.1f' % r['err_pct'] if r.get('err_pct') is not None else '-'))
    print('saved', OUT)


if __name__ == '__main__':
    main()