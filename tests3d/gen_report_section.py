"""Generate the CAMPAIGN_REPORT.md section-6 results block from
final_table.json (run after final_table.py). Prints markdown to stdout.

Convention: cells are CL error % (TSD vs OF); CONV OF rows plain,
still-settling OF rows (t<0.2) in parentheses, missing data '-'.
"""
import json
import os

D = r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign'
GEOS = ['naca0012', 'rae2822', 'biconvex', 'dfvlr_r4']
BAND_M = [0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975,
          1.0, 1.025, 1.05]


def main():
    rows = json.load(open(os.path.join(D, 'final_table.json')))
    print('| geo | ' + ' | '.join('%.3f' % m for m in BAND_M) + ' |')
    print('|' + '---|' * (len(BAND_M) + 1))
    for g in GEOS:
        cells = []
        for M in BAND_M:
            r = next((x for x in rows if x['geo'] == g and x['M'] == M), None)
            if r is None or r['err_pct'] is None:
                cells.append('-')
            elif r['status'] == 'CONV':
                cells.append('%+.1f' % r['err_pct'])
            else:
                cells.append('(%+.1f)' % r['err_pct'])
        print('| %s | %s |' % (g, ' | '.join(cells)))
    print()
    conv = [r for r in rows if r['err_pct'] is not None
            and r['status'] == 'CONV']
    part = [r for r in rows if r['err_pct'] is not None
            and r['status'] != 'CONV']
    miss = [r for r in rows if r['err_pct'] is None]
    n10 = sum(1 for r in conv if abs(r['err_pct']) < 10)
    print('**Coverage:** %d/%d case pairs computed (%d converged, %d '
          'settling, %d missing).' % (len(conv) + len(part), len(rows),
                                      len(conv), len(part), len(miss)))
    print()
    print('**Within 10%% (converged rows):** %d/%d overall; per geometry:'
          % (n10, len(conv)))
    for g in GEOS:
        gr = [r for r in conv if r['geo'] == g]
        if not gr:
            print('- %s: no converged pairs yet' % g)
            continue
        ge = [abs(r['err_pct']) for r in gr]
        g10 = sum(1 for e in ge if e < 10)
        print('- %s: %d/%d (mean |err| %.1f%%, max %.1f%%)'
              % (g, g10, len(gr), sum(ge) / len(ge), max(ge)))


if __name__ == '__main__':
    main()
