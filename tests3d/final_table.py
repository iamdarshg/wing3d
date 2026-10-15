"""Final campaign table: OF Euler vs TSD, band M0.75-1.05, 4 geometries.

OF values: old corpus (container_coeffs) overlaid by new runs
(results2 mirror synced from the container). Column layout of the raw
coefficient.dat rows: 1=Time 2=Cd 5=Cl (7=Cl(r) -- NOT Cd).
Validity rules:
  - All geometries: old data valid (mesh forensics: the rae mesh cache was
    byte-identical pre/post-remesh -> the original corpus always solved the
    correct geometry); new data wins where present.
  - dfvlr_r4 old 0.9/1.075/1.100 are partial runs -> treated as invalid.
Output: printed table + final_table.json + final_table.md.
"""
import os
import json
import glob

D = r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign'
COEFF = os.path.join(D, 'container_coeffs')
RES2 = os.path.join(D, 'results2')
GEOS = ['naca0012', 'rae2822', 'biconvex', 'dfvlr_r4']
M_LO, M_HI = 0.75, 1.05
BAND_M = [0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95,
          0.975, 1.0, 1.025, 1.05]
NEW_ONLY = set()
INVALID_OLD = {('dfvlr_r4', 0.9), ('dfvlr_r4', 1.075), ('dfvlr_r4', 1.1)}


def parse_row(line):
    v = line.split()
    return float(v[0]), float(v[4]), float(v[1])  # t, Cl, Cd


def old_rows():
    out = {}
    for f in glob.glob(os.path.join(COEFF, '*.coeff')):
        name = os.path.basename(f)[:-6]
        geo = name.split('_M')[0]
        M = float(name.split('_M')[1].split('_')[0])
        rows = [l for l in open(f).read().splitlines()
                if l.strip() and not l.strip().startswith('#')]
        if not rows:
            continue
        out.setdefault(geo, {})[M] = parse_row(rows[-1]) + ('old',)
    return out


def new_rows():
    out = {}
    for f in glob.glob(os.path.join(RES2, '*.dat')):
        name = os.path.basename(f)[:-4]
        if '_M' not in name:
            continue
        geo = name.split('_M')[0]
        M = float(name.split('_M')[1].split('_')[0])
        txt = open(f).read().strip()
        if not txt:
            continue
        out.setdefault(geo, {})[M] = parse_row(txt.splitlines()[-1]) + ('new',)
    return out


def tsd_rows():
    out = {}
    for g in GEOS:
        p = os.path.join(D, 'tsd_%s.json' % g)
        if not os.path.exists(p):
            continue
        for r in json.load(open(p)):
            if r.get('CL') is not None:
                out.setdefault(g, {})[round(r['M'], 3)] = r['CL']
    return out


def main():
    old, new, tsd = old_rows(), new_rows(), tsd_rows()
    rows = []
    print('%-9s %-6s %-7s %-9s %-9s %-9s %-8s %s' % (
        'geo', 'M', 'status', 'OF_CL', 'OF_CD', 'TSD_CL', 'err%', 'src'))
    for g in GEOS:
        for M in BAND_M:
            o = old.get(g, {}).get(M)
            if o is not None and (g, M) in INVALID_OLD:
                o = None
            n = new.get(g, {}).get(M)
            if g in NEW_ONLY:
                src = n
            else:
                src = n if n is not None else o
            tv = tsd.get(g, {}).get(round(M, 3))
            st, cl, cd, tag = '-', None, None, '-'
            if src is not None:
                t, cl, cd, tag = src
                st = 'CONV' if t >= 0.2 else 'part'
            err = ((tv - cl) / cl * 100) if (tv is not None and cl) else None
            rows.append({'geo': g, 'M': M, 'status': st, 'of_t': src[0]
                         if src else None, 'of_cl': cl, 'of_cd': cd,
                         'tsd_cl': tv, 'err_pct': err, 'source': tag})
            print('%-9s %-6.3f %-7s %-9s %-9s %-9s %-8s %s' % (
                g, M, st,
                '%.4f' % cl if cl is not None else '-',
                '%.4f' % cd if cd is not None else '-',
                '%.4f' % tv if tv is not None else '-',
                '%+.1f' % err if err is not None else '-', tag))
    json.dump(rows, open(os.path.join(D, 'final_table.json'), 'w'), indent=1)
    # markdown
    lines = ['| geo | M | OF CL | TSD CL | err % | status |',
             '|---|---|---|---|---|---|']
    conv = [r for r in rows if r['err_pct'] is not None
            and r['status'] == 'CONV']
    for r in rows:
        lines.append('| %s | %.3f | %s | %s | %s | %s |' % (
            r['geo'], r['M'],
            '%.4f' % r['of_cl'] if r['of_cl'] is not None else '-',
            '%.4f' % r['tsd_cl'] if r['tsd_cl'] is not None else '-',
            '%+.1f' % r['err_pct'] if r['err_pct'] is not None else '-',
            r['status']))
    open(os.path.join(D, 'final_table.md'), 'w').write('\n'.join(lines))
    # condensed error matrix (geo x M): the report's central deliverable view
    mtx = ['| geo | ' + ' | '.join('%.3f' % m for m in BAND_M) + ' |',
           '|' + '---|' * (len(BAND_M) + 1)]
    for g in GEOS:
        cells = []
        for M in BAND_M:
            r = next((x for x in rows if x['geo'] == g and x['M'] == M),
                     None)
            cells.append('%+.1f' % r['err_pct']
                         if (r and r['err_pct'] is not None) else '-')
        mtx.append('| %s | %s |' % (g, ' | '.join(cells)))
    open(os.path.join(D, 'final_table_matrix.md'), 'w').write('\n'.join(mtx))
    print()
    print('\n'.join(mtx))
    if conv:
        errs = [abs(r['err_pct']) for r in conv]
        n10 = sum(1 for e in errs if e < 10)
        print('CONV rows: %d/%d under 10%% (max |err| %.1f%%, mean %.1f%%)'
              % (n10, len(conv), max(errs), sum(errs) / len(errs)))
        for g in GEOS:
            gr = [r for r in conv if r['geo'] == g]
            if gr:
                ge = [abs(r['err_pct']) for r in gr]
                g10 = sum(1 for e in ge if e < 10)
                print('  %-9s: %d/%d under 10%%  (mean %.1f%%, max %.1f%%)'
                      % (g, g10, len(gr), sum(ge) / len(ge), max(ge)))
    print('saved final_table.json/md (%d rows)' % len(rows))


if __name__ == '__main__':
    main()
