"""CLI validation battery: run the shipped STL CLI (wing3d.tsdcli) on the
delivered 923k-stored-vertex STLs across the M0.75-1.05 band and report
full-wing CL error % vs the OpenFOAM corpus.

Usage:
    python tests3d/cli_battery.py [--jobs 3]

Outputs:
    tests3d/cli_battery_runs/<geo>_M<mach>.json   (per-run raw CLI JSON)
    tests3d/cli_battery_runs/<geo>_M<mach>.log    (per-run stdout)
    tests3d/cli_battery.json                      (assembled battery table)

The battery is the end-to-end check of the shipped CLI path (STL ->
mid-span section -> TSD full-potential) against the same OpenFOAM anchors
used in the campaign report (campaign/final_table.json).
"""
import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STL_DIR = os.path.join(ROOT, 'cases', 'openfoam', 'transonic', 'stl')
CAMPAIGN = os.path.join(ROOT, 'cases', 'openfoam', 'transonic', 'campaign')
OUTDIR = os.path.join(ROOT, 'tests3d', 'cli_battery_runs')

# a dozen pairs spanning all four geometries and the full error spectrum
PAIRS = [
    ('rae2822', 0.750), ('rae2822', 0.850), ('rae2822', 0.925),
    ('rae2822', 0.975), ('rae2822', 1.000),
    ('naca0012', 0.750), ('naca0012', 1.050),
    ('biconvex', 0.850), ('biconvex', 1.050),
    ('dfvlr_r4', 0.850), ('dfvlr_r4', 0.950), ('dfvlr_r4', 1.050),
]


def run_pair(geo, M):
    name = '%s_M%.3f' % (geo, M)
    out = os.path.join(OUTDIR, name + '.json')
    log = os.path.join(OUTDIR, name + '.log')
    stl = os.path.join(STL_DIR, geo + '.stl')
    cmd = [sys.executable, '-m', 'wing3d.tsdcli', stl,
           '--mach', '%.3f' % M, '--epochs', '8', '--json', out]
    t0 = time.time()
    with open(log, 'w') as f:
        r = subprocess.run(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
    return name, r.returncode, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=3)
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)

    of_rows = json.load(open(os.path.join(CAMPAIGN, 'final_table.json')))
    of = {(r['geo'], round(r['M'], 3)): r for r in of_rows}

    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_pair, g, m): (g, m) for g, m in PAIRS}
        n = 0
        for f in as_completed(futs):
            name, rc, dt = f.result()
            n += 1
            print('[%d/%d] %-18s rc=%d  %.0f s' % (n, len(PAIRS), name, rc, dt))

    rows = []
    for g, m in PAIRS:
        name = '%s_M%.3f' % (g, m)
        out = os.path.join(OUTDIR, name + '.json')
        rec = {'geo': g, 'M': m}
        if os.path.isfile(out):
            d = json.load(open(out))
            rec.update(cli_cl=d['CL'], cli_cd=d['CD'],
                       npanels=d['npanels'], nverts=d['nverts'],
                       thick_max=d['thick_max'], camber_max=d['camber_max'],
                       mcrit=d['mcrit'])
            r = of.get((g, round(m, 3)))
            if r is not None and r.get('of_cl'):
                rec['of_cl'] = r['of_cl']
                rec['err_pct'] = (d['CL'] - r['of_cl']) / r['of_cl'] * 100.0
            if r is not None and r.get('tsd_cl'):
                rec['campaign_tsd_cl'] = r['tsd_cl']
                rec['delta_vs_campaign_pct'] = \
                    (d['CL'] - r['tsd_cl']) / r['tsd_cl'] * 100.0
        rows.append(rec)

    res = {'generated': time.strftime('%Y-%m-%d %H:%M'),
           'note': 'wing3d.tsdcli (STL-driven, epochs=8) vs OpenFOAM corpus',
           'rows': rows}
    with open(os.path.join(ROOT, 'tests3d', 'cli_battery.json'), 'w') as f:
        json.dump(res, f, indent=1)

    print()
    print('| geo | M | CLI CL | OF CL | err % | vs campaign % |')
    print('|---|---|---|---|---|---|')
    for r in rows:
        if 'cli_cl' not in r:
            print('| %s | %.3f | - | - | - | - |' % (r['geo'], r['M']))
            continue
        print('| %s | %.3f | %.4f | %s | %s | %s |' % (
            r['geo'], r['M'], r['cli_cl'],
            '%.4f' % r['of_cl'] if 'of_cl' in r else '-',
            ('%+.1f' % r['err_pct']) if 'err_pct' in r else '-',
            ('%+.2f' % r['delta_vs_campaign_pct'])
            if 'delta_vs_campaign_pct' in r else '-'))


if __name__ == '__main__':
    main()
