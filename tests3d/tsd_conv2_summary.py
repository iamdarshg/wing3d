"""Summarize the g120ep16 confirmation runs (tests3d/tsd_conv2_*.json)
against the campaign values and the study's same-grid short-epoch values.

Usage: python tests3d/tsd_conv2_summary.py
"""
import glob
import json
import os

STUDY = os.path.join('tests3d', 'tsd_convergence.json')


def main():
    study = json.load(open(STUDY))
    # same-grid study value per (geo, M, grid)
    key = {(r['geo'], round(r['M'], 6), tuple(r['grid'])): r for r in study}
    files = sorted(glob.glob(os.path.join('tests3d', 'tsd_conv2_*.json')))
    if not files:
        print('no conv2 results yet')
        return
    print('%-9s %6s %10s %11s %11s %11s %9s %9s'
          % ('geo', 'M', 'grid', 'campaign', 'study ep8', 'conv2', 'vs camp',
             'vs ep8'))
    for f in files:
        r = json.load(open(f))
        grid = tuple(r['grid'])
        st = key.get((r['geo'], round(r['M'], 6), grid))
        ep8 = st['CL'] if st else float('nan')
        ref = r['campaign_CL']
        d1 = (r['CL'] - ref) / ref * 100
        d2 = (r['CL'] - ep8) / ep8 * 100 if st else float('nan')
        gs = 'x'.join(str(v) for v in grid)
        print('%-9s %6.3f %10s %11.6f %11.6f %11.6f %+8.2f%% %+8.2f%%'
              % (r['geo'], r['M'], gs, ref, ep8, r['CL'], d1, d2))


if __name__ == '__main__':
    main()
