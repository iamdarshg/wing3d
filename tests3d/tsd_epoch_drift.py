"""Parse tests3d/camp_*.txt epoch logs: per-M TSD convergence diagnostics.

For each (geo, M) block: final CL, mean/sd of last 4 epochs, drift
|last - mean(last4)|, and the worst epoch jump in the last 4.
"""
import glob
import os
import re

for path in sorted(glob.glob(os.path.join('tests3d', 'camp_*.txt'))):
    geo = os.path.basename(path)[5:-4]
    eps = []
    print('== %s ==' % geo)
    print('%-7s %9s %9s %8s %8s %8s' % ('M', 'final', 'mean4', 'sd4', 'drift%', 'maxjump%'))
    rows = []
    for ln in open(path, encoding='utf-16', errors='replace'):
        m = re.match(r'\s*ep\s+(\d+):\s*CL=([0-9.]+)', ln)
        if m:
            eps.append(float(m.group(2)))
            continue
        m = re.match(r'M=([0-9.]+)\s+CL=', ln)
        if m and eps:
            M = float(m.group(1))
            fin = eps[-1]
            last4 = eps[-4:]
            mu = sum(last4) / len(last4)
            sd = (sum((x - mu) ** 2 for x in last4) / len(last4)) ** 0.5
            drift = abs(fin - mu) / fin * 100
            jmp = max(abs(last4[i] - last4[i - 1]) for i in range(1, len(last4))) / fin * 100
            rows.append((M, fin, mu, sd, drift, jmp))
            print('%-7.3f %9.4f %9.4f %8.4f %7.2f%% %7.2f%%' % (M, fin, mu, sd, drift, jmp))
            eps = []
    if rows:
        d_mean = sum(r[4] for r in rows) / len(rows)
        j_mean = sum(r[5] for r in rows) / len(rows)
        print('  --> mean drift %.3f%%  mean last-jump %.3f%%  max drift %.2f%%'
              % (d_mean, j_mean, max(r[4] for r in rows)))
    print()
