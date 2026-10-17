"""Scope: how much of each cell's error is explained by the Mach-label
convention (OF thermodynamic Mach = 1.01683 x label)?

For each CONV cell: interp TSD CL at M*1.01683 (the Mach the OF row
actually ran at), compare vs the same OF CL. Prints delta vs the current
error column.
"""
import json
import os

R = 1.0 / (353.05867 / 347.2235)  # a_label / a_thermo -> thermo factor
# a_label  = sqrt(gamma p/rho) with p=1e5, rho=1.1231  = 353.0587
# a_thermo = sqrt(gamma R_air T) with R=287.05, T=300 = 347.2235
F = 353.05867 / 347.2235

D = os.path.join('cases', 'openfoam', 'transonic', 'campaign', 'final_table.json')
rows = json.load(open(D))
geos = ['naca0012', 'rae2822', 'biconvex', 'dfvlr_r4']

def interp(xs, ys, x):
    if x <= xs[0]:
        i0, i1 = 0, 1
    elif x >= xs[-1]:
        i0, i1 = len(xs) - 2, len(xs) - 1
    else:
        i0 = max(i for i in range(len(xs)) if xs[i] <= x)
        i1 = i0 + 1
    t = (x - xs[i0]) / (xs[i1] - xs[i0])
    return ys[i0] + t * (ys[i1] - ys[i0])

print('%-9s %6s %9s %9s %9s %8s' % ('geo', 'M', 'err_now%', 'err_therm%', 'delta_pp', ''))
tot = []
for g in geos:
    gr = [r for r in rows if r['geo'] == g]
    gm = sorted(gr, key=lambda r: r['M'])
    xs = [r['M'] for r in gm]
    ys = [r['tsd_cl'] for r in gm]
    for r in gm:
        if r['status'] != 'CONV' or r['err_pct'] is None:
            continue
        M = r['M']
        cl_t = interp(xs, ys, M * F)
        err_t = (cl_t - r['of_cl']) / r['of_cl'] * 100.0
        d = err_t - r['err_pct']
        tot.append((g, M, d))
        print('%-9s %6.3f %+9.1f %+9.1f %+9.2f' % (g, M, r['err_pct'], err_t, d))

print()
import statistics
for g in geos:
    ds = [t[2] for t in tot if t[0] == g]
    if ds:
        print('%-9s mean |delta| %.2f pp   max %.2f pp   (n=%d)'
              % (g, sum(abs(x) for x in ds) / len(ds), max(abs(x) for x in ds), len(ds)))
ds = [t[2] for t in tot]
print('ALL       mean |delta| %.2f pp   max %.2f pp' % (sum(abs(x) for x in ds) / len(ds), max(abs(x) for x in ds)))
