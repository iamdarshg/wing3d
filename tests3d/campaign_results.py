"""Aggregate campaign results + compare vs TSD and legacy anchors."""
import os
import glob
import json
import numpy as np

CAMP = r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign'


def fmt(x):
    return ('%.5f' % x) if x is not None else '-'

def parse_case(name):
    """Read coefficient.dat final row -> (CL, CD). Columns:
    time, Cm?, Cd_l, Cd_p, Cd_f?, Cl, ... (12 cols). The line:
    t  Cd?  ..  Cl ..  CD  ..  (need column map: from header we know
    col1=Cd(?) Actually OF forceCoeffs layout: t CmRoll CmYaw CmPitch Cs Cl Cd ...)
    We sample header + last row and pick by known layout:
    OF2406 coefficient.dat fields: t, Cm... the classic:
    time Cm Cm-f Cm-p ? ... safest: use log.run Cl:/Cd: lines.
    """
    log = os.path.join(CAMP, name, 'log.run')
    if not os.path.exists(log):
        return None
    cl = cd = None
    with open(log, 'r', errors='ignore') as f:
        for line in f:
            if line.strip().startswith('Cl:'):
                cl = float(line.split()[1])
            if line.strip().startswith('Cd:'):
                cd = float(line.split()[1])
    return cl, cd


def parse_coeff_dat(name):
    """Parse coefficient.dat directly (last row). Layout (12 cols):
    t  CmRoll  CmPitch  CmYaw  Cl  Cd  ... -> col4=Cl col5=Cd? verify
    with header-free rows: biconvex M0.95: ...1.089e-01 2.953e-02 7.94e-02
    => 0.1089 0.0295 0.0794 -> looks like Cl=0.1089 Cd=0.0794? too high.
    Use log.run Cl/Cd (authoritative)."""
    return parse_case(name)


def main():
    rows = []
    for d in sorted(glob.glob(os.path.join(CAMP, '*_M*.dat'))):
        name = os.path.basename(d)[:-4]
        r = parse_case(name)
        if r is None:
            continue
        geo = name.split('_M')[0]
        M = float(name.split('_M')[1].split('_')[0])
        rows.append((geo, M, r[0], r[1]))
    print('%-11s %-6s %-10s %-10s' % ('geo', 'M', 'CL', 'CD'))
    by_geo = {}
    for geo, M, cl, cd in rows:
        by_geo.setdefault(geo, []).append((M, cl, cd))
    return by_geo


if __name__ == '__main__':
    g = main()
    for geo, lst in g.items():
        print(geo)
        for M, cl, cd in sorted(lst):
            print('  M=%.3f CL=%.4f CD=%.4f' % (M, cl, cd))