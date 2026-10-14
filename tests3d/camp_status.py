"""Clean status of the OF campaign by reading coefficient.dat directly."""
import os, glob

CAMP = '/caseCampaign/campaign'
GEOS = ['naca0012', 'rae2822', 'biconvex', 'dfvlr_r4']

def parse(name):
    f = glob.glob(os.path.join(CAMP, name, 'postProcessing', '*',
                               'coefficient.dat'))
    f = [p for p in f if os.path.getsize(p) > 100]
    if not f:
        return None
    rows = [l.split() for l in open(f[0]).read().splitlines()
            if l.strip() and not l.lstrip().startswith('#')]
    if not rows:
        return None
    t = float(rows[-1][0]); cl = float(rows[-1][4]); cd = float(rows[-1][6])
    return t, cl, cd

print('%-11s %-6s %-8s %-9s %-9s %s' % ('geo', 'M', 't', 'CL', 'CD', 'status'))
for g in GEOS:
    for d in sorted(glob.glob(os.path.join(CAMP, g + '_M*'))):
        name = os.path.basename(d)
        M = name.split('_M')[1].split('_')[0]
        r = parse(name)
        if r is None:
            print('%-11s %-6s %-8s %-9s %-9s %s' % (g, M, '-', '-', '-', 'no-data'))
        else:
            t, cl, cd = r
            st = 'CONV' if t > 0.2 else ('running' if t > 0.01 else 'startup')
            print('%-11s %-6s %-8.3f %-9.5f %-9.5f %s' % (g, M, t, cl, cd, st))