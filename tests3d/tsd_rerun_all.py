import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from tsd_campaign import main

for g in ['naca0012', 'rae2822', 'biconvex', 'dfvlr_r4']:
    main(g)
    print(g, 'DONE', flush=True)
print('ALL TSD DONE', flush=True)