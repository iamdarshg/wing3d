import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from tsd_campaign import main
import json, os

# run remaining geometries one at a time, saving incrementally so
# interruption loses nothing
done = []
for g in ['rae2822', 'biconvex', 'dfvlr_r4']:
    p = os.path.join(r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign',
                     'tsd_%s.json' % g)
    if os.path.exists(p):
        print(g, 'already saved', flush=True)
        continue
    main(g)
    print(g, 'saved', flush=True)
print('ALL TSD DONE', flush=True)