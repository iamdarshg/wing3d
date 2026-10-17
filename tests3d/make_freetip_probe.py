"""Make the naca0012 M0.850 free-tip OF probe case.

Mirrors campaign_run.make_case exactly except:
- blockMesh y: -2.4..2.4, ny=9 (same dy=0.5333 as the corpus)
- the two y-faces join the 'farfield' patch (corpus: separate 'symm'
  symmetry patch = wall-at-tip); this makes a true free-tip slab
- snappy wakebox y +-2.4; maxLocalCells/maxGlobalCells raised for the
  larger base mesh
- controlDict force normalization copied verbatim from the corpus case
  (magUInf/rhoInf/Aref) so CL is directly comparable.
"""
import os
import re
import shutil
import sys

sys.path.insert(0, r'D:\CodeProjects\cfd')
from wing3d.makecase import make_euler_central

HOST = r'D:\CodeProjects\cfd\cases\openfoam\transonic'
ROOT = os.path.join(HOST, 'campaign_freetip', 'naca0012_M0.850_a1.25')
CORPUS = os.path.join(HOST, 'campaign', 'naca0012_M0.850_a1.25')


def main():
    shutil.rmtree(ROOT, ignore_errors=True)
    for d in ('0', 'system', 'constant'):
        os.makedirs(os.path.join(ROOT, d), exist_ok=True)
    make_euler_central(ROOT, 0.85, 1.25, span_in=1.6, end_time=0.25)
    os.makedirs(os.path.join(ROOT, 'constant', 'triSurface'), exist_ok=True)
    shutil.copy(os.path.join(HOST, 'stl', 'naca0012.stl'),
                os.path.join(ROOT, 'constant', 'triSurface', 'naca0012.stl'))
    # controlDict: patch name + Aref, then force-normalization keys verbatim
    cd = os.path.join(ROOT, 'system', 'controlDict')
    txt = open(cd).read()
    txt = txt.replace('patches ("slab.stl")', 'patches ("naca0012.stl")')
    txt = txt.replace('Aref 1.2', 'Aref 1.6')
    open(cd, 'w', newline='\n').write(txt)
    corpus_cd = open(os.path.join(CORPUS, 'system', 'controlDict')).read()
    for key in ('magUInf', 'rhoInf', 'Aref'):
        m = re.search(r'%s\s+([0-9.eE+-]+)' % key, corpus_cd)
        if m:
            txt = open(cd).read()
            txt = re.sub(r'%s\s+[0-9.eE+-]+' % key,
                         '%s %s' % (key, m.group(1)), txt)
            open(cd, 'w', newline='\n').write(txt)
    # 0/ fields: drop the symm entry, keep farfield
    for fn in ('U', 'p', 'T'):
        p = os.path.join(ROOT, '0', fn)
        t = open(p).read()
        t = t.replace('"slab.stl"', '"naca0012.stl"')
        t = t.replace(' symm { type symmetry; }\n', '')
        assert 'symm' not in t, fn
        open(p, 'w', newline='\n').write(t)
    # blockMeshDict: y +-2.4, ny=9, y-faces merged into farfield
    bm = (
        'FoamFile { version 2.0; format ascii; class dictionary; '
        'object blockMeshDict; }\n'
        'convertToMeters 1;\n'
        'vertices (\n'
        ' (-10 -2.4 -10) (21 -2.4 -10) (21 2.4 -10) (-10 2.4 -10)\n'
        ' (-10 -2.4 10) (21 -2.4 10) (21 2.4 10) (-10 2.4 10));\n'
        'blocks (hex (0 1 2 3 4 5 6 7) (100 9 70) simpleGrading (1 1 1));\n'
        'edges ();\n'
        'boundary (\n'
        ' farfield { type patch; faces ((0 4 7 3) (1 2 6 5) (0 1 5 4) '
        '(3 7 6 2) (0 3 2 1) (4 5 6 7)); });\n')
    open(os.path.join(ROOT, 'system', 'blockMeshDict'), 'w',
         newline='\n').write(bm)
    # snappy: reuse the corpus file, widen wakebox + raise caps
    sp_path = os.path.join(CORPUS, 'system', 'snappyHexMeshDict')
    sp = open(sp_path).read()
    sp = sp.replace('maxLocalCells 40000; maxGlobalCells 65000;',
                    'maxLocalCells 80000; maxGlobalCells 95000;')
    sp = sp.replace('min (0 -0.8 -0.6); max (4 0.8 0.6);',
                    'min (0 -2.4 -0.6); max (4 2.4 0.6);')
    open(os.path.join(ROOT, 'system', 'snappyHexMeshDict'), 'w',
         newline='\n').write(sp)
    src = os.path.join(HOST, 'system', 'meshQualityDict')
    if os.path.exists(src):
        shutil.copy(src, os.path.join(ROOT, 'system', 'meshQualityDict'))
    print('made', ROOT)


if __name__ == '__main__':
    main()
