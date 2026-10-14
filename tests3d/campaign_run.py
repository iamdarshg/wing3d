"""M0.85-1.1 high-fidelity campaign runner.

For each geometry (100k+ vertex STL) x Mach point: build a fine-mesh
OpenFOAM Euler case (rhoCentralFoam LTS, snappy from the STL), deploy
to the oftrans container, run in background, reap results into
cases/openfoam/transonic/campaign/.

Mesh: blockMesh 140x4x100 base + snappy level (3 4) -> ~120-150k cells
(host 395MB budget, keep <150k).
"""
import os
import shutil
import subprocess
import json
import numpy as np
import sys

sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.makecase import make_euler_central

HOST_BASE = r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign'
STL_DIR = r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl'
GEOS = ['naca0012', 'rae2822', 'biconvex', 'dfvlr_r4']
MACHES = [0.85, 0.875, 0.90, 0.925, 0.95, 0.975, 1.00, 1.025, 1.05, 1.075, 1.10]
ALPHA = 1.25


def write_blockmesh(root, stl_name):
    os.makedirs(os.path.join(root, 'system'), exist_ok=True)
    x0, x1 = -10, 21
    y0, y1 = -0.8, 0.8
    z0, z1 = -10, 10
    nx, ny, nz = 100, 3, 70
    with open(os.path.join(root, 'system', 'blockMeshDict'), 'w',
              newline='\n') as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; object blockMeshDict; }}
convertToMeters 1;
vertices (
 ({x0} {y0} {z0}) ({x1} {y0} {z0}) ({x1} {y1} {z0}) ({x0} {y1} {z0})
 ({x0} {y0} {z1}) ({x1} {y0} {z1}) ({x1} {y1} {z1}) ({x0} {y1} {z1}));
blocks (hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1));
edges ();
boundary (
 farfield {{ type patch; faces ((0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2)); }}
 symm {{ type symmetry; faces ((0 3 2 1) (4 5 6 7)); }});
""")


def write_snappy(root, stl_name):
    patch = stl_name.replace('.stl', '')
    with open(os.path.join(root, 'system', 'snappyHexMeshDict'), 'w',
              newline='\n') as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; object snappyHexMeshDict; }}
castellatedMesh true; snap true; addLayers true;
geometry {{ {stl_name} {{ type triSurfaceMesh; file "{stl_name}"; }}
 wakebox {{ type searchableBox; min (0 -0.8 -0.6); max (4 0.8 0.6); }} }}
castellatedMeshControls {{
 maxLocalCells 40000; maxGlobalCells 65000; minRefinementCells 10;
 maxLoadUnbalance 0.1; nCellsBetweenLevels 3;
 features ();
 refinementSurfaces {{ {stl_name} {{ level (3 4); patchInfo {{ type wall; }} }} }}
 resolveFeatureAngle 25;
 refinementRegions {{
  {stl_name} {{ mode distance; levels ((0.15 3) (0.4 2)); }}
  wakebox {{ mode inside; levels ((1E15 2)); }} }}
 locationInMesh (5.5 0 0);
 allowFreeStandingZoneFaces true;
}}
snapControls {{ nSmoothPatch 3; tolerance 2.0; nSolveIter 30; nRelaxIter 5; }}
addLayersControls {{
 relativeSizes false;
 layers {{ "{stl_name}" {{ nSurfaceLayers 4; }} }}
 expansionRatio 1.2; finalLayerThickness 2e-4; minThickness 1e-5;
 nGrow 0; featureAngle 60; slipFeatureAngle 30; nRelaxIter 5;
 nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10;
 maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3;
 minMedianAxisAngle 90; minMedialAxisAngle 90; nBufferCellsNoExtrude 0;
 nLayerIter 50; }}
meshQualityControls {{ #include "meshQualityDict" nSmoothScale 4; errorReduction 0.75; }}
writeFlags (scalarLevels layerSets);
mergeTolerance 1e-6;
""")


def make_case(host_dir, geo, M):
    root = os.path.join(host_dir, '%s_M%.3f_a%.2f' % (geo, M, ALPHA))
    shutil.rmtree(root, ignore_errors=True)
    for d in ('0', 'system', 'constant'):
        os.makedirs(os.path.join(root, d), exist_ok=True)
    # 0/ + templates via proven central writer (N2 gas consistent)
    make_euler_central(root, M, ALPHA, span_in=1.6, end_time=0.3)
    stl_name = geo + '.stl'
    # reference the STL: copy into constant/triSurface
    os.makedirs(os.path.join(root, 'constant', 'triSurface'), exist_ok=True)
    shutil.copy(os.path.join(STL_DIR, stl_name),
                os.path.join(root, 'constant', 'triSurface', stl_name))
    # controlDict patch name must match the stl basename
    cd = os.path.join(root, 'system', 'controlDict')
    txt = open(cd).read()
    txt = txt.replace('patches ("slab.stl")', 'patches ("%s")' % stl_name)
    txt = txt.replace('Aref 1.2', 'Aref 1.6')
    open(cd, 'w', newline='\n').write(txt)
    # 0/ field boundary patch names too
    for fn in ('U', 'p', 'T'):
        p = os.path.join(root, '0', fn)
        t = open(p).read()
        t = t.replace('"slab.stl"', '"%s"' % stl_name)
        open(p, 'w', newline='\n').write(t)
    write_blockmesh(root, stl_name)
    write_snappy(root, stl_name)
    # meshQualityDict include (snappy requires it)
    src = r'D:\CodeProjects\cfd\cases\openfoam\transonic\system\meshQualityDict'
    if os.path.exists(src):
        shutil.copy(src, os.path.join(root, 'system', 'meshQualityDict'))
    return root


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'naca0012'
    M = float(sys.argv[2]) if len(sys.argv) > 2 else 0.9
    root = make_case(HOST_BASE, which, M)
    print(root)
