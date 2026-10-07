import os
import shutil
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')

WING = 'D:\\CodeProjects\\cfd\\cases\\openfoam\\wing'
SPH = 'D:\\CodeProjects\\cfd\\cases\\openfoam\\sphere'
os.makedirs(os.path.join(SPH, 'system'), exist_ok=True)
os.makedirs(os.path.join(SPH, '0'), exist_ok=True)
os.makedirs(os.path.join(SPH, 'constant', 'triSurface'), exist_ok=True)

# proven dicts from wing container work
for f in ['snappyHexMeshDict', 'fvSchemes', 'fvSolution', 'controlDict',
          'meshQualityDict', 'decomposeParDict']:
    shutil.copy(os.path.join(WING, 'system_of', f),
                os.path.join(SPH, 'system', f))

# blockMesh: box around D=2 sphere
x0, x1, y0, y1, z0, z1 = -8, 16, -8, 8, -8, 8
nx, ny, nz = 30, 20, 20
open(os.path.join(SPH, 'system', 'blockMeshDict'), 'w').write(
    f"""FoamFile {{ version 2.0; format ascii; class dictionary; object blockMeshDict; }}
convertToMeters 1;
vertices (
 ({x0} {y0} {z0}) ({x1} {y0} {z0}) ({x1} {y1} {z0}) ({x0} {y1} {z0})
 ({x0} {y0} {z1}) ({x1} {y0} {z1}) ({x1} {y1} {z1}) ({x0} {y1} {z1}));
blocks (hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1));
edges ();
boundary (
 inlet {{ type patch; faces ((0 4 7 3)); }}
 outlet {{ type patch; faces ((1 2 6 5)); }}
 walls {{ type wall; faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }});
""")

# laminar transport + turbulence off
open(os.path.join(SPH, 'constant', 'transportProperties'), 'w').write(
    """FoamFile { version 2.0; format ascii; class dictionary; object transportProperties; }
transportModel Newtonian; nu [0 2 -1 0 0 0 0] 0.02;
""")
open(os.path.join(SPH, 'constant', 'turbulenceProperties'), 'w').write(
    """FoamFile { version 2.0; format ascii; class dictionary; object turbulenceProperties; }
simulationType laminar;
""")

# 0/ files (U=1, p, no turbulence)
open(os.path.join(SPH, '0', 'U'), 'w').write(
    """FoamFile { version 2.0; format ascii; class volVectorField; object U; }
dimensions [0 1 -1 0 0 0 0];
internalField uniform (1 0 0);
boundaryField {
 inlet { type fixedValue; value uniform (1 0 0); }
 outlet { type inletOutlet; inletValue uniform (1 0 0); value uniform (1 0 0); }
 walls { type slip; }
 sphere { type noSlip; }
}
""")
open(os.path.join(SPH, '0', 'p'), 'w').write(
    """FoamFile { version 2.0; format ascii; class volScalarField; object p; }
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField {
 inlet { type zeroGradient; }
 outlet { type fixedValue; value uniform 0; }
 walls { type zeroGradient; }
 sphere { type zeroGradient; }
}
""")
print('sphere case skeleton written')
