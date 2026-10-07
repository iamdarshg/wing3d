import os
import shutil
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')

WING = 'D:\\CodeProjects\\cfd\\cases\\openfoam\\wing'
CAR = 'D:\\CodeProjects\\cfd\\cases\\openfoam\\car'

# 1. proven system dicts -> car
for f in ['snappyHexMeshDict', 'fvSchemes', 'fvSolution', 'controlDict',
          'meshQualityDict', 'forceCoeffs', 'decomposeParDict',
          'blockMeshDict']:
    shutil.copy(os.path.join(WING, 'system_of', f),
                os.path.join(CAR, 'system', f))

# 2. 0/ files (SST) -> car, with car inflow (alpha 0, U=30)
for f in ['U', 'p', 'k', 'omega', 'nut']:
    shutil.copy(os.path.join(WING, '0', f), os.path.join(CAR, '0', f))
for extra in ['nuTilda']:
    p = os.path.join(CAR, '0', extra)
    if os.path.exists(p):
        os.remove(p)

U30 = '(30 0 0)'
p = os.path.join(CAR, '0', 'U')
s = open(p).read()
s = s.replace('(29.927 0 2.0927)', U30)
open(p, 'w').write(s)

# 3. adapt dicts for car
p = os.path.join(CAR, 'system', 'snappyHexMeshDict')
s = open(p).read()
s = s.replace('wing.stl', 'car.stl')
s = s.replace('name wing;', 'name car;')
s = s.replace('min  (-1.0 -3.5 -1.0);', 'min  (-1.0 -1.0 -0.5);')
s = s.replace('max  ( 4.0 3.5 1.0);', 'max  (4.0 1.0 1.0);')
s = s.replace('locationInMesh (0.11 6.11 0.13);',
              'locationInMesh (0.0 2.0 0.0);')
s = s.replace('"wing.*"', '"car.*"')
open(p, 'w').write(s)

p = os.path.join(CAR, 'system', 'blockMeshDict')
x0, x1, y0, y1, z0, z1 = -5, 11, -4, 4, -2, 3
nx, ny, nz = 40, 30, 30
open(p, 'w').write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; object blockMeshDict; }}
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

# controlDict forceCoeffs inline values
for f in ['controlDict']:
    p = os.path.join(CAR, 'system', f)
    s = open(p).read()
    s = s.replace('patches ("wing")', 'patches ("car")')
    s = s.replace('CofR (0.25 0 0)', 'CofR (0.5 0 0)')
    s = s.replace('magUInf 30.0', 'magUInf 30')
    s = s.replace('lRef 1.0', 'lRef 1.044')
    s = s.replace('Aref 6.0', 'Aref 0.112')
    open(p, 'w').write(s)
print('car system adapted')
