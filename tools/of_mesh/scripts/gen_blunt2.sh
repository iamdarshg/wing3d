#!/bin/bash
set -e
cd /work/tmr/wing_release/wing_release_072319
rm -f wing_strct.1.p3d wing_strct.1.su2
sed -i "s/tip_type = 'round'/tip_type = 'blunt'/" input.nml
sed -i 's/use_s_blunt = T/use_s_blunt = F/' input.nml
sed -i 's/igrid_type = 4/igrid_type = 5/' input.nml
grep -E 'tip_type|igrid_type' input.nml
./hcf_wing | grep -E 'Root chord|Aspect|Reference area|success'
ls -la wing_strct.1.p3d 2>/dev/null || echo "NO P3D"

source /usr/lib/openfoam/openfoam2406/etc/bashrc
rm -rf /work/naca_blunt
mkdir -p /work/naca_blunt/system /work/naca_blunt/constant
cd /work/naca_blunt
cp /work/tmr/wing_release/wing_release_072319/wing_strct.1.p3d .
cp /work/naca_su2_case3/system/fvSchemes /work/naca_su2_case3/system/fvSolution \
   /work/naca_su2_case3/system/controlDict /work/naca_su2_case3/system/fvOptions system/
plot3dToFoam -noBlank wing_strct.1.p3d > log.convert 2>&1
autoPatch 30 -overwrite > log.patch 2>&1
checkMesh > log.check 2>&1 || true
grep -E 'cells:|Max aspect|non-orthogonality Max|Max skewness|severely non|Mesh OK|Failed' log.check | head -8
