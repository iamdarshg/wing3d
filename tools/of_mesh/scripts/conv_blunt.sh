#!/bin/bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
rm -rf /work/naca_blunt
mkdir -p /work/naca_blunt/system /work/naca_blunt/constant
cd /work/naca_blunt || exit 1
cp /work/tmr/wing_release/wing_release_072319/wing_strct.1.p3d .
cp /work/naca_su2_case3/system/fvSchemes /work/naca_su2_case3/system/fvSolution \
   /work/naca_su2_case3/system/controlDict /work/naca_su2_case3/system/fvOptions system/
plot3dToFoam -noBlank wing_strct.1.p3d > log.convert 2>&1
autoPatch 30 -overwrite > log.patch 2>&1
checkMesh > log.check 2>&1
grep -E 'cells:|Max aspect|non-orthogonality Max|Max skewness|severely non|Mesh OK|Failed' log.check | head -8
