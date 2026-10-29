#!/bin/bash
# (1) orthogonality tuning: production-TMR-like radial/farfield params
source /usr/lib/openfoam/openfoam2406/etc/bashrc >/dev/null 2>&1
G=/work/tmr/wing_release/wing_release_072319
cd "$G" || exit 1
rm -f wing_strct.1.p3d

python3 - "$G/input.nml" <<'PY' 2>/dev/null || true
PY

# edit params with sed (container has no python3)
sed -i 's/^ *nnodes_cylinder_input = .*/           nnodes_cylinder_input = 48/' input.nml
sed -i 's/^ *nre = .*/                             nre = 64/' input.nml
sed -i 's/^ *R_outer = .*/                         R_outer = 100/' input.nml
sed -i 's/^ *nr_gs = .*/                           nr_gs = 8/' input.nml
sed -i 's/^ *target_reynolds_number = .*/          target_reynolds_number = 6.0e6/' input.nml
sed -i "s/^ *tip_type = .*/                        tip_type = 'blunt'/" input.nml
sed -i 's/^ *use_s_blunt = .*/                     use_s_blunt = F/' input.nml
sed -i 's/^ *igrid_type = .*/                      igrid_type = 5/' input.nml
sed -i 's/ugrid_file_unformatted = T/ugrid_file_unformatted = F/' input.nml
sed -i 's/^ *stretching_tanh_towards_lete = .*/    stretching_tanh_towards_lete = 4.5/' input.nml

echo "--- params ---"
grep -E 'nnodes_cylinder_input|nre =|R_outer|tips_type|tip_type|igrid_type|target_reynolds' input.nml

./hcf_wing | grep -E 'Reference area|Aspect|success'
ls -la wing_strct.1.p3d 2>/dev/null || echo "NO P3D"

rm -rf /work/tune1
mkdir -p /work/tune1/system /work/tune1/constant
cd /work/tune1 || exit 1
cp "$G"/wing_strct.1.p3d .
cp /work/naca_su2_case3/system/fvSchemes /work/naca_su2_case3/system/fvSolution \
   /work/naca_su2_case3/system/controlDict /work/naca_su2_case3/system/fvOptions system/
plot3dToFoam -noBlank wing_strct.1.p3d > log.convert 2>&1
autoPatch 30 -overwrite > log.patch 2>&1
checkMesh > log.check 2>&1
echo "--- checkMesh ---"
grep -E 'cells:|Max aspect|non-orthogonality Max|Max skewness|severely non|Mesh OK|Failed' log.check | head -8
