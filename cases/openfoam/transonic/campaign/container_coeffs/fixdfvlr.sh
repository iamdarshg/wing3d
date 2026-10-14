#!/bin/bash
# Fix DFVLR: rebuild the 54k mesh in dfvlr_r4_M0.900_a1.25 (already has
# correct snappy), then copy its polyMesh into ALL dfvlr cases + harden
# numerics (maxCo 0.3, deltaT 1e-5) to survive the sharp LE transient.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign

# 1) regenerate the dfvlr mesh from scratch in a scratch dir with the
#    proven snappy settings + finer LE level (4 5) to resolve the LE.
SCR=/dfvlr_mesh
rm -rf $SCR; mkdir -p $SCR/constant/triSurface $SCR/system $SCR/0
cp dfvlr_r4_M0.900_a1.25/constant/triSurface/dfvlr_r4.stl $SCR/constant/triSurface/
cp dfvlr_r4_M0.900_a1.25/system/blockMeshDict $SCR/system/
cp dfvlr_r4_M0.900_a1.25/system/meshQualityDict $SCR/system/
# refine LE: level (4 5) + tighter distance zones
sed -e 's/level (3 4)/level (4 5)/' \
    -e 's/levels ((0.15 3) (0.4 2))/levels ((0.10 4) (0.3 3))/' \
    dfvlr_r4_M0.900_a1.25/system/snappyHexMeshDict > $SCR/system/snappyHexMeshDict
# need controlDict/fvSchemes/fvSolution for snappy to run
cp dfvlr_r4_M0.900_a1.25/system/controlDict $SCR/system/
cp dfvlr_r4_M0.900_a1.25/system/fvSchemes $SCR/system/
cp dfvlr_r4_M0.900_a1.25/system/fvSolution $SCR/system/
cp dfvlr_r4_M0.900_a1.25/0/* $SCR/0/ 2>/dev/null
cp dfvlr_r4_M0.900_a1.25/constant/thermophysicalProperties $SCR/constant/ 2>/dev/null
cp dfvlr_r4_M0.900_a1.25/constant/turbulenceProperties $SCR/constant/ 2>/dev/null

cd $SCR
blockMesh > log.blockMesh 2>&1
snappyHexMesh -overwrite > log.snappy 2>&1
checkMesh > log.check 2>&1
echo "dfvlr mesh: $(grep -a -m1 'cells:' log.check)"
rm -f 0/cellLevel 0/pointLevel 2>/dev/null

# 2) copy the new polyMesh into every dfvlr case
cd /caseCampaign/campaign
for c in dfvlr_r4_M*; do
  rm -rf "$c/constant/polyMesh"
  cp -r $SCR/constant/polyMesh "$c/constant/polyMesh"
  # harden numerics: maxCo 0.3, deltaT 1e-5, maxDeltaT 5e-5
  sed -i 's/maxCo 0.5/maxCo 0.3/; s/deltaT 2e-5/deltaT 1e-5/; s/maxDeltaT 1e-4/maxDeltaT 5e-5/' "$c/system/controlDict"
  # clear stale results
  rm -f "results/$c.done" "results/$c.dat"
  rm -rf "$c/0.1*" "$c/0.14" "$c/0.28" "$c/postProcessing" 2>/dev/null
done
echo "DFVLR FIXED: $(ls dfvlr_r4_M*/constant/polyMesh | wc -l) cases updated"
echo "numerics now: $(grep -h -E 'maxCo|deltaT|maxDeltaT' dfvlr_r4_M1.000_a1.25/system/controlDict | tr '\n' ' ')"