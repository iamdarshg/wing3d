#!/bin/bash
# Stop the invalid-geometry runs and re-mesh the 3 corrected geometries.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
pkill -9 -f solve_worker 2>/dev/null
pkill -9 -f dfvlr_runner 2>/dev/null
pkill -9 rhoCentralFoam 2>/dev/null
sleep 2
cd /caseCampaign/campaign
# clear done markers for the 3 invalid geometries + naca stays valid
for g in rae2822 biconvex dfvlr_r4; do
  rm -f results/${g}_M*.done results/${g}_M*.dat
  rm -f /meshCache/$g -rf
done
echo "stopped + cleared. re-meshing with corrected STLs"
for g in rae2822 biconvex dfvlr_r4; do
  src="${g}_M0.900_a1.25"
  rm -rf "$src/0.1*" "$src/0.14" "$src/0.28" "$src/postProcessing" "$src/constant/polyMesh" 2>/dev/null
  cd "$src"
  blockMesh > log.blockMesh 2>&1
  snappyHexMesh -overwrite > log.snappy 2>&1
  checkMesh > log.check 2>&1
  echo "$g $(grep -a -m1 'cells:' log.check)"
  rm -f 0/cellLevel 0/pointLevel 2>/dev/null
  mkdir -p /meshCache/$g
  cp -r constant/polyMesh/* /meshCache/$g/ 2>/dev/null
  cd /caseCampaign/campaign
done
echo "RE-MESH DONE $(date +%H:%M:%S)"