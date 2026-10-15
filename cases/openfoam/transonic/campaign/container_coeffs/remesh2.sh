#!/bin/bash
# Fix locationInMesh off the z=0/y=0 faces for rae2822 + biconvex, re-mesh.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
for g in rae2822 biconvex; do
  src="${g}_M0.900_a1.25"
  sed -i 's/locationInMesh (5.5 0 0)/locationInMesh (5.5 0.3 0.25)/' "$src/system/snappyHexMeshDict"
  rm -rf "$src/constant/polyMesh" "$src/0.1*" "$src/0.14" "$src/0.28" "$src/postProcessing" 2>/dev/null
  cd "$src"
  blockMesh > log.blockMesh 2>&1
  snappyHexMesh -overwrite > log.snappy 2>&1
  checkMesh > log.check 2>&1
  echo "$g $(grep -a -m1 'cells:' log.check)"
  rm -f 0/cellLevel 0/pointLevel 2>/dev/null
  rm -rf /meshCache/$g; mkdir -p /meshCache/$g
  cp -r constant/polyMesh/* /meshCache/$g/ 2>/dev/null
  cd /caseCampaign/campaign
done
echo "DONE $(date +%H:%M:%S)"