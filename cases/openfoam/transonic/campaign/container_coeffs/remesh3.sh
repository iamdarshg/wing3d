#!/bin/bash
# Remove legacy 0/polyMesh residue so snappy uses constant/polyMesh.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
for g in rae2822 biconvex; do
  src="${g}_M0.900_a1.25"
  rm -rf "$src/0/polyMesh" "$src/0/cellLevel" "$src/0/pointLevel"
  rm -rf "$src/constant/polyMesh" "$src/0.1*" "$src/0.14" "$src/0.28" "$src/postProcessing" 2>/dev/null
  cd "$src"
  blockMesh > log.blockMesh 2>&1
  echo -n "blockMesh: "; grep -a -c 'Creating block' log.blockMesh
  snappyHexMesh -overwrite > log.snappy 2>&1
  checkMesh > log.check 2>&1
  echo "$g $(grep -a -m1 'cells:' log.check)"
  rm -f 0/cellLevel 0/pointLevel 2>/dev/null
  rm -rf /meshCache/$g; mkdir -p /meshCache/$g
  cp -r constant/polyMesh/* /meshCache/$g/ 2>/dev/null
  cd /caseCampaign/campaign
done
echo "DONE $(date +%H:%M:%S)"