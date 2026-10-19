#!/bin/bash
# lane fix: assemble coherent coarse mesh = good topology + good points (from /dfvlr_coarse)
source /usr/lib/openfoam/openfoam2406/etc/bashrc
FIX=/dfvlr_coarse
for M in 0.925 1.000; do
  L=/caseCampaign/probes/dfvlr_r4_M${M}_coarse45k
  SRC=/caseCampaign/cases/dfvlr_r4_M${M}_a1.25
  echo "===== LANE $M : $L ====="
  if [ ! -d "$L" ]; then
    echo "creating lane from $SRC"
    rm -rf "$L"
    mkdir -p "$L"
    cp -r "$SRC/0" "$L/0"
    cp -r "$SRC/system" "$L/system"
    cp -r "$SRC/constant" "$L/constant"
    rm -rf "$L/constant/polyMesh"
    mkdir -p "$L/constant/polyMesh"
  fi
  cd "$L" || { echo "LANE MISSING"; continue; }
  rm -rf 0/polyMesh
  for f in faces owner neighbour boundary cellZones faceZones pointZones level0Edge surfaceIndex; do
    cp "$FIX/constant/polyMesh/$f" constant/polyMesh/ 2>/dev/null
  done
  cp -r "$FIX/constant/polyMesh/sets" constant/polyMesh/ 2>/dev/null
  cp "$FIX/0/polyMesh/points" constant/polyMesh/points
  cp "$FIX/0/cellLevel" constant/polyMesh/cellLevel 2>/dev/null
  cp "$FIX/0/pointLevel" constant/polyMesh/pointLevel 2>/dev/null
  echo "-- checkMesh:"
  checkMesh > log.checkmesh 2>&1
  grep -a -E "points:|cells:|Mesh OK|Failed" log.checkmesh | head -8
done
echo "===== DONE ====="
