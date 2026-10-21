#!/bin/bash
# dfvlr_fine.sh — L5 follow-on: finer dfvlr mesh + high-M probe pair.
#
# L5 established that the delivered dfvlr anchors (uniform 86k mesh) are
# mesh-sensitive: a reconstructed 45.6k probe gives +30.5%/+19.0% CL at
# M0.925/M1.000. This program builds a FINER mesh (~250-450k cells) from
# the M0.900 template and re-runs the two ceiling points to t=0.25, to
# test whether the 86k anchors are mesh-converged.
#
# Run inside oftrans:  bash /tmp/dfvlr_fine.sh > /tmp/dfvlr_fine.log 2>&1
# Results: /caseCampaign/probes/dfvlr_fine_results.txt
#
# RAM-gated: waits for MemAvailable before meshing (>4 GB) and before
# each solve (>2.5 GB), because the host is shared.

source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign

freemb() { awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo; }
wait_ram() {  # $1 = threshold MB, $2 = stage label
  while [ "$(freemb)" -lt "$1" ]; do
    echo "$(date +%H:%M:%S) [$2] low RAM $(freemb)MB < $1 MB - waiting 5 min"
    sleep 300
  done
}

TEMPLATE=campaign/dfvlr_r4_M0.900_a1.25
BASE=probes/dfvlr_fine

echo "=== stage 1: fine mesh from template ==="
wait_ram 4000 mesh
rm -rf "$BASE"; mkdir -p "$BASE"
cp -r "$TEMPLATE/0" "$TEMPLATE/constant" "$TEMPLATE/system" "$BASE/"
cd "$BASE" || exit 2
sed -i -e 's/maxLocalCells 40000/maxLocalCells 250000/' \
       -e 's/maxGlobalCells 65000/maxGlobalCells 450000/' \
       -e 's/level (3 4)/level (4 5)/' \
       -e 's/levels ((0.15 3) (0.4 2))/levels ((0.15 4) (0.4 3))/' \
       -e 's/nCellsBetweenLevels 3/nCellsBetweenLevels 4/' \
    system/snappyHexMeshDict
rm -rf constant/polyMesh 0/polyMesh 0/cellLevel 0/pointLevel 0.[0-9]* postProcessing
blockMesh > log.blockMesh 2>&1
snappyHexMesh -overwrite > log.snappy 2>&1
checkMesh > log.check 2>&1
CELLS=$(grep -a -m1 'cells:' log.check | grep -ao '[0-9]\+' | head -1)
echo "FINE MESH cells: $CELLS  (86k reference)"
if [ -z "$CELLS" ] || [ "$CELLS" -lt 150000 ] || [ "$CELLS" -gt 700000 ]; then
  echo "cell count out of band ($CELLS) - aborting"; exit 2
fi

echo "=== stage 2: per-M cases (U/p/T copied from the delivered same-M cases) ==="
cd /caseCampaign/probes || exit 2
for M in 0.925 1.000; do
  rm -rf "dfvlr_fine_M${M}"
  cp -r dfvlr_fine "dfvlr_fine_M${M}"
  cd "dfvlr_fine_M${M}" || exit 2
  rm -rf 0.[0-9]* postProcessing log.*
  cp "/caseCampaign/campaign/dfvlr_r4_M${M}_a1.25/0/U" 0/U
  cp "/caseCampaign/campaign/dfvlr_r4_M${M}_a1.25/0/p" 0/p
  cp "/caseCampaign/campaign/dfvlr_r4_M${M}_a1.25/0/T" 0/T
  cd /caseCampaign/probes
done

echo "=== stage 3: solves (sequential, RAM-gated) ==="
: > /tmp/dfvlr_fine_batch.log
run_case() {
  c=$1
  cd "/caseCampaign/probes/$c" || exit 2
  rm -rf 0.[0-9]* postProcessing
  wait_ram 2500 "$c"
  echo "$(date +%H:%M:%S) SOLVE START $c" >> /tmp/dfvlr_fine_batch.log
  rhoCentralFoam > log.run 2>&1
  f=$(find postProcessing -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ]; then
    last=$(tail -n 1 "$f"); set -- $last
    echo "$c done t=$1 CL=$5 CD=$2" >> /tmp/dfvlr_fine_batch.log
  else
    echo "$c NOFORCES" >> /tmp/dfvlr_fine_batch.log
  fi
}
run_case dfvlr_fine_M0.925
run_case dfvlr_fine_M1.000

echo "=== stage 4: comparison vs the 86k anchors ==="
{
  echo "# dfvlr fine-mesh probe results $(date -u +%FT%TZ)"
  echo "# mesh: $(cd /caseCampaign/probes/dfvlr_fine && grep -a -m1 'cells:' log.check)"
  echo "# anchors (86k): M0.925 CL=0.2779525  M1.000 CL=0.2126116"
  for M in 0.925 1.000; do
    f=$(find "/caseCampaign/probes/dfvlr_fine_M${M}/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
    if [ -n "$f" ]; then
      cl=$(tail -n 1 "$f" | awk '{print $5}')
      t=$(tail -n 1 "$f" | awk '{print $1}')
      echo "dfvlr_fine_M${M}  t=${t}  CL=${cl}"
    else
      echo "dfvlr_fine_M${M}  NOFORCES"
    fi
  done
} > /caseCampaign/probes/dfvlr_fine_results.txt
cat /caseCampaign/probes/dfvlr_fine_results.txt
echo "DFVLR-FINE DONE $(date +%H:%M:%S)"
