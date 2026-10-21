#!/bin/bash
# dfvlr_fine.sh — L5 follow-on: finer dfvlr mesh + high-M probe pair.
#
# L5 established that the delivered dfvlr anchors (uniform 86k mesh) are
# mesh-sensitive: a reconstructed 45.6k probe gives +30.5%/+19.0% CL at
# M0.925/M1.000. This program builds a FINER mesh (~250-450k cells) from
# the M0.900 template and re-runs the WHOLE dfvlr band M0.850-1.050
# (9 points: the campaign regime that has delivered 86k anchors) to
# t=0.25, to test whether the 86k anchors are mesh-converged and to
# give the mesh-correction band as a function of M.
#
# Run inside oftrans:  bash /tmp/dfvlr_fine.sh > /tmp/dfvlr_fine.log 2>&1
# Results: /caseCampaign/probes/dfvlr_fine_results.txt
#
# RAM-gated: waits for MemAvailable before meshing (>2 GB) and before
# each solve (>1.8 GB). Sized for the WSL VM the host set at 3 GB RAM
# + 4 GB swap (2026-10-21): full snappy caps (450k cells) restored,
# solve RSS ~5 GB with the overflow paged.

# single-instance lock: a second invocation exits immediately
LOCK=/tmp/dfvlr_fine.lock
if [ -f "$LOCK" ] && kill -0 "$(cat "$LOCK" 2>/dev/null)" 2>/dev/null; then
  echo "dfvlr_fine already running (pid $(cat "$LOCK")) - exiting"
  exit 0
fi
echo $$ > "$LOCK"
trap 'rm -f "$LOCK"' EXIT

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
wait_ram 2000 mesh
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
for M in 0.850 0.875 0.900 0.925 0.950 0.975 1.000 1.025 1.050; do
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
  wait_ram 1800 "$c"
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
for M in 0.850 0.875 0.900 0.925 0.950 0.975 1.000 1.025 1.050; do
  run_case "dfvlr_fine_M${M}"
done

echo "=== stage 4: comparison vs the 86k anchors + TSD ==="
anchor_for() { case $1 in
  0.850) echo 0.2674053;; 0.875) echo 0.2729032;; 0.900) echo 0.2791015;;
  0.925) echo 0.2779525;; 0.950) echo 0.2688365;; 0.975) echo 0.2460678;;
  1.000) echo 0.2126116;; 1.025) echo 0.1725623;; 1.050) echo 0.1334340;; esac; }
tsd_for() { case $1 in
  0.850) echo 0.5036366438470022;; 0.875) echo 0.49591162975187797;;
  0.900) echo 0.4515527059358044;; 0.925) echo 0.47761650712854;;
  0.950) echo 0.4764804179702282;; 0.975) echo 0.4935046159921045;;
  1.000) echo 0.5036980989060582;; 1.025) echo 0.4823874197719977;;
  1.050) echo 0.5018454086784541;; esac; }
{
  echo "# dfvlr fine-mesh band results $(date -u +%FT%TZ)"
  echo "# mesh: $(cd /caseCampaign/probes/dfvlr_fine && grep -a -m1 'cells:' log.check)"
  echo "# columns: M t fine_CL anchor86k mesh_corr% tsd_cl err_vs_fine% err_asrun%"
  for M in 0.850 0.875 0.900 0.925 0.950 0.975 1.000 1.025 1.050; do
    f=$(find "/caseCampaign/probes/dfvlr_fine_M${M}/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
    if [ -n "$f" ]; then
      tail -n 1 "$f" | awk -v m="$M" -v a="$(anchor_for "$M")" -v s="$(tsd_for "$M")" \
        '{t=$1; cl=$5; mc=(cl-a)/a*100; e1=(s-cl)/cl*100; e2=(s-a)/a*100;
          printf "%s  t=%s  fine_CL=%.7f  anchor86k=%.7f  mesh_corr=%+.2f%%  tsd_cl=%.7f  err_vs_fine=%+.2f%%  err_asrun=%+.2f%%\n", m, t, cl, a, mc, s, e1, e2}'
    else
      echo "dfvlr_fine_M${M}  NOFORCES"
    fi
  done
} > /caseCampaign/probes/dfvlr_fine_results.txt
cat /caseCampaign/probes/dfvlr_fine_results.txt
echo "DFVLR-FINE DONE $(date +%H:%M:%S)"
