#!/bin/bash
# Per-case mesh audit for the dfvlr_r4 campaign row.
# Reads the authoritative nCells from each case's constant/polyMesh header note,
# then checks whether a time-dir polyMesh/points exists and what its checksum is.
cd /caseCampaign/campaign || exit 1

count_pts() { awk 'f==0 && /^[0-9]+$/ {print; exit}' "$1" 2>/dev/null; }

for M in 0.750 0.775 0.800 0.825 0.850 0.875 0.900 0.925 0.950 0.975 1.000 1.025 1.050 1.075 1.100; do
  C="dfvlr_r4_M${M}_a1.25"
  echo "===== $C"
  if [ -d "$C" ]; then
    echo "casemtime: $(date -r "$C" '+%Y-%m-%d %H:%M:%S')"
  fi
  OWN="$C/constant/polyMesh/owner"
  if [ -f "$OWN" ]; then
    echo "const_owner_mtime: $(date -r "$OWN" '+%Y-%m-%d %H:%M:%S')"
    grep -m1 -i 'note' "$OWN" | sed 's/^ *//'
    grep -m1 -oE 'nPoints:[0-9]+' "$OWN"
    grep -m1 -oE 'nCells:[0-9]+' "$OWN"
    grep -m1 -oE 'nFaces:[0-9]+' "$OWN"
    echo "const_owner_md5: $(md5sum "$OWN" | cut -d' ' -f1)"
  else
    echo "NO constant/polyMesh/owner"
  fi
  if [ -f "$C/constant/polyMesh/boundary" ]; then
    echo "const_boundary_note: $(grep -m1 -i 'note' "$C/constant/polyMesh/boundary" | sed 's/^ *//')"
  fi
  CP="$C/constant/polyMesh/points"
  if [ -f "$CP" ]; then
    echo "const_points: pts=$(count_pts "$CP") md5=$(md5sum "$CP" | cut -d' ' -f1) mtime=$(date -r "$CP" '+%Y-%m-%d %H:%M:%S')"
  else
    echo "const_points: MISSING"
  fi
  for TD in 0 0.1 0.15 0.2 0.25; do
    TP="$C/$TD/polyMesh/points"
    if [ -f "$TP" ]; then
      echo "time_${TD}_points: pts=$(count_pts "$TP") md5=$(md5sum "$TP" | cut -d' ' -f1) mtime=$(date -r "$TP" '+%Y-%m-%d %H:%M:%S')"
    fi
  done
  echo "polyMesh_dirs: $(ls -d "$C"/*/polyMesh 2>/dev/null | tr '\n' ' ')"
  # log lines
  for L in "$C"/log.checkMesh "$C"/log.*; do
    [ -f "$L" ] || continue
    CELLS=$(grep -m1 -E '^[[:space:]]*cells:' "$L" 2>/dev/null | sed 's/^ *//')
    if [ -n "$CELLS" ]; then
      echo "log: $(basename "$L"): $CELLS"
    fi
  done
  echo ""
done

echo "===== meshCache references"
ls -d /meshCache/*/ 2>/dev/null
for D in /meshCache/*/; do
  if [ -f "$D/points" ]; then
    echo "$D points md5=$(md5sum "$D/points" | cut -d' ' -f1) pts=$(count_pts "$D/points")"
  fi
done
echo "===== dfvlr_coarse (if present)"
ls -d /dfvlr_coarse* 2>/dev/null
for D in /dfvlr_coarse*/constant/polyMesh/points /dfvlr_coarse*/0/polyMesh/points; do
  [ -f "$D" ] && echo "$D md5=$(md5sum "$D" | cut -d' ' -f1) pts=$(count_pts "$D")"
done
echo AUDIT_DONE
