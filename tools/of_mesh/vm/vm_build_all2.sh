#!/bin/bash
# Build all free-tip cases (absolute paths), then solve in parallel.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
REPO=$HOME/cfd
FT=$REPO/cases/openfoam/transonic/freetip
STL=$REPO/cases/openfoam/transonic/stl

CASES="dfvlr_r4_M1.000_a1.25 rae2822_M0.850_a1.25 rae2822_M0.925_a1.25 rae2822_M1.000_a1.25 naca0012_M0.800_a1.25_ny60 naca0012_M0.800_a1.25_ny30 dfvlr_r4_M0.850_a1.25 dfvlr_r4_M0.925_a1.25"

echo "=== BUILD ==="
for d in $CASES; do
  cd "$FT/$d" 2>/dev/null || { echo "MISSING $d"; continue; }
  stl=$(grep -oE 'file +"[^"]+"' system/snappyHexMeshDict | head -1 | sed 's/.*"\(.*\)"/\1/')
  mkdir -p constant/triSurface
  [ -f "$STL/$stl" ] && cp -f "$STL/$stl" constant/triSurface/
  if [ ! -f constant/polyMesh/points ]; then
    blockMesh > log.blockMesh 2>&1
    snappyHexMesh -overwrite > log.snappy 2>&1
  fi
  checkMesh > log.checkMesh 2>&1
  n=$(grep -m1 'cells:' log.checkMesh | tr -s ' ' | sed 's/^ *//')
  ok=$(grep -c 'Mesh OK' log.checkMesh)
  printf '%-34s %-20s meshOK=%s\n' "$d" "$n" "$ok"
done

echo "=== SOLVE (parallel) ==="
for d in $CASES; do
  cd "$FT/$d" 2>/dev/null || continue
  rm -f log.gcp
  nohup rhoCentralFoam > log.gcp 2>&1 &
  printf 'launched %-34s pid=%s\n' "$d" "$!"
done
sleep 20
echo "running: $(pgrep -c rhoCentralFoam 2>/dev/null || echo 0)"
for d in $CASES; do
  printf '%-34s %s\n' "$d" "$(tail -1 $FT/$d/log.gcp 2>/dev/null | cut -c1-70)"
done
