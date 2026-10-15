#!/bin/bash
# v2 worker: dfvlr + biconvex on the PROVEN .old meshes (same sliced geometry).
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
mkdir -p results2
CASES="dfvlr_r4_M0.750_a1.25 dfvlr_r4_M0.775_a1.25 dfvlr_r4_M0.800_a1.25 dfvlr_r4_M0.825_a1.25 dfvlr_r4_M0.900_a1.25 dfvlr_r4_M1.075_a1.25 dfvlr_r4_M1.100_a1.25 biconvex_M0.750_a1.25 biconvex_M0.775_a1.25 biconvex_M0.800_a1.25 biconvex_M0.825_a1.25 biconvex_M0.850_a1.25"
for c in $CASES; do
  if [ -f "results2/$c.done2" ]; then echo "skip $c" >> batchv2.log; continue; fi
  g=$(echo "$c" | sed 's/_M.*//')
  CPATH="/meshCache/$g.old"
  if [ ! -d "$CPATH" ]; then echo "NOCACHE $c" >> batchv2.log; continue; fi
  echo "$(date +%H:%M:%S) start $c (mesh=$g.old)" >> batchv2.log
  rm -rf "$c/0.1"* "$c/0.2"* "$c/postProcessing" "$c/constant/polyMesh" "$c/0/polyMesh" 2>/dev/null
  cp -r "$CPATH" "$c/constant/polyMesh"
  cd "$c"
  rhoCentralFoam > log.run 2>&1
  cd ..
  f=$(find "$c/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    tail -n 1 "$f" > "results2/$c.dat"
    echo done > "results2/$c.done"
    echo done > "results2/$c.done2"
    echo "$(date +%H:%M:%S) DONE $c $(tail -1 $f | cut -f1,5,7)" >> batchv2.log
  else
    touch "results2/$c.dat"; echo no-forces > "results2/$c.done2"
    echo "$(date +%H:%M:%S) NOFORCES $c" >> batchv2.log
  fi
done
echo "V2 WORKER FINISHED $(date +%H:%M:%S)" >> batchv2.log
