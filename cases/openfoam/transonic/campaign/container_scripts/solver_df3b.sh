#!/bin/bash
# v3b worker: dfvlr in-band cases on the golden 86k mesh, numerics matched
# to the sibling corpus (maxCo 0.5, maxDeltaT 1e-4). Out-of-band 1.075/1.100
# dropped (old partials, outside deliverable band).
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
mkdir -p results2
CASES="dfvlr_r4_M0.750_a1.25 dfvlr_r4_M0.775_a1.25 dfvlr_r4_M0.800_a1.25 dfvlr_r4_M0.825_a1.25 dfvlr_r4_M0.900_a1.25"
for c in $CASES; do
  if [ -f "results2/$c.done2" ]; then echo "skip2 $c" >> batchv3b.log; continue; fi
  if [ -f "results2/$c.done" ] && grep -q "^done" "results2/$c.done" 2>/dev/null; then echo "skipdone $c" >> batchv3b.log; continue; fi
  if [ ! -d "/meshCache/dfvlr_r4" ]; then echo "NOCACHE $c" >> batchv3b.log; continue; fi
  echo "$(date +%H:%M:%S) start $c" >> batchv3b.log
  rm -rf "$c/0.0"* "$c/0.1"* "$c/0.2"* "$c/postProcessing" "$c/constant/polyMesh" "$c/0/polyMesh" 2>/dev/null
  cp -r /meshCache/dfvlr_r4 "$c/constant/polyMesh"
  cd "$c"
  rhoCentralFoam > log.run 2>&1
  cd ..
  f=$(find "$c/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    tail -n 1 "$f" > "results2/$c.dat"
    echo done > "results2/$c.done"
    echo done > "results2/$c.done2"
    echo "$(date +%H:%M:%S) DONE $c $(tail -1 $f | cut -f1,5,7)" >> batchv3b.log
  else
    touch "results2/$c.dat"; echo no-forces > "results2/$c.done2"
    echo "$(date +%H:%M:%S) NOFORCES $c" >> batchv3b.log
  fi
done
echo "V3B WORKER FINISHED $(date +%H:%M:%S)" >> batchv3b.log
