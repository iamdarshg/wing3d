#!/bin/bash
# w5 worker: steals the naca0012 low-M cases (claims them so the main
# queue skips them). Same proven flow as the other workers.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
mkdir -p results2
CASES="naca0012_M0.775_a1.25 naca0012_M0.800_a1.25 naca0012_M0.825_a1.25"
for c in $CASES; do
  if [ -f "results2/$c.done2" ]; then echo "skip2 $c" >> batchw5.log; continue; fi
  echo claimed-w5 > "results2/$c.done"
  if [ ! -d "/meshCache/naca0012" ]; then echo "NOCACHE $c" >> batchw5.log; continue; fi
  echo "$(date +%H:%M:%S) start $c" >> batchw5.log
  rm -rf "$c/0.0"* "$c/0.1"* "$c/0.2"* "$c/postProcessing" "$c/constant/polyMesh" "$c/0/polyMesh" 2>/dev/null
  cp -r /meshCache/naca0012 "$c/constant/polyMesh"
  cd "$c"
  rhoCentralFoam > log.run 2>&1
  cd ..
  f=$(find "$c/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    tail -n 1 "$f" > "results2/$c.dat"
    echo done > "results2/$c.done"
    echo done > "results2/$c.done2"
    echo "$(date +%H:%M:%S) DONE $c $(tail -1 $f | cut -f1,5,7)" >> batchw5.log
  else
    touch "results2/$c.dat"; echo no-forces > "results2/$c.done2"
    echo "$(date +%H:%M:%S) NOFORCES $c" >> batchw5.log
  fi
done
echo "W5 WORKER FINISHED $(date +%H:%M:%S)" >> batchw5.log
