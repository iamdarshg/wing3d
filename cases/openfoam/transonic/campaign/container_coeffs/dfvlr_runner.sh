#!/bin/bash
# Dedicated DFVLR runner: runs all 11 dfvlr cases sequentially with the
# new 86k mesh + hardened numerics. Doesn't touch other geometries.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
mkdir -p results
for c in dfvlr_r4_M*; do
  if [ -f "results/$c.done" ] && [ -s "results/$c.dat" ]; then
    echo "skip $c (done)"; continue
  fi
  echo "$(date +%H:%M:%S) DFVLRSOLVE START $c" >> batch.log
  rm -rf "$c/0.1*" "$c/0.14" "$c/0.28" "$c/postProcessing" 2>/dev/null
  cd "$c"
  rhoCentralFoam > log.run 2>&1
  cd ..
  f=$(find "$c/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ]; then
    last=$(tail -n 1 "$f"); set -- $last
    echo "$c done t=$1 CL=$5 CD=$7" >> batch.log
    echo "$last" > "results/$c.dat"
    echo "done" > "results/$c.done"
  else
    echo "$c NOFORCES" >> batch.log
    touch "results/$c.dat"; echo "no-forces" > "results/$c.done"
  fi
done
echo "DFVLR-RUNNER DONE $(date +%H:%M:%S)" >> batch.log