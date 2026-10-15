#!/bin/bash
# Solver worker for the 3 CORRECTED geometries (rae2822, biconvex, dfvlr_r4).
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
W=$1; NP=$2
mkdir -p results
CASES=$(ls -d rae2822_M* biconvex_M* dfvlr_r4_M* 2>/dev/null | sort)
i=0
for c in $CASES; do
  if [ $((i % NP)) -ne $W ]; then i=$((i+1)); continue; fi
  i=$((i+1))
  # skip cases already converged (t>0.2 in coefficient.dat)
  f=$(find "$c/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    t=$(tail -1 "$f" | cut -d' ' -f1)
    ok=$(awk -v t="$t" 'BEGIN{print (t>0.2)?1:0}')
    if [ "$ok" = "1" ]; then
      echo "skip $c (converged t=$t)"
      tail -1 "$f" > "results/$c.dat"; echo done > "results/$c.done"
      continue
    fi
  fi
  g=$(echo "$c" | sed 's/_M.*//')
  if [ ! -d "/meshCache/$g" ]; then echo "NOCACHE $c"; continue; fi
  echo "$(date +%H:%M:%S) W$W START $c" >> batch.log
  rm -rf "$c/0.1*" "$c/0.14" "$c/0.28" "$c/postProcessing" "$c/0/polyMesh" 2>/dev/null
  cp -r "/meshCache/$g" "$c/constant/polyMesh"
  cd "$c"
  rhoCentralFoam > log.run 2>&1
  cd ..
  f=$(find "$c/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    tail -n 1 "$f" > "results/$c.dat"
    echo done > "results/$c.done"
    echo "$(date +%H:%M:%S) W$W DONE $c" >> batch.log
  else
    touch "results/$c.dat"; echo no-forces > "results/$c.done"
    echo "$(date +%H:%M:%S) W$W NOFORCES $c" >> batch.log
  fi
done
echo "W$W FINISHED $(date +%H:%M:%S)" >> batch.log