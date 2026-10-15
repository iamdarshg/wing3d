#!/bin/bash
# Worker: re-run campaign cases on rebuilt caches. Priority = sanity cases first.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
W=$1; NP=$2
mkdir -p results2
PRS="rae2822_M0.850_a1.25 dfvlr_r4_M0.850_a1.25 biconvex_M0.850_a1.25"
LIST="$PRS"
for M in 0.750 0.775 0.800 0.825 0.875 0.900 0.925 0.950 0.975 1.000 1.025 1.050; do
  LIST="$LIST rae2822_M${M}_a1.25"
done
for M in 0.750 0.775 0.800 0.825 0.875 0.900 0.925 0.950 0.975 1.000 1.025 1.050 1.075 1.100; do
  LIST="$LIST dfvlr_r4_M${M}_a1.25"
done
for M in 0.750 0.775 0.800 0.825; do
  LIST="$LIST biconvex_M${M}_a1.25 naca0012_M${M}_a1.25"
done
i=0
for c in $LIST; do
  if [ $((i % NP)) -ne $W ]; then i=$((i+1)); continue; fi
  i=$((i+1))
  if [ -f "results2/$c.done" ]; then continue; fi
  if [ ! -d "$c" ]; then echo "MISSDIR $c" >> batch2.log; continue; fi
  g=$(echo "$c" | sed 's/_M.*//')
  if [ ! -d "/meshCache/$g" ]; then echo "NOCACHE $c" >> batch2.log; continue; fi
  echo "$(date +%H:%M:%S) W$W START $c" >> batch2.log
  rm -rf "$c/0.1"* "$c/0.2"* "$c/postProcessing" "$c/constant/polyMesh" "$c/0/polyMesh" 2>/dev/null
  cp -r "/meshCache/$g" "$c/constant/polyMesh"
  cd "$c"
  rhoCentralFoam > log.run 2>&1
  cd ..
  f=$(find "$c/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    tail -n 1 "$f" > "results2/$c.dat"
    echo done > "results2/$c.done"
    echo "$(date +%H:%M:%S) W$W DONE $c $(cut -f1,5,7 $f | tail -1)" >> batch2.log
  else
    touch "results2/$c.dat"; echo no-forces > "results2/$c.done"
    echo "$(date +%H:%M:%S) W$W NOFORCES $c" >> batch2.log
  fi
done
echo "W$W FINISHED $(date +%H:%M:%S)" >> batch2.log
