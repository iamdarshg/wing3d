#!/bin/bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
mkdir -p results2
C=dfvlr_r4_M0.900_a1.25
if [ -f "results2/$C.done2" ] && ! grep -q '^claimed' "results2/$C.done2"; then
  echo "skip $C (already done)" >> batchv4.log; exit 0; fi
if [ ! -d "/meshCache/dfvlr_r4" ]; then echo NOCACHE >> batchv4.log; exit 0; fi
echo "$(date +%H:%M:%S) start $C" >> batchv4.log
rm -rf "$C/0.0"* "$C/0.1"* "$C/0.2"* "$C/postProcessing" "$C/constant/polyMesh" "$C/0/polyMesh" 2>/dev/null
cp -r /meshCache/dfvlr_r4 "$C/constant/polyMesh"
cd "$C"
sed -i 's/^startFrom .*/startFrom startTime;/; s/^startTime .*/startTime 0;/' system/controlDict
sed -i 's/^endTime .*/endTime 0.25;/; s/^deltaT .*/deltaT 2e-5;/' system/controlDict
sed -i 's/^maxCo .*/maxCo 0.5;/; s/^maxDeltaT .*/maxDeltaT 1e-4;/' system/controlDict
rhoCentralFoam > log.run 2>&1
cd ..
f=$(find "$C/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
if [ -n "$f" ] && [ -s "$f" ]; then
  tail -n 1 "$f" > "results2/$C.dat"
  echo done > "results2/$C.done"
  echo done > "results2/$C.done2"
  echo "$(date +%H:%M:%S) DONE $C $(tail -1 $f | cut -f1,5,7)" >> batchv4.log
else
  touch "results2/$C.dat"; echo no-forces > "results2/$C.done2"
  echo "$(date +%H:%M:%S) NOFORCES $C" >> batchv4.log
fi
