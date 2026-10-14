#!/bin/bash
# Retry dfvlr M0.900 with gentler startup numerics (maxCo 0.15).
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign/dfvlr_r4_M0.900_a1.25
rm -f ../results/dfvlr_r4_M0.900_a1.25.done ../results/dfvlr_r4_M0.900_a1.25.dat
rm -rf 0.1* 0.14 0.28 postProcessing 2>/dev/null
sed -i 's/maxCo 0.3/maxCo 0.15/; s/deltaT 1e-5/deltaT 5e-6/; s/maxDeltaT 5e-5/maxDeltaT 2e-5/' system/controlDict
rhoCentralFoam > log.run 2>&1
f=$(find postProcessing -name coefficient.dat 2>/dev/null | head -1)
if [ -n "$f" ]; then
  last=$(tail -1 "$f"); set -- $last
  echo "M0.900 RETRY done t=$1 CL=$5 CD=$7" >> ../batch.log
  echo "$last" > ../results/dfvlr_r4_M0.900_a1.25.dat
  echo done > ../results/dfvlr_r4_M0.900_a1.25.done
fi