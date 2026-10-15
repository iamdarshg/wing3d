#!/bin/bash
# Launch parallel solvers for the queued dfvlr cases (4 cores available).
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
for c in dfvlr_r4_M1.000_a1.25 dfvlr_r4_M1.025_a1.25 dfvlr_r4_M1.050_a1.25; do
  if [ -f "results/$c.done" ]; then echo "$c already done"; continue; fi
  # ensure standard numerics
  sed -i 's/maxCo 0.15/maxCo 0.3/; s/maxCo 0.2/maxCo 0.3/; s/deltaT 5e-6/deltaT 1e-5/; s/maxDeltaT 2e-5/maxDeltaT 5e-5/' "$c/system/controlDict"
  rm -rf "$c/0.1*" "$c/0.14" "$c/0.28" "$c/postProcessing" 2>/dev/null
  setsid nohup bash -c "source /usr/lib/openfoam/openfoam2406/etc/bashrc; cd /caseCampaign/campaign/$c; rhoCentralFoam > log.run 2>&1; f=\$(find postProcessing -name coefficient.dat 2>/dev/null | head -1); if [ -n \"\$f\" ]; then tail -1 \"\$f\" > ../results/$c.dat; echo done > ../results/$c.done; fi" > /caseCampaign/par_$c.log 2>&1 &
  echo "launched $c"
done
sleep 3
echo "running: $(pgrep -c rhoCentralFoam)"