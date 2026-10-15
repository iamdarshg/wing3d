#!/bin/bash
# Restart dfvlr M0.900 with the proven working numerics (maxCo 0.3,
# deltaT 1e-5) now that thermo corruption is fixed.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign/dfvlr_r4_M0.900_a1.25
# kill the over-slowed solver (PID 710 = fix900 run)
kill -9 710 2>/dev/null
pkill -9 -f 'dfvlr_r4_M0.900_a1.25' 2>/dev/null
sleep 2
# restore proven numerics
sed -i 's/maxCo 0.2/maxCo 0.3/; s/deltaT 5e-6/deltaT 1e-5/; s/maxDeltaT 2e-5/maxDeltaT 5e-5/' system/controlDict
grep -a -E 'maxCo|deltaT|maxDeltaT' system/controlDict
rm -rf 0.1* 0.14 0.28 postProcessing 2>/dev/null
rm -f ../results/dfvlr_r4_M0.900_a1.25.done ../results/dfvlr_r4_M0.900_a1.25.dat
setsid nohup bash -c 'source /usr/lib/openfoam/openfoam2406/etc/bashrc; cd '"$PWD"'; rhoCentralFoam > log.run 2>&1; f=$(find postProcessing -name coefficient.dat 2>/dev/null | head -1); if [ -n "$f" ]; then tail -1 "$f" > ../results/dfvlr_r4_M0.900_a1.25.dat; echo "M0.900 NORMAL done" >> ../batch.log; echo done > ../results/dfvlr_r4_M0.900_a1.25.done; fi' > /caseCampaign/m900v2.log 2>&1 &
echo "M0.900 restarted with standard numerics"