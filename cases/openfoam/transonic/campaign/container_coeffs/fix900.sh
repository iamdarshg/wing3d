#!/bin/bash
# Repair dfvlr M0.900: known-good thermo/turbulence, clean snappy residue,
# gentle numerics, fresh run.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
c=dfvlr_r4_M0.900_a1.25
# borrow known-good constant files from a converged dfvlr case
for f in thermophysicalProperties turbulenceProperties; do
  cp dfvlr_r4_M0.850_a1.25/constant/$f $c/constant/$f
done
# borrow 0/ fields too, then fix the patch name
for f in U p T; do
  sed 's/"dfvlr_r4.stl"/"dfvlr_r4.stl"/' dfvlr_r4_M0.850_a1.25/0/$f > $c/0/$f
done
# verify 0/U has the right patch name (dfvlr_r4.stl) and check vel
grep -a 'dfvlr_r4' $c/0/U | head -2
# remove snappy residue inside 0/
rm -rf $c/0/polyMesh $c/0/cellLevel $c/0/pointLevel
# gentle numerics
sed -i 's/maxCo 0.5/maxCo 0.2/; s/deltaT 2e-5/deltaT 5e-6/; s/maxDeltaT 1e-4/maxDeltaT 2e-5/' $c/system/controlDict
echo "--- fixed controlDict:"
grep -a -E 'maxCo|deltaT' $c/system/controlDict
cd $c
rm -rf 0.1* 0.14 0.28 postProcessing 2>/dev/null
setsid nohup bash -c 'source /usr/lib/openfoam/openfoam2406/etc/bashrc; cd '"$c"'; rhoCentralFoam > log.run 2>&1; f=$(find postProcessing -name coefficient.dat 2>/dev/null | head -1); if [ -n "$f" ]; then tail -1 "$f" > ../results/dfvlr_r4_M0.900_a1.25.dat; echo "M0.900 FIXED t=$(tail -1 $f | cut -d\" \" -f1)" >> ../batch.log; echo done > ../results/dfvlr_r4_M0.900_a1.25.done; fi' > /caseCampaign/m900run. 2>&1 &
echo "M0.900 RERUN LAUNCHED"