#!/bin/bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
c=rae2822_M0.850_a1.25
rm -rf $c/0.1* $c/0.14 $c/0.28 $c/postProcessing $c/0/polyMesh $c/constant/polyMesh 2>/dev/null
cp -r /meshCache/rae2822 $c/constant/polyMesh
cd $c
rhoCentralFoam > log.run 2>&1
f=$(find postProcessing -name coefficient.dat 2>/dev/null | head -1)
echo "rae2822 M0.85: $(tail -1 $f | cut -d' ' -f1,5,7)"