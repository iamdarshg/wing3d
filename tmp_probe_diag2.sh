#!/bin/bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
echo "=A: dfvlr_coarse 0-dir ="
ls -la /dfvlr_coarse/0/
echo "=B: dfvlr_coarse 0/polyMesh ="
ls -la /dfvlr_coarse/0/polyMesh/
echo "=C: cmA mesh lines ="
grep -a -E "Exec   :|Create|mesh for time|points:|cells:" /tmp/cmA.log | head -12
echo "=D: checkMesh -constant on coarse ="
checkMesh -case /dfvlr_coarse -constant > /tmp/cmE.log 2>&1
tail -n 3 /tmp/cmE.log
grep -a -c "Failed" /tmp/cmE.log
echo "=E: combo x5 test (constant topology + 0/polyMesh points) ="
rm -rf /tmp/x5
mkdir -p /tmp/x5/constant/polyMesh /tmp/x5/system
cp /dfvlr_coarse/constant/polyMesh/faces /dfvlr_coarse/constant/polyMesh/owner /dfvlr_coarse/constant/polyMesh/neighbour /dfvlr_coarse/constant/polyMesh/boundary /tmp/x5/constant/polyMesh/
cp /dfvlr_coarse/0/polyMesh/points /tmp/x5/constant/polyMesh/points
cp /dfvlr_coarse/system/controlDict /tmp/x5/system/
checkMesh -case /tmp/x5 > /tmp/x5.log 2>&1
echo "combo-verdict:"
tail -n 3 /tmp/x5.log
echo "=F: good points vs stale points md5 ="
md5sum /dfvlr_coarse/0/polyMesh/points /dfvlr_coarse/constant/polyMesh/points
echo "=G: end="
