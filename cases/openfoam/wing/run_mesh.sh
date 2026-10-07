source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /case
blockMesh > log.blockMesh 2>&1
snappyHexMesh -overwrite > log.snappy6 2>&1
echo DONE6 >> log.snappy6
