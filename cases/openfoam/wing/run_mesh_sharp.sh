source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /case
blockMesh > log.blockMeshS 2>&1
snappyHexMesh -overwrite > log.snappyS 2>&1
echo DONESHARP >> log.snappyS
