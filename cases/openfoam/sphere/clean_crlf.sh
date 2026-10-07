cd /case
sed -i 's/\r$//' system/snappyHexMeshDict system/fvSchemes system/fvSolution system/controlDict system/blockMeshDict system/decomposeParDict system/forceCoeffs system/meshQualityDict
sed -i 's/\r$//' 0/U 0/p constant/transportProperties constant/turbulenceProperties
echo CRLF-STRIPPED
