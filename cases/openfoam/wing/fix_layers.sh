cd /case && sed -i 's/(lowerWall|wing)/wing/' system/snappyHexMeshDict && grep -n '"wing' system/snappyHexMeshDict | head -3
