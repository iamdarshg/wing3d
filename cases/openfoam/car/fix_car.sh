cd /case && sed -i 's/wing\r$/car\r/' system/snappyHexMeshDict && sed -i 's/inGroups (motorBikeGroup);//' system/snappyHexMeshDict && grep -n 'car' system/snappyHexMeshDict | head -6
