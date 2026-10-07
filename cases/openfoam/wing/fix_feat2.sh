cd /case && perl -0pi -e 's/(nCellsBetweenLevels 3;\n)/$1\n    features\n    (\n    );\n/' system/snappyHexMeshDict && sed -n '80,92p' system/snappyHexMeshDict
