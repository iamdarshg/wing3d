source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign/rae2822_M0.850_a1.25
grep -a -A2 -E '^    [a-zA-Z]' constant/polyMesh/boundary | grep -a -E 'nFaces|patch|type|startFace' | head -20
echo "=== patch names ==="
grep -a -E '^    [a-zA-Z_][a-zA-Z0-9_. -]*$' constant/polyMesh/boundary | head
echo "=== naca for comparison ==="
cd ../naca0012_M0.850_a1.25
grep -a -E '^    [a-zA-Z_][a-zA-Z0-9_. -]*$' constant/polyMesh/boundary | head