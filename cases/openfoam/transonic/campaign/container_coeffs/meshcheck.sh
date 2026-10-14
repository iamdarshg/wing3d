echo "=== dfvlr M0.85 current polyMesh points/owner ==="
ls -la /caseCampaign/campaign/dfvlr_r4_M0.850_a1.25/constant/polyMesh/ 2>/dev/null | grep -E 'owner|points|faces' | head
echo "=== cells from owner file size (4 bytes/cell + header) ==="
sz=$(stat -c%s /caseCampaign/campaign/dfvlr_r4_M0.850_a1.25/constant/polyMesh/owner 2>/dev/null)
echo "owner bytes: $sz  (~$(( (sz-300)/4 )) cells)"
echo "=== is the running solver writing to 0.1* dirs? ==="
ls -d /caseCampaign/campaign/dfvlr_r4_M0.850_a1.25/0.1* 2>/dev/null | head -3
echo "=== compare: the scratch fixed mesh ==="
ls -la /dfvlr_mesh/constant/polyMesh/owner 2>/dev/null