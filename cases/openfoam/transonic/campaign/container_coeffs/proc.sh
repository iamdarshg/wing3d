echo "=== processes ==="
pgrep -af 'rhoCentralFoam|retry900|dfvlr_runner' | head -8
echo ""
echo "=== retry900 log ==="
cat /caseCampaign/retry900.log 2>/dev/null
echo ""
echo "=== M0.900 log.run tail ==="
tail -6 /caseCampaign/campaign/dfvlr_r4_M0.900_a1.25/log.run 2>/dev/null