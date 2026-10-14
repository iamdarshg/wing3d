echo "=== M0.900 FATAL message ==="
grep -a -B4 -A4 'FOAM FATAL\|-->\|^[[:space:]]*[A-Za-z].*error' /caseCampaign/campaign/dfvlr_r4_M0.900_a1.25/log.run 2>/dev/null | head -30
echo ""
echo "=== does the case have 0/ files? ==="
ls /caseCampaign/campaign/dfvlr_r4_M0.900_a1.25/0/ 2>/dev/null
echo "=== controlDict patch name ==="
grep -a 'endTime' /caseCampaign/campaign/dfvlr_r4_M0.900_a1.25/system/controlDict 2>/dev/null