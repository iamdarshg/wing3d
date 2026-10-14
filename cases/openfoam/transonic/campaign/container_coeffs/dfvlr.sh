c=dfvlr_r4_M1.000_a1.25
echo "=== full FATAL ERROR ==="
grep -a -A15 'FOAM FATAL' /caseCampaign/campaign/$c/log.run 2>/dev/null | head -25
echo ""
echo "=== near the end ==="
tail -20 /caseCampaign/campaign/$c/log.run 2>/dev/null
echo ""
echo "=== checkMesh stats ==="
grep -a -E 'cells:|Mesh OK|non-orthogonality|skewness' /caseCampaign/campaign/$c/log.check 2>/dev/null | head -6