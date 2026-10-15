for g in naca0012 rae2822; do
  echo "=== $g M0.85 0/U ==="
  grep -a 'internalField' /caseCampaign/campaign/${g}_M0.850_a1.25/0/U 2>/dev/null
  grep -a -A3 'magUInf' /caseCampaign/campaign/${g}_M0.850_a1.25/system/controlDict 2>/dev/null | head -3
  grep -a 'Aref' /caseCampaign/campaign/${g}_M0.850_a1.25/system/controlDict 2>/dev/null
done