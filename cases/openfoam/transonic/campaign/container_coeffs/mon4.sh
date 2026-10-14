echo "=== container time ==="
date -u
echo "=== running case progress (Time steps) ==="
for c in biconvex_M0.850_a1.25 biconvex_M0.875_a1.25 biconvex_M0.900_a1.25 biconvex_M0.925_a1.25; do
  t=$(grep '^Time = ' /caseCampaign/campaign/$c/log.run 2>/dev/null | tail -1)
  echo "  $c: $t"
done
echo "=== any forces yet ==="
for c in biconvex_M0.850_a1.25 biconvex_M0.875_a1.25 biconvex_M0.900_a1.25 biconvex_M0.925_a1.25; do
  f=$(find /caseCampaign/campaign/$c/postProcessing -name coefficient.dat 2>/dev/null | head -1)
  [ -n "$f" ] && echo "  $c: force data present" || echo "  $c: no forces yet"
done