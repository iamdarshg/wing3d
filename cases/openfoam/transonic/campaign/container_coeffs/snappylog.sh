cd /caseCampaign/campaign
for g in rae2822 biconvex; do
  echo "=== $g snappy tail ==="
  tail -12 ${g}_M0.900_a1.25/log.snappy 2>/dev/null
  echo "=== $g snappy cells line ==="
  grep -a -i -m3 'Finished meshing\|cells\b\|number of cells\|Castellated\|terminating' ${g}_M0.900_a1.25/log.snappy 2>/dev/null | head -5
done