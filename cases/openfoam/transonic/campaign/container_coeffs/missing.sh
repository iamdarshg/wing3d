cd /caseCampaign/campaign
for c in naca0012_M* rae2822_M* biconvex_M* dfvlr_r4_M*; do
  if [ ! -f "results/$c.done" ]; then echo "MISSING $c"; fi
done
echo "=== done count: $(ls results/*.done 2>/dev/null | wc -l) / 44 ==="