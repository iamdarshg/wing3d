cd /caseCampaign/campaign
for g in rae2822 biconvex dfvlr_r4; do
  echo -n "$g STL: "
  ls -la ${g}_M0.900_a1.25/constant/triSurface/${g}.stl 2>/dev/null || echo MISSING
  # deploy to ALL cases of that geom (cases each carry their own copy)
done
echo "--- all case triSurface dirs ---"
ls */constant/triSurface/ 2>/dev/null | head -5