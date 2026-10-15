cd /caseCampaign/campaign
for c in dfvlr_r4_M0.850_a1.25 dfvlr_r4_M0.875_a1.25 dfvlr_r4_M0.900_a1.25 dfvlr_r4_M0.925_a1.25 dfvlr_r4_M0.950_a1.25 dfvlr_r4_M0.975_a1.25; do
  f=$(find $c/postProcessing -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    echo "$c"
    tail -1 "$f"
  else
    echo "$c: NO-DATA"
  fi
done
echo "--- Time steps: $(grep -c 'Time = ' dfvlr_r4_M0.975_a1.25/log.run 2>/dev/null)"