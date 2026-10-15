echo "=== M0.900 progress ==="
grep 'Time = ' /caseCampaign/campaign/dfvlr_r4_M0.900_a1.25/log.run | tail -1
echo "=== M0.975 progress (if running) ==="
grep 'Time = ' /caseCampaign/campaign/dfvlr_r4_M0.975_a1.25/log.run | tail -1
echo "=== dfvlr cases with converged data ==="
cd /caseCampaign/campaign
for c in dfvlr_r4_M*; do
  f=$(find $c/postProcessing -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    echo "$c : $(tail -1 $f)"
  fi
done