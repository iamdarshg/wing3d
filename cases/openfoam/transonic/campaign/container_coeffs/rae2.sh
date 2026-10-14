echo "=== rae2822 .dat sizes ==="
for c in rae2822_M*; do
  sz=$(stat -c%s /caseCampaign/campaign/results/$c.dat 2>/dev/null || echo missing)
  echo "$c: dat=$sz done=$(test -f /caseCampaign/campaign/results/$c.done && echo Y || echo N)"
done
echo ""
echo "=== rae2822 with real force data (t=0.25) ==="
for c in /caseCampaign/campaign/rae2822_M*; do
  f=$(find $c/postProcessing -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    last=$(tail -1 "$f"); set -- $last
    echo "$(basename $c): t=$1 CL=$5 CD=$7"
  fi
done