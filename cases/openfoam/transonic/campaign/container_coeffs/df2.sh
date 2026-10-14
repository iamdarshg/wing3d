echo "=== dfvlr runner log (full) ==="
cat /caseCampaign/dfvlr_runner2.log 2>/dev/null
echo "=== batch.log dfvlr tail ==="
grep dfvlr /caseCampaign/campaign/batch.log | tail -6
echo "=== dfvlr data ==="
cd /caseCampaign/campaign
for c in dfvlr_r4_M*; do
  f=$(find $c/postProcessing -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    last=$(tail -1 "$f"); set -- $last
    echo "$c: t=$1 CL=$5 CD=$7"
  else
    echo "$c: none"
  fi
done