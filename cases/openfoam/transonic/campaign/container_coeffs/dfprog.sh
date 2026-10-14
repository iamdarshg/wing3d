echo "=== dfvlr runner progress ==="
tail -6 /caseCampaign/dfvlr_runner.log 2>/dev/null
echo "=== dfvlr cases with data now ==="
for c in /caseCampaign/campaign/dfvlr_r4_M*; do
  f=$(find $c/postProcessing -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    last=$(tail -1 "$f"); set -- $last
    echo "$(basename $c): t=$1 CL=$5 CD=$7"
  fi
done
echo "=== rhoCentral running ==="
pgrep -c rhoCentralFoam || echo 0