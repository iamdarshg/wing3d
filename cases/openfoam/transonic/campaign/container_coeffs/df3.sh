echo "=== dfvlr data now ==="
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
echo ""
echo "=== retry900 (M0.900) ==="
c=dfvlr_r4_M0.900_a1.25
f=$(find $c/postProcessing -name coefficient.dat 2>/dev/null | head -1)
[ -n "$f" ] && echo "M0.900 retry: $(tail -1 $f | awk '{print $1, $5, $7}')" || echo "M0.900 retry: no data"
echo "=== batch tail ==="
tail -4 batch.log
echo "=== solvers ==="
pgrep -af rhoCentralFoam | wc -l