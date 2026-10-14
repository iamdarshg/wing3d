echo "done=$(ls /caseCampaign/campaign/results/*.done 2>/dev/null | wc -l)"
echo "=== batch tail ==="
tail -6 /caseCampaign/campaign/batch.log 2>/dev/null
echo "=== rhoCentral ==="
pgrep -c rhoCentralFoam || echo 0
echo "=== completed OF rows (t,CL,CD) ==="
for c in $(ls -d /caseCampaign/campaign/*_M* 2>/dev/null); do
  f=$(find $c/postProcessing -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    last=$(tail -1 "$f")
    set -- $last
    echo "$(basename $c): t=$1 CL=$5 CD=$7"
  fi
done | sort | head -12