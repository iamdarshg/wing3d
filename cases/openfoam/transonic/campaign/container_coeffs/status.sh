echo "===== OF CAMPAIGN ====="
echo "done markers: $(ls /caseCampaign/campaign/results/*.done 2>/dev/null | wc -l)"
echo "rhoCentral running: $(pgrep -c rhoCentralFoam 2>/dev/null || echo 0)"
echo "--- batch.log tail ---"
tail -8 /caseCampaign/campaign/batch.log 2>/dev/null
echo ""
echo "--- completed cases (t=0.25, real data) ---"
for c in $(ls -d /caseCampaign/campaign/*_M* 2>/dev/null); do
  f=$(find $c/postProcessing -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    last=$(tail -1 "$f"); set -- $last
    echo "$(basename $c): t=$1 CL=$5 CD=$7"
  fi
done | sort
echo ""
echo "----- TSD files on host -----"
ls -la /caseCampaign/campaign/tsd_*.json 2>/dev/null || echo "none in container"