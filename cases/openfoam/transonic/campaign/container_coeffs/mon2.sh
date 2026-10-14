sleep 240
echo "=== rhoCentral alive ==="
pgrep -c rhoCentralFoam || echo 0
echo "=== Killed events (new) ==="
for w in 0 1 2 3; do
  n=$(grep -c Killed /caseCampaign/w$w.log 2>/dev/null || echo 0)
  echo "  w$w: $n"
done
echo "=== batch tail ==="
tail -8 /caseCampaign/campaign/batch.log 2>/dev/null
echo "=== done ==="
ls /caseCampaign/campaign/results/*.done 2>/dev/null | wc -l