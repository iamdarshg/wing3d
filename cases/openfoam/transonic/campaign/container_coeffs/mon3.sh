echo "=== done so far ==="
ls /caseCampaign/campaign/results/*.done 2>/dev/null | wc -l
echo "=== batch tail ==="
tail -5 /caseCampaign/campaign/batch.log 2>/dev/null
echo "=== running ==="
pgrep -c rhoCentralFoam || echo 0
echo "=== killed total ==="
for w in 0 1 2 3; do grep -c Killed /caseCampaign/w$w.log 2>/dev/null; done | awk '{s+=$1} END{print s+0}'
echo "=== mem ==="
free -m | head -2