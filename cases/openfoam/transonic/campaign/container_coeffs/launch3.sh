#!/bin/bash
# launch 3 solver workers cleanly (no outer quoting traps)
cd /caseCampaign/campaign
for w in 0 1 2; do
  setsid nohup bash /tmp/solve_worker.sh "$w" 3 > "/caseCampaign/w$w.log" 2>&1 &
done
sleep 2
echo "workers: $(pgrep -fc solve_worker || true)"
tail -3 /caseCampaign/campaign/batch.log 2>/dev/null