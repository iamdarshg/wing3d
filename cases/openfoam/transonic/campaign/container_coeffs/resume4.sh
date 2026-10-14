#!/bin/bash
cd /caseCampaign/campaign
pkill -9 -f solve_worker 2>/dev/null
pkill -9 rhoCentralFoam 2>/dev/null
sleep 2
rm -f results/*.done results/*.dat batch.log
for w in 0 1 2 3; do
  setsid nohup bash /tmp/solve_worker.sh "$w" 4 > "/caseCampaign/w$w.log" 2>&1 &
done
sleep 3
echo "workers=$(pgrep -c -f solve_worker)"
free -m | head -2