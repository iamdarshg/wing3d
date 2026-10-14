#!/bin/bash
# Clean restart: clear all stale results, relaunch 2 workers (split NP=2).
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
pkill -9 -f solve_worker 2>/dev/null
pkill -9 rhoCentralFoam 2>/dev/null
sleep 2
rm -f results/*.done results/*.dat 2>/dev/null
rm -f batch.log
for w in 0 1; do
  setsid nohup bash /tmp/solve_worker.sh "$w" 2 > "/caseCampaign/w$w.log" 2>&1 &
done
sleep 3
echo "workers: $(pgrep -c -f solve_worker)"
echo "rhoCentral: $(pgrep -c rhoCentralFoam)"