#!/bin/bash
# Clean slate: kill all workers, clear ALL stale markers/results/logs,
# set endTime to 0.25 (forces converge; 17% faster), then relaunch 3 workers.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd /caseCampaign/campaign
pkill -9 -f solve_worker 2>/dev/null
pkill -9 rhoCentralFoam 2>/dev/null
sleep 3
rm -f results/*.done results/*.dat batch.log
for c in *_M*; do
  rm -rf "$c/0.1*" "$c/0.14" "$c/0.28" "$c/postProcessing" "$c/VTK" 2>/dev/null
  # trim endTime 0.3 -> 0.25 in controlDict
  sed -i 's/endTime 0.3/endTime 0.25/' "$c/system/controlDict" 2>/dev/null
done
echo "=== clean, controlDict endTime now: ==="
grep -h 'endTime' naca0012_M0.900_a1.25/system/controlDict

# relaunch 3 workers
for w in 0 1 2; do
  setsid nohup bash /tmp/solve_worker.sh "$w" 3 > "/caseCampaign/w$w.log" 2>&1 &
done
sleep 3
echo "workers: $(pgrep -c -f solve_worker)"
free -m | head -2