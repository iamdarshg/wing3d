#!/bin/bash
# Deploy corrected STLs to all case dirs + launch 4 corrected-geometry workers.
for g in rae2822 biconvex dfvlr_r4; do
  for d in /caseCampaign/campaign/${g}_M*; do
    cp "/caseCampaign/campaign/${g}_M0.900_a1.25/constant/triSurface/${g}.stl" \
       "$d/constant/triSurface/${g}.stl" 2>/dev/null
  done
done
echo "STLs synced"
for w in 0 1 2 3; do
  setsid nohup bash /tmp/solver2.sh "$w" 4 > "/caseCampaign/v2w$w.log" 2>&1 &
done
sleep 3
echo "workers: $(pgrep -c -f solver2.sh)"
echo "solver running: $(pgrep -c rhoCentralFoam)"