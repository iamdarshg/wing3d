#!/bin/bash
cd /caseCampaign/campaign
echo "== naca0.775 controlDict (running) =="
cat naca0012_M0.775_a1.25/system/controlDict
echo ""
echo "== grep attempts on pending case =="
grep -n "endTime" naca0012_M0.800_a1.25/system/controlDict
grep -n "maxCo\|maxDeltaT\|deltaT" naca0012_M0.800_a1.25/system/controlDict
echo ""
echo "== endTime across all 16 pending (with tr) =="
for c in rae2822_M0.775_a1.25 rae2822_M0.800_a1.25 rae2822_M0.825_a1.25 \
         naca0012_M0.750_a1.25 naca0012_M0.775_a1.25 naca0012_M0.800_a1.25 naca0012_M0.825_a1.25 \
         biconvex_M0.750_a1.25 biconvex_M0.775_a1.25 biconvex_M0.800_a1.25 biconvex_M0.825_a1.25 \
         dfvlr_r4_M0.750_a1.25 dfvlr_r4_M0.775_a1.25 dfvlr_r4_M0.800_a1.25 dfvlr_r4_M0.825_a1.25 dfvlr_r4_M0.900_a1.25; do
  printf "%s  " $c
  grep -H "endTime" $c/system/controlDict | tr -s ' ' | head -1
done
