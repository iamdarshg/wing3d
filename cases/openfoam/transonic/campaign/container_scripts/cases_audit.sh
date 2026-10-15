#!/bin/bash
cd /caseCampaign/campaign
echo "== clock/disk/mem =="
date
df -h / | tail -1
free -m | head -2
echo ""
echo "== pending cases audit =="
for c in rae2822_M0.775_a1.25 rae2822_M0.800_a1.25 rae2822_M0.825_a1.25 \
         naca0012_M0.750_a1.25 naca0012_M0.775_a1.25 naca0012_M0.800_a1.25 naca0012_M0.825_a1.25 \
         biconvex_M0.750_a1.25 biconvex_M0.775_a1.25 biconvex_M0.800_a1.25 biconvex_M0.825_a1.25 \
         dfvlr_r4_M0.750_a1.25 dfvlr_r4_M0.775_a1.25 dfvlr_r4_M0.800_a1.25 dfvlr_r4_M0.825_a1.25 dfvlr_r4_M0.900_a1.25; do
  [ -d $c ] || { echo "$c MISSING"; continue; }
  et=$(grep -E '^endTime' $c/system/controlDict | awk '{print $2}' | tr -d ';')
  mc=$(grep -E '^maxCo' $c/system/controlDict | awk '{print $2}' | tr -d ';')
  wc=$(grep -E '^writeControl' $c/system/controlDict | awk '{print $2}' | tr -d ';')
  wi=$(grep -E '^writeInterval' $c/system/controlDict | awk '{print $2}' | tr -d ';')
  nd=$(ls -d $c/[0-9]* 2>/dev/null | wc -l)
  sz=$(du -sh $c 2>/dev/null | awk '{print $1}')
  echo "$c end=$et co=$mc write=$wc/$wi timeDirs=$nd size=$sz"
done
echo ""
echo "== functions block (naca0.800) =="
sed -n '/^functions/,/^}/p' naca0012_M0.800_a1.25/system/controlDict
echo ""
echo "== running cases: latest time dirs =="
for c in rae2822_M0.775_a1.25 rae2822_M0.800_a1.25 naca0012_M0.775_a1.25 biconvex_M0.750_a1.25 dfvlr_r4_M0.750_a1.25; do
  echo "$c: $(ls -d $c/[0-9]* 2>/dev/null | xargs -n1 basename 2>/dev/null | sort -g | tail -4 | tr '\n' ' ')"
done
