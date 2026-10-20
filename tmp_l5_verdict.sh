#!/bin/bash
# L5 verdict calculator: coarse-lane final CL vs refined-mesh (86k) anchors.
# Reads the last line of each lane's postProcessing/forces/0/coefficient.dat
# (12-column; CL is column 5), compares vs the bit-exact 86k probe values
# (l5_probe_dfvlr_M<M>_86k.dat), and emits the verdict per the decision rule:
#   <1%  -> refined-vs-coarse dCL flat -> OF saturation genuine
#   1-3% -> marginally resolved
#   >=3% -> resolution matters (ceiling story merges with L2)
# Usage: docker exec oftrans bash /probes/l5_verdict.sh
set -u
L1=/caseCampaign/probes/dfvlr_r4_M0.925_coarse45k
L2=/caseCampaign/probes/dfvlr_r4_M1.000_coarse45k
A1=0.2779525   # refined-mesh anchor for M0.925 (probe2, bit-exact vs corpus)
A2=0.2126116   # refined-mesh anchor for M1.000
printf '%-8s %-12s %-12s %8s %s\n' M 'CL_coarse' 'CL_refined' dCL verdict
for spec in "$L1 $A1 0.925" "$L2 $A2 1.000"; do
  set -- $spec
  cfile=$(find "$1/postProcessing/forces" -name coefficient.dat -printf '%T@ %p\n' | sort -n | tail -1 | cut -d' ' -f2-)
  line=$(tail -n 1 "$cfile")
  cl=$(echo "$line" | awk '{print $5}')
  d=$(awk -v c="$cl" -v a="$2" 'BEGIN{printf "%.2f", (c-a)/a*100}')
  ad=$(awk -v v="$d" 'BEGIN{if (v<0) v=-v; print v}')
  if awk -v v="$ad" 'BEGIN{exit !(v<1)}'; then
    verdict="flat: saturation genuine"
  elif awk -v v="$ad" 'BEGIN{exit !(v>=3)}'; then
    verdict="resolution matters (>=3%)"
  else
    verdict="marginally resolved (1-3%)"
  fi
  printf '%-8s %-12s %-12s %+7s%% %s\n' "$3" "$cl" "$2" "$d" "$verdict"
done