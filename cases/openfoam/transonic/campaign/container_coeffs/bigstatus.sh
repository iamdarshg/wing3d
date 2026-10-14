#!/bin/bash
echo "=== REAL DATA SUMMARY (coefficient.dat directly, t>0.1) ==="
echo "geo            M      t      CL        CD"
for d in /caseCampaign/campaign/{naca0012,rae2822,biconvex,dfvlr_r4}_M*; do
  [ -d "$d" ] || continue
  f=$(find "$d/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    last=$(tail -1 "$f"); set -- $last
    t=$1; cl=$5; cd=$7
    # only print if t > 0.1 (converged) or the dfvlr/rae tails
    # print everything; mark converged
    geo=$(basename $d | sed 's/_M.*//')
    M=$(basename $d | sed 's/_.*M\([0-9.]*\)_.*/\1/')
    flag=""
    if [ $(echo "$t > 0.2" | bc) -eq 1 ]; then flag="CONV"; fi
    echo "$geo  $M  $t  $cl  $cd  $flag"
  fi
done