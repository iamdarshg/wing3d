#!/bin/bash
# Full corpus mesh map: nCells per (geo, M) case from its polyMesh note.
cd /caseCampaign/campaign
MS="0.750 0.775 0.800 0.825 0.850 0.875 0.900 0.925 0.950 0.975 1.000 1.025 1.050 1.075 1.100"
for g in naca0012 rae2822 biconvex dfvlr_r4; do
  echo "== $g"
  for M in $MS; do
    C="${g}_M${M}_a1.25"
    N=""
    if [ -f "$C/constant/polyMesh/owner" ]; then
      N=$(grep -m1 -oE 'nCells:[0-9]+' "$C/constant/polyMesh/owner")
    fi
    D="--"
    if [ -f "results2/$C.dat" ] && [ -s "results2/$C.dat" ]; then
      D=$(awk 'NR==1{printf "t=%s CL=%s", $1, $5}' "results2/$C.dat")
    elif [ -f "results2/$C.done2" ]; then
      D=$(cat "results2/$C.done2" 2>/dev/null | head -c 20)
    fi
    printf "  M%-6s %-14s %s\n" "$M" "${N:-no-mesh}" "$D"
  done
done
echo ""
echo "===== container_coeffs location ====="
ls -d /caseCampaign/*/ 2>/dev/null
find /caseCampaign -maxdepth 2 -name '*coeff*' -o -maxdepth 2 -name 'container_coeffs' 2>/dev/null | head
ls /caseCampaign/container_coeffs 2>/dev/null | head -8
echo ""
echo "===== old dfvlr 0.750/0.775 (if in container_coeffs / old dats) ====="
for f in /caseCampaign/container_coeffs/dfvlr_r4_M0.750_a1.25.dat /caseCampaign/container_coeffs/dfvlr_r4_M0.775_a1.25.dat; do
  if [ -f "$f" ]; then echo "--- $f"; cat "$f" | head -3; ls -la --time-style=long-iso "$f"; fi
done
echo ""
echo "===== results2 dfvlr 0.750/0.775/0.800/0.825/0.900 + mtimes ====="
ls -la --time-style=long-iso results2/dfvlr_r4_M0.7*.dat results2/dfvlr_r4_M0.8*.dat results2/dfvlr_r4_M0.9*.dat 2>/dev/null
for f in results2/dfvlr_r4_M0.750_a1.25.dat results2/dfvlr_r4_M0.775_a1.25.dat; do
  echo "--- $f"; head -1 "$f"
done
echo DONE
