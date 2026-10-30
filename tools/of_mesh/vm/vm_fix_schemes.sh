#!/bin/bash
# The committed freetip fvSchemes is the rhoSimpleFoam variant (steadyState ddt),
# which gives rho a ZERO diagonal in rhoCentralFoam -> FPE. Patch to the
# rhoCentralFoam (transient localEuler) set and relaunch.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
FT=$HOME/cfd/cases/openfoam/transonic/freetip
CASES="dfvlr_r4_M0.850_a1.25 dfvlr_r4_M0.925_a1.25 dfvlr_r4_M1.000_a1.25 rae2822_M0.850_a1.25 rae2822_M0.925_a1.25 rae2822_M1.000_a1.25 naca0012_M0.800_a1.25_ny30 naca0012_M0.800_a1.25_ny60"

for d in $CASES; do
  cd "$FT/$d" || continue
  S=system/fvSchemes
  sed -i 's/ddtSchemes *{ *default *steadyState; *}/ddtSchemes { default localEuler; }/' $S
  # ensure the central-upwind knobs rhoCentralFoam needs
  grep -q 'fluxScheme' $S || echo 'fluxScheme     Kurganov;' >> $S
  grep -q 'reconstruct(rho)' $S || echo 'reconstruct(rho) vanLeer;' >> $S
  grep -q 'reconstruct(U)' $S || echo 'reconstruct(U)   vanLeer;' >> $S
  grep -q 'reconstruct(T)' $S || echo 'reconstruct(T)   vanLeer;' >> $S
  # rhoCentralFoam has no SIMPLE loop; drop relaxation use by removing rhoSimple/SIMPLE blocks
  printf '%-34s %s\n' "$d" "$(grep -E 'ddtSchemes' $S | head -1)"
done

echo "=== launching ==="
for d in $CASES; do
  cd "$FT/$d" || continue
  rm -f log.gcp
  setsid nohup rhoCentralFoam > log.gcp 2>&1 < /dev/null &
  printf 'launched %-34s pid=%s\n' "$d" "$!"
done
sleep 30
echo "running: $(pgrep -c rhoCentralFoam 2>/dev/null || echo 0)"
for d in $CASES; do
  printf '%-34s %s\n' "$d" "$(grep -E '^Time = ' $FT/$d/log.gcp 2>/dev/null | tail -1)"
done
