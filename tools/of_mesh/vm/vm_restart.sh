#!/bin/bash
# Restore 0/ from git (a previous rm -rf [0-9]*/ deleted it), re-apply the
# rhoCentralFoam scheme + gentle time controls, relaunch.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd "$HOME/cfd" || exit 1
git checkout -- cases/openfoam/transonic/freetip/ 2>/dev/null
echo "restored 0/ dirs: $(ls cases/openfoam/transonic/freetip/*/0/U 2>/dev/null | wc -l)"

FT=$HOME/cfd/cases/openfoam/transonic/freetip
CASES="dfvlr_r4_M0.850_a1.25 dfvlr_r4_M0.925_a1.25 dfvlr_r4_M1.000_a1.25 rae2822_M0.850_a1.25 rae2822_M0.925_a1.25 rae2822_M1.000_a1.25 naca0012_M0.800_a1.25_ny30 naca0012_M0.800_a1.25_ny60"

for d in $CASES; do
  cd "$FT/$d" || continue
  S=system/fvSchemes
  sed -i 's/ddtSchemes *{ *default *steadyState; *}/ddtSchemes { default localEuler; }/' $S
  grep -q 'fluxScheme' $S || echo 'fluxScheme     Kurganov;' >> $S
  C=system/controlDict
  foamDictionary -entry deltaT    -set 1e-07 $C >/dev/null 2>&1
  foamDictionary -entry maxCo     -set 0.1   $C >/dev/null 2>&1
  foamDictionary -entry maxDeltaT -set 5e-05 $C >/dev/null 2>&1
  foamDictionary -entry endTime   -set 0.12  $C >/dev/null 2>&1
  foamDictionary -entry adjustTimeStep -set yes $C >/dev/null 2>&1
  rm -f log.gcp
done

echo "=== launching ==="
for d in $CASES; do
  cd "$FT/$d" || continue
  setsid nohup rhoCentralFoam > log.gcp 2>&1 < /dev/null &
done
sleep 60
echo "running: $(pgrep -c rhoCentralFoam 2>/dev/null || echo 0)"
for d in $CASES; do
  t=$(grep -E '^Time = ' "$FT/$d/log.gcp" 2>/dev/null | tail -1 | awk '{print $3}')
  f=$(grep -c 'FOAM FATAL' "$FT/$d/log.gcp" 2>/dev/null)
  printf '%-34s t=%-12s fatal=%s\n' "$d" "${t:-none}" "$f"
done
