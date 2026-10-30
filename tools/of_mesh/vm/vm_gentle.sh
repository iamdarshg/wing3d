#!/bin/bash
# Gentle start for rhoCentralFoam: small maxCo/maxDeltaT to survive the
# freestream -> shocked-flow transient (it died on the janaf T inversion).
source /usr/lib/openfoam/openfoam2406/etc/bashrc
FT=$HOME/cfd/cases/openfoam/transonic/freetip
CASES="dfvlr_r4_M0.850_a1.25 dfvlr_r4_M0.925_a1.25 dfvlr_r4_M1.000_a1.25 rae2822_M0.850_a1.25 rae2822_M0.925_a1.25 rae2822_M1.000_a1.25 naca0012_M0.800_a1.25_ny30 naca0012_M0.800_a1.25_ny60"

pkill -f rhoCentralFoam 2>/dev/null
sleep 2

for d in $CASES; do
  cd "$FT/$d" || continue
  foamDictionary -entry deltaT    -set 1e-07 system/controlDict >/dev/null 2>&1
  foamDictionary -entry maxCo     -set 0.1   system/controlDict >/dev/null 2>&1
  foamDictionary -entry maxDeltaT -set 5e-05 system/controlDict >/dev/null 2>&1
  foamDictionary -entry endTime   -set 0.12  system/controlDict >/dev/null 2>&1
  rm -rf [0-9]*/ 0/uniform 2>/dev/null
  # keep 0/ only
  rm -f log.gcp
done

echo "=== launching (gentle) ==="
for d in $CASES; do
  cd "$FT/$d" || continue
  setsid nohup rhoCentralFoam > log.gcp 2>&1 < /dev/null &
  printf 'launched %-34s pid=%s\n' "$d" "$!"
done
sleep 45
echo "running: $(pgrep -c rhoCentralFoam 2>/dev/null || echo 0)"
for d in $CASES; do
  t=$(grep -E '^Time = ' "$FT/$d/log.gcp" 2>/dev/null | tail -1 | awk '{print $3}')
  printf '%-34s t=%s\n' "$d" "${t:-none}"
done
