#!/bin/bash
# Fix time controls (deltaT was 1 -> zero-step loop) and launch all cases.
source /usr/lib/openfoam/openfoam2406/etc/bashrc
FT=$HOME/cfd/cases/openfoam/transonic/freetip
CASES="dfvlr_r4_M0.850_a1.25 dfvlr_r4_M0.925_a1.25 dfvlr_r4_M1.000_a1.25 rae2822_M0.850_a1.25 rae2822_M0.925_a1.25 rae2822_M1.000_a1.25 naca0012_M0.800_a1.25_ny30 naca0012_M0.800_a1.25_ny60"

for d in $CASES; do
  cd "$FT/$d" 2>/dev/null || { echo "MISSING $d"; continue; }
  cd "$FT/$d"
  foamDictionary -entry startFrom        -set startTime system/controlDict >/dev/null 2>&1
  foamDictionary -entry startTime        -set 0         system/controlDict >/dev/null 2>&1
  foamDictionary -entry deltaT           -set 1e-05     system/controlDict >/dev/null 2>&1
  foamDictionary -entry writeControl     -set adjustable system/controlDict >/dev/null 2>&1
  foamDictionary -entry writeInterval    -set 0.02      system/controlDict >/dev/null 2>&1
  foamDictionary -entry adjustTimeStep   -set yes       system/controlDict >/dev/null 2>&1
  foamDictionary -entry maxCo            -set 0.5       system/controlDict >/dev/null 2>&1
  foamDictionary -entry maxDeltaT        -set 1e-04     system/controlDict >/dev/null 2>&1
  foamDictionary -entry endTime          -set 0.12      system/controlDict >/dev/null 2>&1
  echo "$d: $(grep -E '^(deltaT|endTime|adjustTimeStep|maxCo)' system/controlDict | tr -s ' ' | tr '\n' ' ')"
done

echo "=== launching ==="
for d in $CASES; do
  cd "$FT/$d" 2>/dev/null || continue
  rm -f log.gcp
  setsid nohup rhoCentralFoam > log.gcp 2>&1 < /dev/null &
  printf 'launched %-34s pid=%s\n' "$d" "$!"
done
sleep 25
echo "running: $(pgrep -c rhoCentralFoam 2>/dev/null || echo 0)"
for d in $CASES; do
  printf '%-34s %s\n' "$d" "$(grep -E '^Time = ' $FT/$d/log.gcp 2>/dev/null | tail -1)"
done
