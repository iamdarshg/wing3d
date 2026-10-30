#!/bin/bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
FT=$HOME/cfd/cases/openfoam/transonic/freetip
cd "$FT" || exit 1
echo "running: $(pgrep -c rhoCentralFoam 2>/dev/null || echo 0)"
printf '%-34s %-12s %-8s %s\n' CASE TIME FATAL CL
for d in */; do
  d=${d%/}
  [ -f "$d/log.gcp" ] || continue
  t=$(grep -E '^Time = ' "$d/log.gcp" 2>/dev/null | tail -1 | awk '{print $3}')
  f=$(grep -cE 'FOAM FATAL|Floating point|Aborted' "$d/log.gcp" 2>/dev/null)
  cl=$(grep -E '^\s+Cl:' "$d/log.gcp" 2>/dev/null | tail -1 | awk '{print $2}')
  printf '%-34s %-12s %-8s %s\n' "$d" "${t:-none}" "$f" "${cl:-}"
done
echo "--- any FATAL text ---"
grep -h -A2 -E 'FOAM FATAL' */log.gcp 2>/dev/null | head -8
uptime
