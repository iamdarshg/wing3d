#!/bin/bash
cd /caseCampaign/campaign
echo "== solvers:"
for p in $(pgrep -x rhoCentralFoam); do readlink /proc/$p/cwd; done
echo "== newest coefficients:"
for c in rae2822_M0.850_a1.25 rae2822_M0.750_a1.25 biconvex_M0.750_a1.25 dfvlr_r4_M0.750_a1.25 naca0012_M0.775_a1.25; do
  f="$c/postProcessing/forces/0/coefficient.dat"
  if [ -s "$f" ]; then
    echo "$c: $(tail -1 $f | cut -f1,5)"
  else
    echo "$c: no data yet"
  fi
done
echo "== worker logs:"
tail -2 batchv3b.log 2>/dev/null
tail -2 batchw5.log 2>/dev/null
tail -2 batchv2.log 2>/dev/null
tail -2 batch2.log 2>/dev/null
echo "== results2 .dat files:"
ls -la results2/*.dat 2>/dev/null | awk '{print $5, $9}' | sort -k2
echo "== memory:"
free -m | head -2
