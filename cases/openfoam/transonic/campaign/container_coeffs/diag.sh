echo "=== done markers ==="
ls /caseCampaign/campaign/results/ | grep done | sort
echo ""
echo "=== dfvlr_r4 log tails (killed or ran?) ==="
for c in dfvlr_r4_M0.875_a1.25 dfvlr_r4_M1.000_a1.25; do
  echo "--- $c"
  grep -a 'Killed\|FOAM exiting\|FOAM FATAL\|Time = ' /caseCampaign/campaign/$c/log.run 2>/dev/null | tail -3
  echo "  last force row: $(find /caseCampaign/campaign/$c -name coefficient.dat 2>/dev/null | head -1 | xargs tail -1 2>/dev/null)"
done
echo ""
echo "=== how many .dat have real t=0.25 data ==="
grep -l '^0.25' /caseCampaign/campaign/results/*.dat 2>/dev/null | wc -l
echo "=== rhoCentral running ==="
pgrep -c rhoCentralFoam || echo 0