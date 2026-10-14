echo "=== rae2822 mesh cells ==="
for c in /caseCampaign/campaign/rae2822_M*; do
  echo "$(basename $c): $(grep -a -m1 'cells:' $c/log.check 2>/dev/null || echo no-check)"
done
echo ""
echo "=== rae2822 done markers ==="
ls /caseCampaign/campaign/results/rae2822_*.done 2>/dev/null || echo none
echo ""
echo "=== current rhoCentral (4 workers + dfvlr runner = 5?) ==="
pgrep -c rhoCentralFoam || echo 0
echo "=== dfvlr runner progress ==="
tail -3 /caseCampaign/dfvlr_runner.log 2>/dev/null