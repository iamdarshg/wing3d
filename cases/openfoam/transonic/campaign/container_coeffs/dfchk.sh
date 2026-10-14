echo "=== dfvlr M0.85 mesh check (should be 86145) ==="
grep -a -m1 'cells:' /caseCampaign/campaign/dfvlr_r4_M0.850_a1.25/log.check 2>/dev/null || echo "no check"
echo "=== runner alive? ==="
pgrep -af dfvlr_runner || echo "runner gone"
echo "=== dfvlr runner stdout ==="
cat /caseCampaign/dfvlr_runner.log 2>/dev/null | tail -8
echo "=== all rhoCentral ==="
pgrep -af rhoCentral | awk '{print $1}' | wc -l