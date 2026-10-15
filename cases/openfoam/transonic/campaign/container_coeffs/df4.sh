echo "=== dfvlr state ==="
cd /caseCampaign/campaign
for c in dfvlr_r4_M0.900_a1.25 dfvlr_r4_M0.950_a1.25 dfvlr_r4_M0.975_a1.25 dfvlr_r4_M1.000_a1.25; do
  f=$(find $c/postProcessing -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then
    echo "$c: $(tail -1 $f | awk '{print \"t=\"$1\" CL=\"$5\" CD=\"$7}')"
  else
    echo "$c: no-data ($(ls $c/log.run >/dev/null 2>&1 && echo has-log))"
  fi
done
echo "=== running ==="
pgrep -af rhoCentralFoam | wc -l
echo "=== done total ==="
ls results/*.done | wc -l