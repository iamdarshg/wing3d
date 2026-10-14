echo "=== controlDict numerics (naca0012, known-good) ==="
grep -a -E 'maxCo|maxDeltaT|deltaT|adjustTimeStep|endTime' /caseCampaign/campaign/naca0012_M0.900_a1.25/system/controlDict 2>/dev/null
echo ""
echo "=== fvSchemes flux ==="
grep -a -E 'fluxScheme|ddtSchemes|default' /caseCampaign/campaign/naca0012_M0.900_a1.25/system/fvSchemes 2>/dev/null | head -8
echo ""
echo "=== dfvlr snappy levels ==="
grep -a -E 'refinementSurfaces|level|maxLocalCells|maxGlobalCells' /caseCampaign/campaign/dfvlr_r4_M1.000_a1.25/system/snappyHexMeshDict 2>/dev/null | head -8
echo ""
echo "=== dfvlr mesh cells all cases ==="
for c in /caseCampaign/campaign/dfvlr_r4_M*; do
  echo "$(basename $c): $(grep -a -m1 'cells:' $c/log.check 2>/dev/null || echo no-check)"
done