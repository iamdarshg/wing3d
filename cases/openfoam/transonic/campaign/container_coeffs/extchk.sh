source /usr/lib/openfoam/openfoam2406/etc/bashrc
for g in naca0012 rae2822; do
  cd /caseCampaign/campaign/${g}_M0.850_a1.25
  echo "=== $g airfoil patch extent ==="
  surfaceCheck constant/triSurface/${g}.stl 2>/dev/null | grep -a -E 'Bounding|box' | head -3
done
echo "=== naca STL span ==="
cd /caseCampaign/campaign/naca0012_M0.850_a1.25
source /usr/lib/openfoam/openfoam2406/etc/bashrc
surfaceCheck constant/triSurface/naca0012.stl 2>&1 | grep -a -A2 Bounding | head -4