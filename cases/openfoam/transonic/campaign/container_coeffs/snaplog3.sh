echo "=== rae snappy: first 40 lines ==="
head -40 /caseCampaign/campaign/rae2822_M0.900_a1.25/log.snappy
echo ""
echo "=== STL validity check ==="
for g in rae2822 biconvex dfvlr_r4; do
  sz=$(stat -c%s /caseCampaign/campaign/${g}_M0.900_a1.25/constant/triSurface/${g}.stl)
  cnt=$(head -c 84 /caseCampaign/campaign/${g}_M0.900_a1.25/constant/triSurface/${g}.stl | tail -c 4 | od -An -tu4 | tr -d ' ')
  echo "$g size=$sz declared_tris=$cnt expected=$((84+$cnt*50))"
done