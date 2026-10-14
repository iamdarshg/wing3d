#!/bin/bash
cd /caseCampaign/campaign
for c in naca0012_M* rae2822_M* biconvex_M* dfvlr_r4_M*; do
  f=$(find "$c/postProcessing" -name coefficient.dat 2>/dev/null | head -1)
  if [ -n "$f" ]; then
    cp "$f" "/tmp/${c}.coeff"
  fi
done
echo "coeff files copied: $(ls /tmp/*.coeff 2>/dev/null | wc -l)"