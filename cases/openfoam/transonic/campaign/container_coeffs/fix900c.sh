#!/bin/bash
# Ensure M0.900 uses maxCo 0.3 (matches the converged dfvlr cases).
sed -i 's/maxCo 0.15/maxCo 0.3/' /caseCampaign/campaign/dfvlr_r4_M0.900_a1.25/system/controlDict
grep -a -E 'maxCo|deltaT|maxDeltaT' /caseCampaign/campaign/dfvlr_r4_M0.900_a1.25/system/controlDict
pgrep -af rhoCentralFoam | head -3