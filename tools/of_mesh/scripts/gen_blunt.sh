#!/bin/bash
set -e
cd /work/tmr/wing_release/wing_release_072319
rm -f wing_strct.1.p3d wing_strct.1.su2
python3 - <<'EOF' 2>/dev/null || sed -i "s/tip_type = 'round'/tip_type = 'blunt'/" input.nml
EOF
sed -i "s/tip_type = 'round'/tip_type = 'blunt'/" input.nml
sed -i 's/use_s_blunt = T/use_s_blunt = F/' input.nml
grep -E 'tip_type|igrid_type' input.nml
./hcf_wing | grep -E 'Root chord|Aspect|Reference area|success'
ls -la wing_strct.1.p3d 2>/dev/null || echo "NO P3D"
