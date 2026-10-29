#!/bin/bash
cd /work/tmr/wing_release/wing_release_072319 || exit 1
rm -f wing_strct.1.p3d
sed -i 's/ugrid_file_unformatted = T/ugrid_file_unformatted = F/' input.nml
grep -E 'ugrid_file_unformatted|igrid_type|tip_type' input.nml
./hcf_wing | grep -E 'Reference area|Aspect|success'
ls -la wing_strct.1.p3d 2>/dev/null || echo "NO P3D"
