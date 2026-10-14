c=biconvex_M0.850_a1.25
f=$(find /caseCampaign/campaign/$c/postProcessing -name coefficient.dat 2>/dev/null | head -1)
echo "Cl trajectory (t, Cl) for $c:"
python3 - 2>/dev/null << EOF || awk 'NR>2 && NR%20==1 {print $1, $5}' "$f" | tail -12
EOF
python3 -c "
import re
f='''$(find /caseCampaign/campaign/$c/postProcessing -name coefficient.dat 2>/dev/null | head -1)'''
print('using', f)
lines = open(f).read().splitlines()
hdr = 0
import sys
for i,ln in enumerate(lines):
    if ln.strip() and not ln.lstrip().startswith('#'):
        hdr = i
        break
rows = [ln.split() for ln in lines[hdr:] if ln.split()]
print('last rows (t Cl Cd):')
for r in rows[-6:]:
    print('  t=%.4f Cl=%+.4f Cd=%.4f' % (float(r[0]), float(r[4]), float(r[6])))
print('first rows:')
for r in rows[1:6]:
    print('  t=%.4f Cl=%+.4f Cd=%.4f' % (float(r[0]), float(r[4]), float(r[6])))
" 2>&1 | tail -15