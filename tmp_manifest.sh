#!/bin/bash
A=/dfvlr_coarse
B=/caseCampaign/probes/dfvlr_r4_M0.925_coarse45k
echo "=du bytes="
du -sb "$A/constant/polyMesh" "$B/constant/polyMesh"
echo "=sets dir A="
ls -la "$A/constant/polyMesh/sets/" 2>&1 | head -14
echo "=sets dir B="
ls -la "$B/constant/polyMesh/sets/" 2>&1 | head -14
cd "$A" || exit 1
find . -type f -size -50M | sort | xargs -d '\n' md5sum | sort -k2 > /tmp/manA.txt
cd "$B" || exit 1
find . -type f -size -50M | sort | xargs -d '\n' md5sum | sort -k2 > /tmp/manB.txt
echo "=manifest counts A/B="
wc -l /tmp/manA.txt /tmp/manB.txt
echo "=path-set diff (< A only, > B only)="
diff <(awk '{print $NF}' /tmp/manA.txt) <(awk '{print $NF}' /tmp/manB.txt) | head -30
echo "=hash diff on common paths="
join -1 2 -2 2 <(sort -k2 /tmp/manA.txt) <(sort -k2 /tmp/manB.txt) | awk '$2 != $3 {print $1}' | head -30
echo "=controlDict diff="
diff "$A/system/controlDict" "$B/system/controlDict" | head -20
echo "=fvSchemes diff="
diff "$A/system/fvSchemes" "$B/system/fvSchemes" | head -20
echo "=x3 verdict="
grep -a -E "Mesh OK|Failed|FOAM exiting" /tmp/x3.log
echo "=x4 verdict="
grep -a -E "Mesh OK|Failed|FOAM exiting" /tmp/x4.log
echo "=cmA/cmB exec lines="
grep -a -E "Exec :" /tmp/cmA.log /tmp/cmB.log
