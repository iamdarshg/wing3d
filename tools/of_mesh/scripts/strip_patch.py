import re
import subprocess
import tempfile
import os

CONT = 'oftrans'
remote = '/work/tune1/constant/polyMesh/boundary'
fd, tmp = tempfile.mkstemp()
os.close(fd)
subprocess.run(['docker', 'cp', f'{CONT}:{remote}', tmp], check=True)
txt = open(tmp, errors='ignore').read()

# parse: N\n(\n blocks \n)
head, rest = txt.split('\n(', 1)
body, tail = rest.rsplit('\n)', 1)
blocks = re.findall(r'\n    (\w+)\s*\n    \{([^}]*)\}', body)
kept = []
for name, b in blocks:
    nf = re.search(r'nFaces\s+(\d+)', b)
    if name == 'defaultFaces' or (nf and int(nf.group(1)) == 0):
        continue
    kept.append((name, b))
out = head + '\n(\n'
for name, b in kept:
    out += '    %s\n    {%s}\n' % (name, b)
out += ')\n' + tail
out = re.sub(r'^\d+$', str(len(kept)), out, count=1, flags=re.M)
# fix the count line (first standalone integer after the header block)
lines = out.split('\n')
for i, ln in enumerate(lines):
    if ln.strip().isdigit():
        lines[i] = str(len(kept))
        break
out = '\n'.join(lines)
open(tmp, 'w').write(out)
subprocess.run(['docker', 'cp', tmp, f'{CONT}:{remote}'], check=True)
os.unlink(tmp)
print('patches kept:', [n for n, _ in kept])
