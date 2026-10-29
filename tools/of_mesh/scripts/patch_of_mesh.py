import io, os
p = r'D:\CodeProjects\cfd\tests3d\of_mesh.py'
s = open(p, encoding='utf-8').read()

if 'def cp(*parts):' not in s:
    s = s.replace(
        "GEN_BIN = os.path.join(GEN_DIR, 'hcf_wing')",
        "GEN_BIN = GEN_DIR + '/hcf_wing'\n\n\n"
        "def cp(*parts):\n"
        "    \"\"\"POSIX join for container paths.\"\"\"\n"
        "    return '/'.join(str(x).rstrip('/') for x in parts)")

repl = [
    ("os.path.join(GEN_DIR, 'input.nml')", "cp(GEN_DIR, 'input.nml')"),
    ("os.path.join(GEN_DIR, a.dat_name)", "cp(GEN_DIR, a.dat_name)"),
    ("os.path.join(out, 'system/controlDict')", "cp(out, 'system/controlDict')"),
    ("os.path.join(out, 'constant/polyMesh/points')", "cp(out, 'constant/polyMesh/points')"),
    ("os.path.join(out, 'constant/polyMesh/faces')", "cp(out, 'constant/polyMesh/faces')"),
    ("os.path.join(out, 'constant/polyMesh/boundary')", "cp(out, 'constant/polyMesh/boundary')"),
    ("os.path.join(out, 'system', k)", "cp(out, 'system', k)"),
    ("os.path.join(out, 'constant', k)", "cp(out, 'constant', k)"),
    ("os.path.join(out, '0', k)", "cp(out, '0', k)"),
]
for a, b in repl:
    s = s.replace(a, b)

open(p, 'w', encoding='utf-8').write(s)
print('patched ok')
