"""Camber/thickness of master vs deployed case STLs (section slice)."""
import sys
sys.path.insert(0, r'D:\CodeProjects\cfd')
from wing3d.geometry import read_stl
from wing3d.tsdcli import extract_section

paths = [
    (r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl\rae2822.stl', 'master rae'),
    (r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl\dfvlr_r4.stl', 'master dfvlr'),
    (r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl\naca0012.stl', 'master naca'),
    (r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl\biconvex.stl', 'master biconvex'),
    (r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign\rae2822_M0.850_a1.25\constant\triSurface\rae2822.stl', 'deployed rae'),
    (r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign\dfvlr_r4_M0.850_a1.25\constant\triSurface\dfvlr_r4.stl', 'deployed dfvlr'),
    (r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign\biconvex_M0.850_a1.25\constant\triSurface\biconvex.stl', 'deployed biconvex'),
    (r'D:\CodeProjects\cfd\cases\openfoam\transonic\campaign\naca0012_M0.850_a1.25\constant\triSurface\naca0012.stl', 'deployed naca'),
]
for p, tag in paths:
    try:
        m = read_stl(p)
        s = extract_section(m)
        print('%-16s t/c=%.4f camber=%.4f chord=%.3f station=%.2f extent=%s'
              % (tag, s['thick_max'], s['camber_max'], s['chord'],
                 s['station'], ['%.1f' % e for e in s['mesh_extent']]))
    except Exception as e:
        print('%-16s FAIL %s' % (tag, e))
