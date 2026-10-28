"""Direct SU2 ASCII grid -> OpenFOAM polyMesh converter.

Why: `plot3dToFoam` folds the cells on the NASA TMR grids (measured max
skewness 183, aspect 212769, even with NASA's own input.nml), while the
generator's SU2 export carries explicit element connectivity -- no grid
ordering ambiguity.  This writes points/faces/owner/neighbour/boundary
straight from the SU2 hexahedra.

Usage:
  python tests3d/su2_to_of.py <grid.su2> <case_dir> [--check]

Only hexahedra (type 12) and quads (type 9) are supported, which is what
hcf_wing igrid_type=5 (structured) emits.
"""
import os
import re
import sys
from collections import defaultdict


def parse_su2(path):
    pts, elems, markers = [], [], {}
    cur_marker = None
    with open(path, 'r', errors='ignore') as f:
        lines = f.read().splitlines()
    i, n = 0, len(lines)
    while i < n:
        ln = lines[i].strip()
        if ln.startswith('NPOIN='):
            npt = int(ln.split('=')[1])
            i += 1
            for _ in range(npt):
                v = lines[i].split()
                pts.append((float(v[0]), float(v[1]), float(v[2])))
                i += 1
            continue
        if ln.startswith('NELEM='):
            ne = int(ln.split('=')[1])
            i += 1
            for _ in range(ne):
                v = [int(t) for t in lines[i].split()]
                elems.append((v[0], v[1:-1]))   # (type, nodes)
                i += 1
            continue
        if ln.startswith('NMARK='):
            nm = int(ln.split('=')[1])
            i += 1
            for _ in range(nm):
                tag = lines[i].split('=')[-1].strip()
                i += 1
                m = lines[i].split('=')[-1].strip()
                i += 1
                lst = []
                for _ in range(int(m)):
                    v = [int(t) for t in lines[i].split()]
                    lst.append((v[0], v[1:]))
                    i += 1
                markers[tag] = lst
            continue
        i += 1
    return pts, elems, markers


# VTK/SU2 hexahedron faces (outward normals)
HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
QUAD_FACES = [(0, 1, 2, 3)]


def patch_type(tag):
    t = tag.lower()
    if 'wall' in t or 'solid' in t or 'wing' in t or 'viscous' in t:
        return 'wall'
    if 'symm' in t:
        return 'symmetry'
    return 'patch'


def convert(su2_path, case_dir):
    pts, elems, markers = parse_su2(su2_path)
    print('points %d, elements %d, markers %s'
          % (len(pts), len(elems), {k: len(v) for k, v in markers.items()}))

    # boundary faces: sorted-key -> marker name
    bnd_key = {}
    for tag, lst in markers.items():
        for typ, nodes in lst:
            faces = HEX_FACES if typ == 12 else QUAD_FACES
            # a marker polygon may be a quad directly
            if typ == 9:
                bnd_key[tuple(sorted(nodes))] = tag
            else:
                bnd_key[tuple(sorted(nodes))] = tag

    face_id = {}          # sorted key -> face index
    faces = []            # list of vertex tuples (orientation as stored)
    owner = []
    neigh = []
    bnd_face_tag = []

    def add_face(vs, cell, key):
        if key in face_id:
            fi = face_id[key]
            if len(neigh) > fi and neigh[fi] == -1:
                neigh[fi] = cell
            return
        fi = len(faces)
        face_id[key] = fi
        faces.append(vs)
        owner.append(cell)
        neigh.append(-1)
        bnd_face_tag.append(None)

    for cell, (typ, nodes) in enumerate(elems):
        if typ != 12:
            continue
        for f in HEX_FACES:
            vs = tuple(nodes[j] for j in f)
            add_face(vs, cell, tuple(sorted(vs)))

    # boundary faces: those with neighbour == -1
    n_bound = sum(1 for x in neigh if x == -1)
    print('faces %d, interior %d, boundary %d'
          % (len(faces), len(faces) - n_bound, n_bound))

    patches = defaultdict(list)
    for fi, nb in enumerate(neigh):
        if nb != -1:
            continue
        key = tuple(sorted(faces[fi]))
        tag = bnd_key.get(key, 'defaultFaces')
        patches[tag].append(fi)

    # reorder faces: interior first (owner<neighbour), then boundary per patch
    order = [fi for fi, nb in enumerate(neigh) if nb != -1]
    patch_ranges = {}
    for tag, lst in patches.items():
        start = len(order)
        order.extend(lst)
        patch_ranges[tag] = (start, len(lst))

    os.makedirs(os.path.join(case_dir, 'constant', 'polyMesh'), exist_ok=True)
    pm = os.path.join(case_dir, 'constant', 'polyMesh')

    with open(os.path.join(pm, 'points'), 'w') as f:
        f.write('FoamFile { version 2.0; format ascii; class vectorField; '
                'object points; }\n%d\n(\n' % len(pts))
        for p in pts:
            f.write('(%r %r %r)\n' % p)
        f.write(')\n')

    with open(os.path.join(pm, 'faces'), 'w') as f:
        f.write('FoamFile { version 2.0; format ascii; class faceList; '
                'object faces; }\n%d\n(\n' % len(order))
        for fi in order:
            vs = faces[fi]
            f.write('%d(%s)\n' % (len(vs), ' '.join(str(v) for v in vs)))
        f.write(')\n')

    with open(os.path.join(pm, 'owner'), 'w') as f:
        f.write('FoamFile { version 2.0; format ascii; class labelList; '
                'note "nPoints: %d nCells: %d"; object owner; }\n%d\n(\n'
                % (len(pts), len(elems), len(order)))
        for fi in order:
            f.write('%d\n' % owner[fi])
        f.write(')\n')

    with open(os.path.join(pm, 'neighbour'), 'w') as f:
        f.write('FoamFile { version 2.0; format ascii; class labelList; '
                'note "nPoints: %d nCells: %d"; object neighbour; }\n%d\n(\n'
                % (len(pts), len(elems), len(order) - n_bound))
        for fi in order:
            if neigh[fi] != -1:
                f.write('%d\n' % neigh[fi])
        f.write(')\n')

    with open(os.path.join(pm, 'boundary'), 'w') as f:
        f.write('FoamFile { version 2.0; format ascii; class polyBoundaryMesh; '
                'object boundary; }\n%d\n(\n' % len(patch_ranges))
        for tag, (start, cnt) in patch_ranges.items():
            f.write('    %s\n    {\n        type            %s;\n'
                    '        nFaces          %d;\n        startFace       %d;\n'
                    '    }\n' % (tag, patch_type(tag), cnt, start))
        f.write(')\n')
    print('wrote', pm)
    return pm


if __name__ == '__main__':
    su2 = sys.argv[1]
    case = sys.argv[2]
    convert(su2, case)
