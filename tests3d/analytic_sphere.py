"""Shared sphere mesh builder for the 3D suite."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
from wing3d.geometry import Mesh


def sphere_mesh(nlat=12, nlon=24):
    verts = []
    for i in range(nlat + 1):
        th = np.pi * i / nlat
        for j in range(nlon):
            ph = 2 * np.pi * j / nlon
            verts.append([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph),
                          np.cos(th)])
    verts = np.array(verts)
    faces = []
    for i in range(nlat):
        for j in range(nlon):
            a = i * nlon + j
            b = i * nlon + (j + 1) % nlon
            c = (i + 1) * nlon + j
            d = (i + 1) * nlon + (j + 1) % nlon
            if i > 0:
                faces.append([a, c, b])
            if i < nlat - 1:
                faces.append([b, c, d])
    return Mesh(verts, faces)
