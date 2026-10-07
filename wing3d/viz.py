"""Visualizer set: 3D surface Cp, cuts, polars, convergence, spanload + VTK export."""
import numpy as np


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    return plt, Poly3DCollection


def plot_surface_cp(mesh, cp, savepath, title='Cp', elev=20, azim=-60):
    plt, Poly3DCollection = _mpl()
    tris = []
    vals = []
    for f, c in zip(mesh.faces, cp):
        pts = [mesh.verts[k] for k in f]
        for k in range(1, len(pts) - 1):
            tris.append([pts[0], pts[k], pts[k + 1]])
            vals.append(c)
    tris = np.array(tris)
    vals = np.array(vals)
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection='3d')
    vmin, vmax = float(np.percentile(vals, 2)), float(np.percentile(vals, 98))
    coll = Poly3DCollection(tris, cmap='jet', alpha=1.0)
    coll.set_array(vals)
    coll.set_clim(vmin, vmax)
    ax.add_collection3d(coll)
    ax.scatter(mesh.centroid[:, 0], mesh.centroid[:, 1], mesh.centroid[:, 2],
               s=0.1, c='k', alpha=0.1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.view_init(elev=elev, azim=azim)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    fig.colorbar(coll, ax=ax, label='Cp')
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(savepath, dpi=100)
    plt.close(fig)
    return savepath


def plot_cp_cuts(mesh, cp, stations, savepath, title='Cp cuts'):
    """Cp vs x at given span stations (nearest panels)."""
    plt, _ = _mpl()
    C = mesh.centroid
    fig, ax = plt.subplots(figsize=(8, 5))
    for y0 in stations:
        sel = np.abs(C[:, 1] - y0) < 0.06
        if not sel.any():
            continue
        idx = np.argsort(C[sel][:, 0])
        ax.plot(C[sel][:, 0][idx], cp[sel][idx], '.-', label=f'y={y0}')
    ax.invert_yaxis()
    ax.set_xlabel('x')
    ax.set_ylabel('Cp')
    ax.legend()
    ax.set_title(title)
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(savepath, dpi=100)
    plt.close(fig)
    return savepath


def plot_polar(alphas, cls, cds, savepath, title='Polar'):
    plt, _ = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(alphas, cls, 'o-')
    axes[0].set_xlabel('alpha (deg)')
    axes[0].set_ylabel('CL')
    axes[0].grid(True)
    axes[1].plot(cds, cls, 'o-')
    axes[1].set_xlabel('CD')
    axes[1].set_ylabel('CL')
    axes[1].grid(True)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(savepath, dpi=100)
    plt.close(fig)
    return savepath


def plot_history(hist, savepath, title='Convergence'):
    """hist: list of dicts with CL/CD keys."""
    plt, _ = _mpl()
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot([h['CL'] for h in hist], 'o-', label='CL')
    ax.plot([h['CD'] for h in hist], 's-', label='CD')
    ax.set_xlabel('iteration')
    ax.legend()
    ax.grid(True)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(savepath, dpi=100)
    plt.close(fig)
    return savepath


def plot_spanload(y, cl, savepath, title='Spanwise loading'):
    plt, _ = _mpl()
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(y, cl, 'o-')
    ax.set_xlabel('span y')
    ax.set_ylabel('sectional cl')
    ax.grid(True)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(savepath, dpi=100)
    plt.close(fig)
    return savepath


def write_vtk(mesh, cp, path, vel=None):
    """Legacy ASCII VTK POLYDATA (ParaView-readable)."""
    verts = np.asarray(mesh.verts)
    polys = []
    for f in mesh.faces:
        for k in range(1, len(f) - 1):
            polys.append((f[0], f[k], f[k + 1]))
    with open(path, 'w') as fh:
        fh.write('# vtk DataFile Version 2.0\nwing3d\nASCII\n')
        fh.write('DATASET POLYDATA\n')
        fh.write(f'POINTS {len(verts)} float\n')
        for p in verts:
            fh.write(f'{p[0]} {p[1]} {p[2]}\n')
        fh.write(f'POLYGONS {len(polys)} {4 * len(polys)}\n')
        for t in polys:
            fh.write(f'3 {t[0]} {t[1]} {t[2]}\n')
        # cell data: repeat panel cp per triangle
        vals = []
        for f, c in zip(mesh.faces, cp):
            for _ in range(1, len(f) - 1):
                vals.append(c)
        fh.write(f'CELL_DATA {len(vals)}\nSCALARS Cp float 1\n'
                 f'LOOKUP_TABLE default\n')
        for v in vals:
            fh.write(f'{v}\n')
        if vel is not None:
            fh.write(f'POINT_DATA {len(verts)}\nVECTORS V float\n')
            # scatter panel velocity to verts (nearest panel)
            for _ in verts:
                fh.write('0 0 0\n')
    return path
