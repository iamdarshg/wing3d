"""Final validation plots: wing3d vs public data vs OpenFOAM RANS."""
import numpy as np
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- 1. wing polar ---
alphas = np.array([0, 2, 4, 6, 8])
# wing3d inviscid doublet (24x12, te_gap 0.002) — from validate runs
w_inv_cl = np.array([0.0, 0.1844, 0.3654, 0.5432, 0.7181])
w_inv_cd = np.array([0.00773, 0.01254, 0.01280, 0.00821, -0.00149])
# wing3d viscous-coupled (24x12, Re 3e6) — measured at 4deg only
w_vis_cl = np.array([0.3636])
w_vis_cd = np.array([0.0264])
w_vis_a = np.array([4.0])
# lifting line AR6 (a0=2pi)
a0 = 2 * np.pi
AR = 6.0
ll_cl = a0 * np.radians(alphas) / (1 + a0 / (np.pi * AR))
# VLM (computed)
vlm_cl = np.array([0.0, 0.196, 0.3921, 0.588, 0.784])
# OpenFOAM RANS (SST, this study; wing mesh-limited -> stalled branch)
rans_cl, rans_cd = 0.113, 0.0224

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(alphas, w_inv_cl, 'o-', label='wing3d inviscid')
axes[0].plot(w_vis_a, w_vis_cl, 's', ms=8, label='wing3d viscous-coupled')
axes[0].plot(alphas, ll_cl, '--', label='lifting line AR6')
axes[0].plot(alphas, vlm_cl, ':', label='VLM (in-repo)')
axes[0].plot([4], [rans_cl], 'rx', ms=10, mew=2,
             label='OpenFOAM RANS (mesh-limited)')
axes[0].set_xlabel('alpha (deg)')
axes[0].set_ylabel('CL')
axes[0].legend(fontsize=8)
axes[0].grid(True)
axes[0].set_title('NACA0012 AR6 lift')
axes[1].plot(w_inv_cd, w_inv_cl, 'o-', label='wing3d inviscid')
axes[1].plot(w_vis_cd, w_vis_cl, 's', ms=8, label='wing3d viscous-coupled')
axes[1].plot([rans_cd], [rans_cl], 'rx', ms=10, mew=2, label='OF RANS')
axes[1].set_xlabel('CD')
axes[1].set_ylabel('CL')
axes[1].legend(fontsize=8)
axes[1].grid(True)
axes[1].set_title('drag polar')
fig.tight_layout()
fig.savefig('D:\\CodeProjects\\cfd\\cases\\openfoam\\polar_compare.png', dpi=110)
print('polar saved')

# --- 2. sphere Cp ---
from tests3d.analytic_sphere import sphere_mesh
from wing3d.solver import solve
m = sphere_mesh(12, 24)
res = solve(m, [], [1, 0, 0])
C = m.centroid
th = np.arccos(np.clip(C[:, 0], -1, 1))
o = np.argsort(th)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(th[o] * 180 / np.pi, res['cp'][o], '.', ms=2, label='wing3d')
ax.plot(th[o] * 180 / np.pi, 1 - 2.25 * np.sin(th[o]) ** 2, 'k-',
        label='analytic')
ax.set_xlabel('theta from stagnation (deg)')
ax.set_ylabel('Cp')
ax.legend()
ax.grid(True)
ax.set_title('sphere Cp (12x24, rms 0.08)')
fig.tight_layout()
fig.savefig('D:\\CodeProjects\\cfd\\cases\\openfoam\\sphere_compare.png', dpi=110)
print('sphere saved')
