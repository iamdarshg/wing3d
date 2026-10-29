"""A17 estimator lab: all variants from the cached fields (no re-solve).

    python tmp_a17_lab.py [tag]     # tag default nz16, dumps one case
"""
import glob
import json
import os
import sys

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)


def thickness_slope(xc, t=0.12):
    x = np.clip(np.asarray(xc, float), 0, 1)
    return 5 * t * (0.2969 / (2 * np.sqrt(np.maximum(x, 1e-9))) - 0.1260
                    - 2 * 0.3516 * x + 3 * 0.2843 * x ** 2 - 4 * 0.1036 * x ** 3)


class C:
    def __init__(self, path):
        d = np.load(path, allow_pickle=False)
        self.__dict__.update({k: d[k] for k in d.files})
        self.free_mu = bool(self.free_mu[0])
        self.tag = str(self.tag[0])
        self.CL_fp_stored = float(self.CL_fp[0])
        self.nx, self.ny, self.nz = self.phi.shape
        self.wing_ix = [int(v) for v in self.wing_ix]
        self.wing_jy = [int(v) for v in self.wing_jy]

    # --- the jump the residual applies, on wing nodes -------------------
    def jump(self):
        if self.free_mu:
            xc = np.clip(self.xc[self.wing_ix] / self.chord, 0, 1)
            J = np.empty((len(self.wing_ix), len(self.wing_jy)))
            for p in range(len(self.wing_jy)):
                J[:, p] = np.interp(xc, self.x_ctrl, self.mu_b[:, p])
            return J
        return np.asarray(self.gamma_b)[None, :] * 1.0

    # --- the model's own analytic loading dJ/dx at the ctrl stations ---
    def loading(self):
        """gamma(x) = d(mu)/dx per chordwise control station, (nc, ng)."""
        return np.gradient(self.mu_b, self.x_ctrl, axis=0)

    def gradx(self, k):
        g = np.zeros((self.nx, self.ny))
        g[1:-1] = (self.phi[2:, :, k] - self.phi[:-2, :, k]) / (
            self.x[2:] - self.x[:-2])[:, None]
        return g

    def extrap(self, k1, k2):
        z1, z2 = self.z[k1], self.z[k2]
        w = abs(z2) / max(abs(z2 - z1), 1e-12)
        return w * self.gradx(k1) + (1 - w) * self.gradx(k2)

    def cp_of_g(self, gx):
        gm = self.gamma
        q2 = np.maximum((1.0 + gx) ** 2, 1e-9)
        rho = np.maximum(np.maximum(
            1 + 0.5 * (gm - 1) * self.Minf ** 2 * (1 - q2), 1e-9)
            ** (1 / (gm - 1)), 0.2)
        return 2.0 / (gm * self.Minf ** 2) * (rho ** gm - 1.0)

    def integ(self, cpu, cpl):
        CL = CD = 0.0
        for i in self.wing_ix:
            xc = self.xc[i] / self.chord
            st = thickness_slope(xc, self.thick)
            for j in self.wing_jy:
                dA = self.dx[min(i, len(self.dx) - 1)] * self.dy[min(j,
                                                                    len(self.dy) - 1)]
                CL += (cpl[i, j] - cpu[i, j]) * dA
                CD += (cpu[i, j] - cpl[i, j]) * st * dA
        S = self.chord * self.span
        return CL / S, CD / S


def dJ_full(c, mode):
    """d(jump)/dx on the (nx, ny) grid, zero outside the wing."""
    Jf = np.zeros((c.nx, c.ny))
    Jf[np.ix_(c.wing_ix, c.wing_jy)] = c.jump()
    d = np.zeros_like(Jf)
    if mode == 'central':           # same operator as gradx
        d[1:-1] = (Jf[2:, :] - Jf[:-2, :]) / (c.x[2:] - c.x[:-2])[:, None]
    elif mode == 'face':            # face differences -> telescopes exactly
        iu = np.asarray(c.wing_ix)
        d[iu, :] = (Jf[np.minimum(iu + 1, c.nx - 1), :]
                    - Jf[np.maximum(iu - 1, 0), :]) / (
            c.x[np.minimum(iu + 1, c.nx - 1)]
            - c.x[np.maximum(iu - 1, 0)])[:, None]
    elif mode == 'onesided':
        iu = np.asarray(c.wing_ix)
        for n, i in enumerate(iu):
            if n + 1 < len(iu):
                d[i, :] = (Jf[iu[n + 1], :] - Jf[i, :]) / (c.x[iu[n + 1]]
                                                          - c.x[i])
            else:
                d[i, :] = (Jf[i, :] - Jf[iu[n - 1], :]) / (c.x[i]
                                                           - c.x[iu[n - 1]])
    else:
        raise SystemExit(mode)
    return d


def upper_local(c):
    """Upper-surface dphi/dx at z=0, extrapolated ONLY from the two planes
    nearest the sheet (k0 is the upper side of the cut face k0-1).

    Bounded: |0 - z[k0]| <= one cell, so the extrapolation weight on the
    near plane is <= 1 and there is no dz-dependent amplification.
    """
    k0 = c.k0
    k1 = min(k0 + 1, c.nz - 1)
    g0, g1 = c.gradx(k0), c.gradx(k1)
    dz = c.z[k1] - c.z[k0]
    t = 0.0 if abs(dz) < 1e-14 else (0.0 - c.z[k0]) / dz
    t = float(np.clip(t, -1.0, 1.0))
    return g0 + t * (g1 - g0)


def gam_nodal(c):
    """d(jump)/dx on the wing nodes -- the exact sheet identity's RHS.

    No zero padding (that put a -0.18 spike on the last station and cost
    80% of CL); central in the interior, one-sided at the wing ends.
    """
    J = c.jump()                                    # (nix, njy)
    x = c.x[np.asarray(c.wing_ix)]
    return np.gradient(J, x, axis=0, edge_order=1)


def gam_analytic(c):
    """gamma(x) = dmu/dx of the model's own sheet at the ctrl stations."""
    load = c.loading()
    xc = np.clip(c.xc[c.wing_ix] / c.chord, 0, 1)
    out = np.empty((len(c.wing_ix), len(c.wing_jy)))
    for p in range(len(c.wing_jy)):
        out[:, p] = np.interp(xc, c.x_ctrl / c.chord, load[:, p])
    return out


def gam_cell(c):
    """Cell-integrated jump derivative -- telescopes EXACTLY to Gamma.

    gamma_i = (J_face(i+) - J_face(i-)) / dx_i with the faces the cell
    averages of the nodal jump and the end faces pinned to the physical
    sheet ends (0 ahead of the LE, J_TE behind the TE).  Then
    sum_i gamma_i dx_i = J_TE - 0 = Gamma by construction, which removes
    the leading-edge cell that a nodal finite difference drops (44% of
    the circulation on this 6-node chord).
    """
    J = c.jump()                                    # (nix, njy)
    n = len(c.wing_ix)
    Jf = np.empty((n + 1, J.shape[1]))
    Jf[0] = 0.0
    Jf[1:n] = 0.5 * (J[:-1] + J[1:])
    Jf[n] = J[-1]
    x = c.x[np.asarray(c.wing_ix)]
    dxi = np.empty(n)
    dxi[:-1] = np.diff(x)
    dxi[-1] = dxi[-2] if n > 1 else 1.0
    return (Jf[1:] - Jf[:-1]) / dxi[:, None]


def est(c, how):
    k0 = c.k0
    ku1, ku2 = min(k0 + 1, c.nz - 1), min(k0 + 2, c.nz - 1)
    kl1, kl2 = max(k0 - 1, 0), max(k0 - 2, 0)
    gU = c.extrap(ku1, ku2)
    gL = c.extrap(kl1, kl2)
    if how == 'fp':
        pass
    elif how in ('jump_c', 'jump_f', 'jump_o'):
        m = {'jump_c': 'central', 'jump_f': 'face', 'jump_o': 'onesided'}[how]
        d = dJ_full(c, m)
        gU = c.gradx(k0)                 # node k0 is ON the sheet (upper)
        gL = gU - d
    elif how == 'load':
        # the model's own sheet strength: dCp = 2*gamma(x), gamma from the
        # mu sheet's own derivative, interpolated onto the wing nodes
        gU, gj = upper_local(c), gam_analytic(c)
        gL = np.array(gU)
        gL[np.ix_(c.wing_ix, c.wing_jy)] = gU[np.ix_(c.wing_ix,
                                                    c.wing_jy)] - gj
    elif how == 'decont2':
        gU = upper_local(c)
        gj = gam_cell(c)
        gL = np.array(gU)
        gL[np.ix_(c.wing_ix, c.wing_jy)] = gU[np.ix_(c.wing_ix,
                                                    c.wing_jy)] - gj
    elif how == 'decont':
        # local, bounded 2-point extrapolation of the UPPER branch to z=0
        # (node k0 is on the sheet side of the cut, k0+1 one cell further),
        # and the LOWER one-sided value from the exact sheet identity
        #   dphi/dx(0-) = dphi/dx(0+) - dJ/dx
        gU = upper_local(c)
        gj = gam_nodal(c)
        gL = np.array(gU)
        gL[np.ix_(c.wing_ix, c.wing_jy)] = gU[np.ix_(c.wing_ix,
                                                    c.wing_jy)] - gj
    elif how == 'node':
        gU, gL = c.gradx(k0), c.gradx(kl1)
    else:
        raise SystemExit(how)
    return c.cp_of_g(gU), c.cp_of_g(gL)


HOWS = ('fp', 'decont2', 'decont', 'node')
ORDER = ('nz12', 'nz16', 'nz24', 'nz32', 'zoff25', 'zoff50')


def main(argv):
    cases = {}
    for f in sorted(glob.glob(os.path.join(REPO, 'tmp_a17_case_*.npz'))):
        c = C(f)
        cases[c.tag] = c
    if argv and argv[0] in cases:                 # single-case dump
        c = cases[argv[0]]
        k0 = c.k0
        iu = np.asarray(c.wing_ix)
        j = c.wing_jy[len(c.wing_jy) // 2]
        xcf = c.xc[iu] / c.chord
        print('case %s  k0=%d z[k0]=%.5f  dz=%.5f/%.5f' % (
            c.tag, k0, c.z[k0], c.dz[k0 - 1], c.dz[k0]))
        print('jump at nodes  :', np.round(c.jump()[:, len(c.wing_jy) // 2], 5))
        print('dphi[k0]-[k0-1]:',
              np.round(c.phi[iu, j, k0] - c.phi[iu, j, k0 - 1], 5))
        for how in HOWS:
            cpu, cpl = est(c, how)
            dcp = (cpl[iu, j] - cpu[iu, j])
            print('%-8s CL=%.5f  dCp=%s' % (
                how, c.integ(cpu, cpl)[0], np.round(dcp, 4)))
        print('xc     :', np.round(xcf, 3))
        print('2*load :', np.round(2 * np.interp(
            np.clip(xcf, 0, 1), c.x_ctrl / c.chord,
            c.loading()[:, len(c.wing_jy) // 2]), 4))
        return
    rows = {}
    for tag in ORDER:
        c = cases[tag]
        rows[tag] = {h: float(c.integ(*est(c, h))[0]) for h in HOWS}
        rows[tag]['CL_fp_stored'] = c.CL_fp_stored
    print('\n%-10s' % 'CL' + ''.join('%9s' % t for t in ORDER))
    for h in HOWS:
        v = [rows[t][h] for t in ORDER]
        sp = 100 * (max(v[:4]) - min(v[:4])) / np.mean(v[:4])
        d1 = 100 * (v[4] - v[1]) / v[1]
        d2 = 100 * (v[5] - v[1]) / v[1]
        print('%-10s' % h + ''.join('%9.5f' % x for x in v) +
              '  dz=%5.1f%%  zoff=%.1f%%/%.1f%%' % (sp, d1, d2))
    with open(os.path.join(REPO, 'tmp_a17_lab.json'), 'w') as f:
        json.dump(rows, f, indent=1)


if __name__ == '__main__':
    main(sys.argv[1:])
