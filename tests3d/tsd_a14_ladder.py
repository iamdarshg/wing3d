"""A14 ladder: BL symmetry, B' = B + w_thick (circulation-off reference),
and the dz_core/nx scaling of the discrete offset.

  python tsd_a14_ladder.py BL     dfvlr M0.85 +/- WING3D_BL (test C redo)
  python tsd_a14_ladder.py THK    M0.3 symmetric, mode 7, WING3D_KUTTA_THICK
                                  on/off  (does the -26% level restore?)
  python tsd_a14_ladder.py B      transonic dfvlr 0.85/0.925/1.00, B'
  python tsd_a14_ladder.py GRID   M0.3 symmetric grid sweep (dz_core, nx)

Recorder `one14` mirrors `one7` (tsd_a13_ladder) and additionally stores
w_thick (the reference field mapped by L), the post-reference residual and
the BL wav decomposition, so the B' effect is readable from the record.
"""
import json
import os
import sys
import time

import numpy as np

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from tsd_a11_ladder import ANCHORS, OUT, tag              # noqa: E402
from tsd_coupled import build_case                       # noqa: E402
from tsd_m_sup import run_sup                            # noqa: E402
from tsd_campaign import build_profile, build_section    # noqa: E402


def one14(M, geo='dfvlr_r4', epochs=2, span=1.6, use_profile=True,
          grid=(80, 32, 32), note='', prefix='a14', alpha=1.25):
    prof = build_profile(geo) if use_profile else None
    sec = build_section(geo)
    nx, ny, nz = grid
    t0 = time.time()
    base, _r, _m = build_case(0.3, alpha, False, True, None, nx, ny, nz,
                              span, section=sec, profile=prof)
    L, xc, ix, j0 = run_sup(M, span=span, epochs=epochs, free_k=False,
                            profile=prof, section=sec, base=base)
    rec = {'geo': geo, 'M': M, 'grid': '%dx%dx%d' % (nx, ny, nz),
           'epochs': epochs, 'span': span,
           'mode': os.environ.get('WING3D_KUTTA', '0'),
           'entropy': os.environ.get('WING3D_ENTROPY', '0'),
           'bl': os.environ.get('WING3D_BL', '0'),
           'thick': os.environ.get('WING3D_KUTTA_THICK', '0'),
           'thick_ep': os.environ.get('WING3D_THICK_EP', '1'),
           'nc': os.environ.get('WING3D_KUTTA_NC', 'all'),
           'alpha_deg': alpha, 'free_k': False, 'profile': use_profile,
           'CL': float(L['CL']), 'CD': float(L['CD']),
           'secs': time.time() - t0, 'note': note}
    a = ANCHORS.get('%.3f' % M)
    if a:
        rec['anchor'] = a
        rec['anchor_old'] = a
        rec['err_pct'] = 100.0 * (rec['CL'] - a) / a
    K = L['kutta']
    rec['w'] = K.get('w')
    rec['hist'] = K['hist']
    if 'grid' in K:
        rec['gridinfo'] = K['grid']
    y = np.asarray(K['y'])
    j0 = int(np.argmin(np.abs(y)))
    mu = np.asarray(K['mu'])
    Gsec = np.asarray(K['Gamma_mu'])
    load = np.asarray(K['loading'])
    rec['x_ctrl'] = np.asarray(K['x_ctrl']).tolist()
    rec['mu_mid'] = mu[:, j0].tolist()
    rec['loading_mid'] = load[:, j0].tolist()
    rec['Gamma_mu'] = Gsec.tolist()
    rec['Gamma_mu_mid'] = float(Gsec[j0])
    rec['tang_rows_mid'] = np.asarray(K['tang_rows'])[:, j0].tolist()
    rec['loading_TE_ratio'] = float(
        abs(load[-1, j0]) / max(abs(load[:, j0]).max(), 1e-12))
    G = np.asarray(K['Gamma'])
    rec['Gamma'] = G.tolist()
    if len(G) > 2:
        rec['gamma_mid'] = float(G[len(G) // 2])
        rec['gamma_tip_ratio'] = float(
            np.mean([abs(G[0]), abs(G[-1])]) / max(abs(G[len(G) // 2]), 1e-12))
    al = np.radians(alpha)
    if 'w_field' in K:
        wf = np.asarray(K['w_field'])
        we = np.asarray(K['w_eq'])
        bi = np.asarray(K['bias'])
        rec['w_field'] = wf.tolist()
        rec['w_eq'] = we.tolist()
        rec['bias'] = bi.tolist()
        rec['resid_phys'] = (wf - we).tolist()
        rec['resid_phys_mid_alpha'] = ((wf[:, j0] - we[:, j0]) / al).tolist()
        # A13's algebraic reading: a uniform residual offset r (in alpha
        # units) is an extra downwash, so mu is scaled by 1/(1 - r). Record
        # it so the grid study can track the "1.36x" factor directly.
        _r = np.asarray(rec['resid_phys_mid_alpha'])
        rec['resid_phys_mean_alpha'] = float(_r.mean())
        rec['norm_factor'] = float(1.0 / max(1.0 - _r.mean(), 1e-6))
        rec['coupling'] = K['coupling']
        if 'w_thick' in K:
            wt = np.asarray(K['w_thick'])
            rec['w_thick'] = wt.tolist()
            rec['w_thick_alpha_mid'] = (wt[:, j0] / al).tolist()
            # A14 residual: field - analytic sheet - reference field
            rec['resid_phys2'] = (wf - we - np.asarray(K['w_thick_node'])
                                  ).tolist()
            rec['resid_phys2_mid_alpha'] = (
                (wf[:, j0] - we[:, j0]
                 - np.asarray(K['w_thick_node'])[:, j0]) / al).tolist()
            _r2 = np.asarray(rec['resid_phys2_mid_alpha'])
            rec['resid_phys2_mean_alpha'] = float(_r2.mean())
    if 'bl' in L:
        rec['bl_hist'] = [{'dstar_max': h['dstar_max'],
                           'change': h['change'],
                           'n_sep': h['n_sep']} for h in L['bl']['hist']]
        rec['bl_dstar_max'] = float(max(np.max(L['bl']['dstar_u']),
                                        np.max(L['bl']['dstar_l'])))
        if 'wav_terms' in L['bl']:
            rec['bl_wav_terms'] = L['bl']['wav_terms']
    suff = ''
    if grid != (80, 32, 32):
        suff = '_g%dx%dx%d' % grid
    if os.environ.get('WING3D_KUTTA_THICK', '0') == '1':
        suff += '_thk'
    if os.environ.get('WING3D_BL_SCALE', '1') not in ('1', ''):
        suff += '_bls' + os.environ['WING3D_BL_SCALE']
    name = '%s_%s_M%.3f_%s_ep%d%s%s.json' % (
        prefix, geo, M, tag(), epochs, suff,
        '' if use_profile else '_noprof')
    path = os.path.join(OUT, name)
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1)
    extra = ''
    if 'err_pct' in rec:
        extra = '  err=%+.2f%% vs %.5f' % (rec['err_pct'], rec['anchor'])
    if 'gamma_mid' in rec:
        extra += '  G_mid=%.4f  tip/mid=%.2f' % (
            rec['gamma_mid'], rec['gamma_tip_ratio'])
    if 'resid_phys_mid_alpha' in rec:
        r1 = rec['resid_phys_mid_alpha']
        extra += '  res1/alpha[%.2f..%.2f]' % (min(r1), max(r1))
    if 'resid_phys2_mid_alpha' in rec:
        r2 = rec['resid_phys2_mid_alpha']
        extra += '  res2/alpha[%.2f..%.2f]' % (min(r2), max(r2))
    if 'bl_wav_terms' in rec:
        extra += '  BL wav/alpha=%.4f (source/2=%.3e)' % (
            rec['bl_wav_terms']['wav_alpha_absmax'],
            rec['bl_wav_terms']['half_absmax'])
    print('RESULT %s  M=%.3f CL=%.5f%s  (%s) -> %s' % (
        tag(), M, rec['CL'], extra, '%.0fs' % rec['secs'], name), flush=True)
    return rec


def main(ladder):
    ep = int(os.environ.get('A14_EP', '2'))
    epa = int(os.environ.get('A14_EPA', '4'))
    if ladder == 'BL':
        for M in (0.85, 0.3):
            one14(M, epochs=ep, note='BL symmetry M%.2f' % M)
    elif ladder == 'BLS':
        # liveness/ablation of the viscous term actually reaching wav:
        # WING3D_BL_SCALE multiplies the displacement slope, so the sign of
        # dCL/d(displacement) is measurable even though the physical term is
        # small.  SCALE=0 is bit-identical to WING3D_BL=0.
        for M, sc in ((0.3, '1'), (0.3, '50'), (0.85, '50')):
            os.environ['WING3D_BL_SCALE'] = sc
            one14(M, epochs=ep, use_profile=False, note='BL scale %s' % sc)
    elif ladder == 'THK':
        one14(0.3, epochs=epa, use_profile=False, note='B-prime M0.3 sym')
    elif ladder == 'B':
        for M in (0.85, 0.925, 1.0):
            one14(M, epochs=ep, note="B' transonic")
    elif ladder == 'GRID':
        for g in ((80, 32, 32), (120, 40, 40), (120, 48, 48),
                  (160, 56, 56)):
            one14(0.3, epochs=ep, use_profile=False, grid=g,
                  note='grid %s' % (g,))
    else:
        raise SystemExit('unknown ladder %s' % ladder)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'THK')
