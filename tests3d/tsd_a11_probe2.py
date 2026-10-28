"""A11 probe 2: LOCAL structure of the mode-5 tangency operator.

probe 1 (tsd_a11_probe.py) scaled the whole mu matrix and measured a
uniform d wav / d mu ~ -1.05. That says the rows respond to the sheet
AMPLITUDE, but not whether they can resolve the chordwise distribution.
Here a single chordwise panel is perturbed (all span cells together) and
the response is read at every chordwise station, at mid span:

  * a DELTA response (only stations next to the bumped panel move) means
    the row set is a difference/dipole operator -> nc rows resolve shape
  * a FLAT response (~ same everywhere) means every row reads the same
    global quantity -> nc rows are rank-deficient in the chordwise
    direction and the shape is NOT resolvable (the honest failure mode)

Also reports the response to a CONSTANT shift of mu (add delta to every
panel): a continuum doublet sheet has loading gamma = dmu/dx, so a
constant shift must produce (almost) no downwash -- the size of that
response measures how much of the code's sheet is "mu itself" rather
than "dmu/dx".

  python tsd_a11_probe2.py sym|cam
"""
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import newton_krylov

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

from tsd_a11_probe import setup, phi_solve              # noqa: E402
from tsd_campaign import build_profile                  # noqa: E402
from tsd_coupled import lifting_line_gamma              # noqa: E402
from wing3d.tsd import flat_plate_jump                  # noqa: E402


def probe(kind='sym'):
    prof = None if kind == 'sym' else build_profile('dfvlr_r4')
    t0 = time.time()
    s = setup(0.3, 1.6, profile=prof)
    g_ll = np.asarray(lifting_line_gamma(
        1.25, 0.3, s.span, chord=s.chord, profile=prof,
        y=s.yc[np.asarray(s.wing_jy)]), dtype=float)
    mu0 = flat_plate_jump(s.x_ctrl)[:, None] * g_ll[None, :]
    nc, ng = s.n_ctrl, len(s.wing_jy)
    phi_solve(s, mu0, ep=2)
    phi_base = s.phi.copy()
    _wu, _wl, wav = s.w_plane()
    iu = np.asarray(s.wing_ix)[s.i_ctrl][:, None]
    jj = np.asarray(s.wing_jy)[None, :]
    base = wav[iu, jj]
    print('base solved (%.0fs)  wav(mid, 0.75c)=%.5f' % (
        time.time() - t0, base[9, nc // 2]), flush=True)

    delta = 0.004
    cases = {'panel2': 2, 'panel6': 6, 'panel11': 11, 'const': -1}
    out = {'x_ctrl': s.x_ctrl.tolist(),
           'y': s.yc[np.asarray(s.wing_jy)].tolist(),
           'base': base.tolist(), 'delta': delta, 'nc': nc}
    for name, ip in cases.items():
        mu = mu0.copy()
        if ip < 0:
            mu += delta
        else:
            mu[ip, :] += delta
        s.phi[:] = phi_base
        t = time.time()
        phi_solve(s, mu, ep=1)
        _wu, _wl, wav = s.w_plane()
        resp = (wav[iu, jj] - base) / delta
        out[name] = resp.tolist()
        print('  bump %-8s (%.0fs)  response at mid span, by x_ctrl:' % (
            name, time.time() - t), flush=True)
        print('    ' + ' '.join('%7.2f' % resp[i, nc // 2]
                                for i in range(nc)), flush=True)
        print('    spread over span (min/max of column med): %.2f / %.2f' % (
            np.median(resp, axis=1).min(), np.median(resp, axis=1).max()),
            flush=True)
    out['secs'] = time.time() - t0
    path = os.path.join(REPO, 'tests3d', 'a11_probe2_%s.json' % kind)
    with open(path, 'w') as f:
        json.dump(out, f, indent=1)
    print('  -> %s (%.0fs)' % (path, out['secs']), flush=True)


if __name__ == '__main__':
    probe(sys.argv[1] if len(sys.argv) > 1 else 'sym')
