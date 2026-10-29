"""A15 multi-geometry x Mach TSD sweep with the VALIDATED configuration.

Validated config (do not tune here):
  WING3D_IMPRINT=ll   (Prandtl lifting-line imprint)
  WING3D_ENTROPY=1    (shock entropy correction)
  WING3D_KUTTA        unset/0  -- modes 6/7 are NOT validated for level
                                 or trend, so they are explicitly cleared
  grid g80x32x32, epochs 8, span 1.6, alpha 1.25 deg  (campaign protocol)

The point is to run *one* frozen configuration across many geometries and
report honestly where it fails.  No per-geometry tuning, ever.

Usage
  python sweep_multi.py                      # full 6 geo x 7 Mach set
  python sweep_multi.py --machs 0.8,0.925,1.0   # reduced set first
  python sweep_multi.py --report             # table + consistency checks
  python sweep_multi.py --grid 40x16x16      # LOW-RES fallback (labelled)

Every solve is checkpointed to --out (default tests3d/sweep_multi.json)
with an atomic replace, and an existing ok-record is skipped, so a restart
loses nothing.  RAM is checked before each solve; if the box is below
--min-free-mb the runner waits (up to --max-wait) instead of getting killed.
"""
import argparse
import ctypes
import json
import os
import sys
import time
import traceback

REPO = 'D:\\CodeProjects\\cfd'
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tests3d'))

import numpy as np                                            # noqa: E402

from gen_transonic_stl import (                               # noqa: E402
    naca0012_xy, naca2412_xy, naca4412_xy, rae2822_xy, biconvex_xy,
    dfvlr_r4_xy, section_metrics)

OUT_DEFAULT = os.path.join(REPO, 'tests3d', 'sweep_multi.json')

GEOS = ['naca0012', 'naca2412', 'naca4412', 'rae2822', 'biconvex',
        'dfvlr_r4']
MACHES = [0.800, 0.850, 0.900, 0.925, 0.950, 1.000, 1.050]
GRID = (80, 32, 32)
EPOCHS = 8
SPAN = 1.6
ALPHA = 1.25          # campaign alpha (matches the OF anchors)

# ---------------------------------------------------------------------------
# Anchors
# ---------------------------------------------------------------------------
# NEW free-tip OpenFOAM anchors (rhoCentralFoam, AR 1.6 free-tip), band
# +/-5.7%.  These SUPERSEDE the 86k corpus for dfvlr_r4 -- the corpus is
# 30-50% low there (0.2674 vs 0.32682 at M0.85).
NEW_FREETIP = {
    'dfvlr_r4': {'0.850': 0.32682, '0.925': 0.36580},
    'naca0012': {'0.800': 0.0571},
}
# M1.000 dfvlr free-tip value is a LOWER BOUND (CL still falling at the
# end of the OF run), so it may only be used as a one-sided check.
FREETIP_LOWER_BOUND = {'dfvlr_r4': {'1.000': 0.21032}}
FREETIP_BAND = 0.057          # +/-5.7%

CORPUS_CAVEAT = ('86k corpus value; the new free-tip anchors for this '
                 'geometry were NOT collected, so this anchor carries the '
                 'corpus bias (known 30-50% low for dfvlr) -- treat the '
                 'error as indicative only')


def load_corpus():
    p = os.path.join(REPO, 'tests3d', 'validation', 'known_anchors.json')
    with open(p) as f:
        j = json.load(f)
    return j.get('corpus_86k', {})


# ---------------------------------------------------------------------------
# geometry -> TSD profile / section
# ---------------------------------------------------------------------------
XY = {'naca0012': lambda: _sym(naca0012_xy()),
      'naca2412': naca2412_xy,
      'naca4412': naca4412_xy,
      'rae2822': rae2822_xy,
      'biconvex': lambda: _sym(biconvex_xy()),
      'dfvlr_r4': dfvlr_r4_xy}


def _sym(xy):
    """(x, yt) half-thickness form -> (xu, yu, xl, yl)."""
    x, y = xy
    return x, y, x, -y


def geometry_xy(geo):
    return XY[geo]()


def build_profile(geo):
    """Slope/thickness interpolants for the TSD forcing (WING3D profile)."""
    from wing3d.tsd import section_profile_xy
    from gen_transonic_stl import _monotone_x
    xu, yu, xl, yl = geometry_xy(geo)
    xu, yu = _monotone_x(xu, yu)
    xl, yl = _monotone_x(xl, yl)
    xg = np.linspace(0, 1, 400)
    return section_profile_xy(xg, np.interp(xg, xu, yu),
                              np.interp(xg, xl, yl))


def build_section(geo, n=24):
    from wing3d.geometry import section_loop
    from gen_transonic_stl import _monotone_x
    xu, yu, xl, yl = geometry_xy(geo)
    xu, yu = _monotone_x(xu, yu)
    xl, yl = _monotone_x(xl, yl)
    return section_loop(xu, yu, xl, yl, n)


# ---------------------------------------------------------------------------
# RAM gate
# ---------------------------------------------------------------------------
class _MEMSTATUS(ctypes.Structure):
    _fields_ = [('dwLength', ctypes.c_ulong),
                ('dwMemoryLoad', ctypes.c_ulong),
                ('ullTotalPhys', ctypes.c_ulonglong),
                ('ullAvailPhys', ctypes.c_ulonglong),
                ('ullTotalPageFile', ctypes.c_ulonglong),
                ('ullAvailPageFile', ctypes.c_ulonglong),
                ('ullTotalVirtual', ctypes.c_ulonglong),
                ('ullAvailVirtual', ctypes.c_ulonglong),
                ('ullAvailExtendedVirtual', ctypes.c_ulonglong)]


def free_mb():
    ms = _MEMSTATUS()
    ms.dwLength = ctypes.sizeof(_MEMSTATUS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms))
    return ms.ullAvailPhys / (1024.0 * 1024.0)


def wait_for_ram(min_free_mb, max_wait_s, label=''):
    t0 = time.time()
    while True:
        f = free_mb()
        if f >= min_free_mb:
            return f
        if time.time() - t0 > max_wait_s:
            print('  RAM %.0f MB < %d MB after %.0fs -- proceeding anyway'
                  % (f, min_free_mb, max_wait_s), flush=True)
            return f
        print('  RAM %.0f MB < %d MB, waiting 60s (%s)'
              % (f, min_free_mb, label), flush=True)
        time.sleep(60)


# ---------------------------------------------------------------------------
# checkpointing
# ---------------------------------------------------------------------------
def key(geo, M):
    return '%s@M%.3f' % (geo, M)


def load_db(path):
    if os.path.exists(path):
        try:
            with open(path) as f:
                return json.load(f)
        except Exception:
            print('WARNING: %s unreadable, starting a fresh registry' % path)
    return {'_meta': {}, 'points': {}}


def save_db(path, db):
    db['_meta']['updated'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    db['_meta']['n_points'] = len(db['points'])
    db['_meta']['n_ok'] = sum(
        1 for r in db['points'].values() if r.get('status') == 'ok')
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(db, f, indent=1)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# one solve
# ---------------------------------------------------------------------------
def apply_validated_env():
    """Frozen validated configuration -- set explicitly, never inherited."""
    os.environ['WING3D_IMPRINT'] = 'll'
    os.environ['WING3D_ENTROPY'] = '1'
    for k in ('WING3D_KUTTA', 'WING3D_BL', 'WING3D_MORINO',
              'WING3D_KUTTA_THICK', 'WING3D_BL_SCALE', 'WING3D_ALPHA'):
        os.environ.pop(k, None)
    os.environ['WING3D_KUTTA'] = '0'
    os.environ['WING3D_BL'] = '0'


def run_one(M, geo, grid, epochs, base=None, note=''):
    """One (geo, M) solve.  Returns (record, base)."""
    from tsd_coupled import build_case
    from tsd_m_sup import run_sup

    prof = build_profile(geo)
    sec = build_section(geo)
    nx, ny, nz = grid
    t0 = time.time()
    built = base is None
    if built:
        base, _r, _m = build_case(0.3, ALPHA, False, True, None, nx, ny, nz,
                                  SPAN, section=sec, profile=prof)
    hist = []

    def cb(ep, cl, nsup):
        hist.append({'ep': int(ep), 'CL': float(cl), 'nsup': int(nsup)})

    L, xc, ix, j0 = run_sup(M, span=SPAN, epochs=epochs, free_k=False,
                            profile=prof, section=sec, base=base,
                            callback=cb, alpha_deg=ALPHA)
    rec = {
        'geo': geo, 'M': M, 'status': 'ok',
        'grid': '%dx%dx%d' % (nx, ny, nz),
        'epochs': epochs, 'span': SPAN, 'alpha_deg': ALPHA,
        'config': {'imprint': 'll', 'entropy': '1', 'kutta': '0',
                   'bl': '0', 'free_k': False},
        'CL': float(L['CL']), 'CD': float(L['CD']),
        'hist': hist,
        'cl_ep_last': hist[-1]['CL'] if hist else None,
        'cl_ep_drift': (hist[-1]['CL'] - hist[0]['CL']) if len(hist) > 1
        else 0.0,
        'nsup_last': hist[-1]['nsup'] if hist else None,
        'built_base': built,
        'secs': time.time() - t0, 'note': note,
        'ram_free_mb_at_start': round(free_mb()),
    }
    if (nx, ny, nz) != GRID:
        rec['resolution'] = 'LOW-RES %s (absolute CL NOT comparable to the ' \
            'g80x32x32 anchors; trends only)' % rec['grid']
    else:
        rec['resolution'] = 'campaign g80x32x32'
    return rec, base


def grade(rec, corpus):
    """Attach anchors.  Free-tip anchors are authoritative where they exist."""
    geo, M = rec['geo'], rec['M']
    k = '%.3f' % M
    rec.pop('anchor', None)
    rec.pop('err_pct', None)
    rec.pop('anchor_src', None)
    rec.pop('anchor_note', None)
    rec.pop('pass_band', None)
    cl = rec.get('CL')
    fa = NEW_FREETIP.get(geo, {})
    lb = FREETIP_LOWER_BOUND.get(geo, {})
    if k in fa:
        a = fa[k]
        rec.update(anchor=a, anchor_src='freetip_rhoCentralFoam',
                   anchor_note='AR1.6 free-tip; band +/-5.7%%')
        if cl is not None:
            rec['err_pct'] = 100.0 * (cl - a) / a
            rec['pass_band'] = bool(abs(rec['err_pct']) <= 100 * FREETIP_BAND)
        return rec
    if k in lb:
        a = lb[k]
        rec.update(anchor=a, anchor_src='freetip_rhoCentralFoam_LOWER_BOUND',
                   anchor_note='OF CL still falling at end of run -- true '
                               'value is >= this anchor, so a TSD CL above '
                               'it is NOT automatically wrong')
        if cl is not None:
            rec['err_pct'] = 100.0 * (cl - a) / a
            rec['pass_band'] = None
        return rec
    c = corpus.get(geo, {}).get(k) or corpus.get(geo, {}).get(
        ('%.1f' % M))
    if c is not None and rec['grid'] == '80x32x32':
        rec.update(anchor=c, anchor_src='corpus_86k',
                   anchor_note=CORPUS_CAVEAT, err_pct=(
                       100.0 * (cl - c) / c if cl is not None else None))
        rec['pass_band'] = None
    return rec


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def report(path, min_free_mb):
    db = load_db(path)
    pts = db['points']
    corpus = load_corpus()
    geos = [g for g in GEOS if any(p['geo'] == g for p in pts.values())]
    machs = sorted({p['M'] for p in pts.values()})
    grid = sorted({p['grid'] for p in pts.values()})
    print('\n=== A15 sweep report =======================================')
    print('file: %s' % path)
    print('grid(s): %s   epochs: %s   alpha: %.2f deg   span: %.1f'
          % (', '.join(grid),
             sorted({p['epochs'] for p in pts.values()}), ALPHA, SPAN))
    if any(g != '80x32x32' for g in grid):
        print('*** LOW-RES RUN: absolute CL is NOT comparable to the '
              'g80x32x32 anchors -- trends only ***')

    def cell(g, M):
        p = pts.get(key(g, M))
        if p is None:
            return '   --   '
        if p.get('status') != 'ok':
            return '  FAIL  '
        return '%8.4f' % p['CL']

    print('\nCL table (geometry x Mach), validated config ll+entropy')
    print('%-10s' % 'geo' + ''.join('%9.3f' % M for M in machs))
    for g in geos:
        print('%-10s' % g + ''.join('%9s' % cell(g, M) for M in machs))

    print('\nCD table')
    for g in geos:
        row = ''
        for M in machs:
            p = pts.get(key(g, M))
            row += ('%9s' % ('%8.4f' % p['CD']
                             if p and p.get('status') == 'ok' else '--'))
        print('%-10s' % g + row)

    print('\nCamber monotonicity (CL should rise 0012 < 2412 < 4412 at '
          'fixed M; t/c is 0.120 for all three)')
    ok_all = True
    for M in machs:
        vals = [pts.get(key(g, M), {}) for g in
                ('naca0012', 'naca2412', 'naca4412')]
        if any(v.get('status') != 'ok' for v in vals):
            print('  M%.3f: incomplete' % M)
            continue
        c = [v['CL'] for v in vals]
        mono = c[0] < c[1] < c[2]
        ok_all &= mono
        print('  M%.3f: 0012 %.4f | 2412 %.4f (%+.1f%%) | 4412 %.4f '
              '(%+.1f%%)  -> %s'
              % (M, c[0], c[1], 100 * (c[1] - c[0]) / c[0], c[2],
                 100 * (c[2] - c[0]) / c[0],
                 'MONOTONE' if mono else '*** VIOLATED ***'))
    print('  verdict: %s' % ('monotone at every Mach' if ok_all
                             else 'NOT monotone at >=1 Mach'))

    print('\nMach trend: dCL/dM over the sampled band, collapse ratio '
          'CL(1.000)/CL(0.850)')
    print('%-10s %10s %10s %12s %s'
          % ('geo', 'dCL/dM', 'CL0.85', 'ratio1.0/0.85', 'shape'))
    for g in geos:
        ms = [M for M in machs if pts.get(key(g, M), {}).get('status') == 'ok']
        if len(ms) < 2:
            print('%-10s  incomplete' % g)
            continue
        cs = [pts[key(g, M)]['CL'] for M in ms]
        sl = float(np.polyfit(ms, cs, 1)[0])
        r = (pts[key(g, 1.000)]['CL'] / pts[key(g, 0.850)]['CL']
             if key(g, 1.000) in pts and key(g, 0.850) in pts
             and pts[key(g, 0.850)]['CL'] else float('nan'))
        imax = int(np.argmax(cs))
        shape = ('peak at M%.3f then falls' % ms[imax]
                 if imax < len(cs) - 1 else 'rising to top of band')
        print('%-10s %+10.4f %10.4f %12s %s'
              % (g, sl, cs[ms.index(0.850)] if 0.850 in ms else float('nan'),
                 ('%.3f' % r) if np.isfinite(r) else '--', shape))

    print('\nAnchored errors')
    hdr = '%-10s %6s %9s %9s %9s %8s  %s' % (
        'geo', 'M', 'CL', 'anchor', 'err%', 'in5.7%', 'anchor source')
    print(hdr)
    graded = []
    for g in geos:
        for M in machs:
            p = grade(dict(pts.get(key(g, M), {})) if key(g, M) in pts
                      else {'geo': g, 'M': M, 'grid': 'x'}, corpus)
            if 'anchor' not in p or p.get('CL') is None:
                continue
            graded.append(p)
            print('%-10s %6.3f %9.4f %9.5f %+8.1f %8s  %s'
                  % (g, M, p['CL'], p['anchor'], p['err_pct'],
                     {True: 'PASS', False: 'FAIL', None: 'n/a'}[
                         p.get('pass_band')],
                     p.get('anchor_src')))
    print('\n  notes: dfvlr M1.000 free-tip anchor is a LOWER BOUND '
          '(OF CL still falling) -> err%% is an upper bound on the '
          'discrepancy; rae2822/naca0012-at-0.85+ rows use the 86k corpus '
          'and carry its bias.')
    fails = [p for p in graded if p.get('pass_band') is False]
    print('  band verdict: %d/%d graded points inside +/-5.7%%'
          % (len([p for p in graded if p.get('pass_band')]), len(graded)))

    bad = {k: v for k, v in pts.items() if v.get('status') != 'ok'}
    if bad:
        print('\nFailed points (%d)' % len(bad))
        for k, v in sorted(bad.items()):
            print('  %-18s %s: %s' % (k, v['status'],
                                      (v.get('error') or '')[:160]))
    return fails


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--geos', default=','.join(GEOS))
    ap.add_argument('--machs', default=','.join('%.3f' % m for m in MACHES))
    ap.add_argument('--grid', default='%dx%dx%d' % GRID)
    ap.add_argument('--epochs', type=int, default=EPOCHS)
    ap.add_argument('--out', default=OUT_DEFAULT)
    ap.add_argument('--min-free-mb', type=int, default=700)
    ap.add_argument('--max-wait', type=int, default=3600)
    ap.add_argument('--force', action='store_true',
                    help='re-run points already recorded ok')
    ap.add_argument('--report', action='store_true')
    ap.add_argument('--report-only', action='store_true')
    a = ap.parse_args(argv)

    if a.report_only:
        report(a.out, a.min_free_mb)
        return 0

    geos = [g.strip() for g in a.geos.split(',') if g.strip()]
    machs = [float(m) for m in a.machs.split(',') if m.strip()]
    grid = tuple(int(v) for v in a.grid.lower().split('x'))
    for g in geos:
        if g not in XY:
            raise SystemExit('unknown geometry %r (have %s)'
                             % (g, ', '.join(XY)))
    if grid != GRID:
        print('!! LOW-RES GRID %s: absolute CL will NOT be comparable to '
              'the g80x32x32 anchors; label results accordingly' % (grid,))

    apply_validated_env()
    db = load_db(a.out)
    db['_meta'].update({
        'agent': 'A15 multi-geometry sweep',
        'config': 'WING3D_IMPRINT=ll WING3D_ENTROPY=1 WING3D_KUTTA=0 '
                  'WING3D_BL=0',
        'protocol': 'grid %s, epochs %d, span %.1f, alpha %.2f deg'
                    % (grid, a.epochs, SPAN, ALPHA),
        'caveat': 'Unanchored geometries (naca2412/naca4412) report raw CL '
                  'only -- no anchor exists.  Anchors grade against the NEW '
                  'free-tip rhoCentralFoam values where they exist, else '
                  'the 86k corpus with its known bias.',
    })
    save_db(a.out, db)

    todo = [(g, M) for g in geos for M in machs
            if a.force or db['points'].get(key(g, M), {}).get(
                'status') != 'ok']
    print('A15 sweep: %d geometry x %d Mach, %d to run (grid %s, ep %d)'
          % (len(geos), len(machs), len(todo), grid, a.epochs), flush=True)

    t_all = time.time()
    base = None
    base_geo = None
    for i, (g, M) in enumerate(todo, 1):
        print('\n[%d/%d] %s M%.3f  (RAM %.0f MB free)'
              % (i, len(todo), g, M, free_mb()), flush=True)
        wait_for_ram(a.min_free_mb, a.max_wait, '%s M%.3f' % (g, M))
        if base_geo != g:
            base = None          # never carry one geometry's state to another
            base_geo = g
        try:
            rec, base = run_one(M, g, grid, a.epochs, base=base)
        except MemoryError as e:
            rec = {'geo': g, 'M': M, 'status': 'memkill', 'error': str(e),
                   'grid': '%dx%dx%d' % grid, 'epochs': a.epochs}
            base = None
            print('  MEMKILL: %s' % e, flush=True)
            time.sleep(60)
        except Exception as e:
            rec = {'geo': g, 'M': M, 'status': 'error',
                   'error': '%s: %s' % (type(e).__name__, e),
                   'traceback': traceback.format_exc()[-1500:],
                   'grid': '%dx%dx%d' % grid, 'epochs': a.epochs}
            base = None          # a failed solve may have poisoned the state
            print('  FAIL: %s' % e, flush=True)
        rec = grade(rec, load_corpus())
        db['points'][key(g, M)] = rec
        save_db(a.out, db)                       # checkpoint every solve
        if rec.get('status') == 'ok':
            extra = ''
            if 'err_pct' in rec:
                extra = '  err=%+.1f%% vs %s' % (rec['err_pct'],
                                                 rec['anchor_src'])
            print('  CL=%.5f CD=%.5f  drift(ep)=%+.4f  %.0fs%s -> %s'
                  % (rec['CL'], rec['CD'], rec['cl_ep_drift'], rec['secs'],
                     extra, os.path.basename(a.out)), flush=True)
        # release the prebuilt base only when the geometry changes
    print('\nall done in %.1f min; %d record(s) in %s'
          % ((time.time() - t_all) / 60.0, len(db['points']), a.out),
          flush=True)
    report(a.out, a.min_free_mb)
    return 0


if __name__ == '__main__':
    sys.exit(main())
