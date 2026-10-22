"""ctest-style validation harness: TSD solver vs published NASA/AGARD refs.

Suites:
  flow      -- run the TSD (run_sup, g80x32x32 ep8 by default, AR6 span)
               at each reference (M, alpha) and grade CL vs cl_pub/cl_band.
               Also grades CL-ratio pairs (cancels finite-AR scale).
  geometry  -- instant data checks: t/c, camber, TE gap of our generated
               sections vs the published airfoil definitions.

Usage:
  python tests3d/validation/ctest_validate.py [--suite flow|geometry|all]
       [--targets k1,k2] [--label NAME] [--ep N] [--grid 80 32 32]
       [--span 6.0]

Env:
  WING3D_ALPHA    default 1.25 (per-target alpha wins)
  WING3D_FREE_K   '1' -> free-k closure (default fixed k=1)
  WING3D_BC_ROWS  '1' -> real farfield rows (default rank-deficient zeros)
  WING3D_TE_GAP   geometry knob (default 0.002)

Exit code 0 = all PASS, 1 = any FAIL (for agent gating).
Results JSON: tests3d/validation/results/<label>.json

RAM: flow targets wait for >=2.0 GB free (they need ~1 GB each) up to
30 min, then FAIL-skips with reason=ram_gated.
"""
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'tests3d'))

REFS = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'published_refs.json')))


def ram_free_gb():
    try:
        out = subprocess.run(
            ['powershell', '-NoProfile', '-Command',
             '[math]::Round((Get-CimInstance Win32_OperatingSystem).'
             'FreePhysicalMemory/1MB,2)'],
            capture_output=True, text=True, timeout=60)
        return float(out.stdout.strip())
    except Exception:
        return 0.0


def wait_ram(min_gb=2.0, timeout_s=1800):
    """Block until >= min_gb free. Returns True if acquired."""
    t0 = time.time()
    while time.time() - t0 < timeout_s:
        if ram_free_gb() >= min_gb:
            return True
        time.sleep(60)
    return False


def run_target(t, grid, ep, span, label, results, free_k, bc_rows):
    key = t['key']
    if t.get('cl_pub') is None or t.get('status', '') not in ('', None):
        results.append({'key': key, 'status': 'PENDING',
                        'reason': t.get('status', 'cl_pub missing')})
        print('  %-24s -> PENDING (%s)' % (key, t.get('status')),
              flush=True)
        return None
    from tsd_campaign import build_profile, build_section
    from tsd_m_sup import run_sup
    geo = t['geometry']
    M, alpha = float(t['M']), float(t['alpha'])
    rec = {'key': key, 'geometry': geo, 'M': M, 'alpha': alpha,
           'cl_pub': t['cl_pub'], 'cl_band': t['cl_band']}
    if not wait_ram():
        rec.update({'status': 'FAIL', 'reason': 'ram_gated',
                    'cl_tsd': None, 'err_pct': None})
        results.append(rec)
        print('RAM-GATED skip %s' % key, flush=True)
        return None
    try:
        prof = build_profile(geo)
        sec = build_section(geo)
        L, xc, ix, j0 = run_sup(M, span=span, epochs=ep, free_k=free_k,
                                profile=prof, section=sec, alpha_deg=alpha)
        cl = float(L['CL'])
        err = (cl - t['cl_pub']) / t['cl_pub'] * 100
        inside = abs(cl - t['cl_pub']) <= t['cl_band']
        rec.update({'status': 'PASS' if inside else 'FAIL',
                    'cl_tsd': round(cl, 6), 'err_pct': round(err, 2),
                    'band': t['cl_band'], 'inside': inside,
                    'nsup': int(float(L.get('nsup', 0) or 0))})
    except Exception as e:
        rec.update({'status': 'FAIL', 'reason': 'exception: %s' % e,
                    'cl_tsd': None, 'err_pct': None})
    results.append(rec)
    print('  %-24s CL=%.5f pub=%.4f err=%+.1f%% -> %s' % (
        key, rec.get('cl_tsd') or 0, t['cl_pub'],
        rec.get('err_pct') or 0, rec['status']), flush=True)
    return rec


def grade_ratios(results):
    bykey = {r['key']: r for r in results if r.get('cl_tsd') is not None}
    for t in REFS['flow']:
        prt = t.get('cl_ratio_partner')
        if not prt or prt not in bykey or t['key'] not in bykey:
            continue
        r_hi = bykey[t['key']]; r_lo = bykey[prt]
        if not r_hi['cl_tsd'] or not r_lo['cl_tsd'] or r_lo['cl_tsd'] == 0:
            continue
        ratio = r_hi['cl_tsd'] / r_lo['cl_tsd']
        pub = t['cl_ratio_pub']
        rec = {'key': 'ratio_%s_over_%s' % (t['key'], prt),
               'cl_ratio': round(ratio, 4), 'cl_ratio_pub': pub,
               'err_pct': round((ratio - pub) / pub * 100, 2),
               'status': 'PASS' if abs(ratio - pub) / pub < 0.15
               else 'FAIL'}
        results.append(rec)
        print('  %-24s ratio=%.3f pub=%.3f err=%+.1f%% -> %s' % (
            rec['key'], ratio, pub, rec['err_pct'], rec['status']),
            flush=True)


def geom_measure(geo):
    """t/c, max camber + its x, TE gap from our generated section.

    Coordinates from gen_transonic_stl are chord-normalized (x,y in [0,1]
    fractions of chord). Reverses descending arrays (TE->LE) before
    interpolation, like build_profile does.
    """
    from gen_transonic_stl import (naca0012_xy, rae2822_xy,
                                   biconvex_xy, dfvlr_r4_xy)
    import numpy as np
    if geo == 'naca0012':
        xa, yt = naca0012_xy(400); xu, yu, xl, yl = xa, yt, xa, -yt
    elif geo == 'biconvex':
        xb, yb = biconvex_xy(400); xu, yu, xl, yl = xb, yb, xb, -yb
    elif geo == 'rae2822' or geo == 'dfvlr_r4':
        xu, yu, xl, yl = rae2822_xy() if geo == 'rae2822' \
            else dfvlr_r4_xy()
    else:
        raise ValueError(geo)
    xu, yu, xl, yl = map(np.asarray, (xu, yu, xl, yl))
    if xl[0] > xl[-1]:
        xl, yl = xl[::-1], yl[::-1]
    if xu[0] > xu[-1]:
        xu, yu = xu[::-1], yu[::-1]
    xg = np.linspace(0, 1, 400)
    yu_i = np.interp(xg, xu, yu)
    yl_i = np.interp(xg, xl, yl)
    th = yu_i - yl_i
    tc = float(np.max(th))
    cam = 0.5 * (yu_i + yl_i)
    camber = float(np.max(np.abs(cam)))
    camber_x = float(xg[np.argmax(np.abs(cam))])
    te_gap = float(np.abs(yu_i[-1] - yl_i[-1]))
    return tc, camber, camber_x, te_gap


def run_geometry(results):
    for g in REFS['geometry']:
        if g['geo'] == 'onera_m6':
            results.append({'key': g['key'], 'status': 'PENDING',
                            'reason': 'M6 geometry lives in OF case '
                                      'generator, not TSD proxies'})
            print('  %-20s -> PENDING (OF-side geometry)' % g['key'],
                  flush=True)
            continue
        try:
            tc, camber, cx, te_gap = geom_measure(g['geo'])
        except Exception as e:
            results.append({'key': g['key'], 'status': 'FAIL',
                            'reason': 'exception: %s' % e})
            continue
        th_ok = abs(tc - g['thickness_tc']) <= 0.01
        cam_ok = abs(camber - g['camber_max']) <= 0.005
        te_ok = g.get('te_gap_chord') is None or \
            abs(te_gap - g['te_gap_chord']) <= 0.003
        inner = th_ok and cam_ok and te_ok
        rec = {'key': g['key'], 't_c': round(tc, 4),
               't_c_pub': g['thickness_tc'],
               'camber_max': round(camber, 4), 'camber_pub': g['camber_max'],
               'camber_x': round(cx, 3),
               'te_gap': round(te_gap, 4),
               'te_gap_pub': g.get('te_gap_chord'),
               'status': 'PASS' if inner else 'FAIL'}
        results.append(rec)
        print('  %-20s t/c=%.4f (pub %.3f)  camber=%.4f (pub %.3f)  '
              'TE=%.4f (pub %s) -> %s' % (
                  g['key'], tc, g['thickness_tc'], camber, g['camber_max'],
                  te_gap, g.get('te_gap_chord'), rec['status']),
              flush=True)


def main():
    args = [a for a in sys.argv[1:]]
    suite = 'all'
    targets = None
    label = 'last'
    ep = int(os.environ.get('WING3D_CTEST_EP', '8'))
    span = float(os.environ.get('WING3D_TSD_SPAN', '6.0'))
    grid = [int(x) for x in os.environ.get('WING3D_CTEST_GRID',
                                           '80 32 32').split()]
    i = 0
    while i < len(args):
        a = args[i]
        if a == '--suite' and i + 1 < len(args):
            suite = args[i + 1]; i += 2
        elif a == '--targets' and i + 1 < len(args):
            targets = args[i + 1].split(','); i += 2
        elif a == '--label' and i + 1 < len(args):
            label = args[i + 1]; i += 2
        elif a == '--ep' and i + 1 < len(args):
            ep = int(args[i + 1]); i += 2
        elif a == '--span' and i + 1 < len(args):
            span = float(args[i + 1]); i += 2
        elif a == '--grid' and i + 3 < len(args):
            grid = [int(x) for x in args[i + 1:i + 4]]; i += 4
        else:
            i += 1
    free_k = os.environ.get('WING3D_FREE_K', '') == '1'
    bc_rows = os.environ.get('WING3D_BC_ROWS', '') == '1'
    results = []
    print('== ctest[%s] label=%s grid=%dx%dx%d ep=%d span=%.1f '
          'free_k=%s bc_rows=%s ==' % (
              label, suite, grid[0], grid[1], grid[2], ep, span,
              free_k, bc_rows), flush=True)
    if suite in ('flow', 'all'):
        for t in REFS['flow']:
            if targets and t['key'] not in targets:
                continue
            run_target(t, grid, ep, span, label, results, free_k, bc_rows)
        grade_ratios(results)
    if suite in ('geometry', 'all'):
        run_geometry(results)
    resdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'results')
    os.makedirs(resdir, exist_ok=True)
    out = os.path.join(resdir, '%s.json' % label)
    json.dump({'label': label, 'results': results}, open(out, 'w'),
              indent=1)
    fails = [r for r in results if r['status'] == 'FAIL']
    pend = [r for r in results if r['status'] == 'PENDING']
    print('== ctest %s: %d/%d PASS, %d PENDING -> %s' % (
        label, len(results) - len(fails) - len(pend), len(results),
        len(pend), out), flush=True)
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()