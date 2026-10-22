"""TSD ablation ladder -- exactly the ctest flow protocol, per-knob.

Protocol mirror of tests3d/validation/ctest_validate.py --suite flow:
  grid 80x32x32, ep 8 (WING3D_CTEST_EP), span from WING3D_TSD_SPAN
  (ctest default 6.0; pass --span 1.6 for the campaign-matched pairing),
  alpha per target from published_refs.json. Knobs read from env exactly
  as ctest reads them:
    WING3D_FREE_K  '1' -> free-k closure (default fixed k=1)
    WING3D_BC_ROWS '1' -> real farfield rows (default rank-deficient zeros)
    WING3D_TE_GAP  geometry knob (default 0.002)

Honest grading: a target is graded PASS/FAIL only when the referee
(published_refs) anchors it (cl_pub set AND status empty). PENDING
targets are still MEASURED at the probe (M, alpha) so the ablation cost
is one subtraction away when A5 lands the 3D RANS anchors -- but their
status stays 'MEASURED', never PASS/FAIL. onera_m6 has no 2D-section
TSD proxy (swept 3D wing; OF-side geometry) -> SKIP, never graded.

dfvlr span=1.6 runs also get the PROVISIONAL fine-branch (dfvlr_fine
125k) error, clearly flagged -- that branch is NOT the registered 3D RANS
anchor and is used for direction only.

RAM: like ctest, waits for >= WING3D_ABLATE_RAM_GB (default 2.5) free
before each solve (the OF container needs ~1.8-2.5 GB for its fine-band
solve; one TSD solve at a time, cooperative).

Usage:
  python tests3d/tsd_ablate.py <tag> [--targets k1,k2] [--span 6.0]
    env: WING3D_FREE_K / WING3D_BC_ROWS / WING3D_TE_GAP / WING3D_ALPHA /
         WING3D_CTEST_EP / WING3D_ABLATE_RAM_GB
  Results: tests3d/validation/results/ablate_<tag>.json
"""
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'tests3d'))

REFS = json.load(open(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'validation',
    'published_refs.json')))

# dfvlr fine-branch 125k values (in-flight dfvlr_fine program; PROVISIONAL
# only -- A5 registers the official 3D RANS anchors in published_refs).
DFVLR_FINE = {0.850: 0.3210697, 0.875: 0.3436987, 0.900: 0.3635240,
              0.925: 0.3629209, 0.950: 0.3464955}


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


def wait_ram(min_gb=None, timeout_s=7200):
    if min_gb is None:
        min_gb = float(os.environ.get('WING3D_ABLATE_RAM_GB', '2.5'))
    t0 = time.time()
    while time.time() - t0 < timeout_s:
        if ram_free_gb() >= min_gb:
            return True
        time.sleep(60)
    return False


def main():
    args = [a for a in sys.argv[1:]]
    tag = args[0] if args else 'ladder'
    targets = None
    span = float(os.environ.get('WING3D_TSD_SPAN', '6.0'))
    i = 1
    while i < len(args):
        if args[i] == '--targets' and i + 1 < len(args):
            targets = args[i + 1].split(','); i += 2
        elif args[i] == '--span' and i + 1 < len(args):
            span = float(args[i + 1]); i += 2
        else:
            i += 1
    free_k = os.environ.get('WING3D_FREE_K', '') == '1'
    bc_rows = os.environ.get('WING3D_BC_ROWS', '') == '1'
    te_gap = float(os.environ.get('WING3D_TE_GAP', '0.002'))
    ep = int(os.environ.get('WING3D_CTEST_EP', '8'))
    print('== ablate[%s] targets=%s span=%.1f ep=%d free_k=%s bc_rows=%s '
          'te_gap=%.4f ==' % (tag, targets, span, ep, free_k, bc_rows,
                              te_gap), flush=True)
    from tsd_campaign import build_profile, build_section
    from tsd_m_sup import run_sup
    results = []
    for t in REFS['flow']:
        if targets and t['key'] not in targets:
            continue
        geo, M, alpha = t['geometry'], float(t['M']), float(t['alpha'])
        anchored = (t.get('cl_pub') is not None
                    and t.get('status', '') in ('', None))
        if geo == 'onera_m6_sharp_te':
            results.append({'key': t['key'], 'status': 'SKIP',
                            'reason': 'no 2D TSD section proxy '
                                      '(swept 3D wing; diagnostic)'})
            print('  %-24s -> SKIP (no 2D section proxy)' % t['key'],
                  flush=True)
            continue
        rec = {'key': t['key'], 'geometry': geo, 'M': M, 'alpha': alpha,
               'cl_pub': t['cl_pub'] if anchored else None,
               'cl_band': t['cl_band'] if anchored else None,
               'status': None}
        if not anchored:
            rec['status'] = 'MEASURED'
            rec['note'] = 'PENDING until A5 lands 3D RANS anchor; ' \
                          'CL recorded for later grading'
        if not wait_ram():
            rec.update({'status': 'FAIL', 'reason': 'ram_gated',
                        'cl_tsd': None})
            results.append(rec)
            print('RAM-GATED skip %s' % t['key'], flush=True)
            continue
        try:
            prof = build_profile(geo)
            sec = build_section(geo)
            t0 = time.time()
            L, xc, ix, j0 = run_sup(M, span=span, epochs=ep, free_k=free_k,
                                    profile=prof, section=sec,
                                    alpha_deg=alpha)
            cl = float(L['CL'])
            secs = round(time.time() - t0, 1)
            rec.update({'cl_tsd': round(cl, 6), 'secs': secs})
            if anchored:
                err = (cl - t['cl_pub']) / t['cl_pub'] * 100
                inside = abs(cl - t['cl_pub']) <= t['cl_band']
                rec.update({'status': 'PASS' if inside else 'FAIL',
                            'err_pct': round(err, 2), 'inside': inside})
            # provisional fine-branch reference (dfvlr only), flagged
            if geo.startswith('dfvlr') and abs(span - 1.6) < 1e-9:
                f = DFVLR_FINE.get(M)
                if f:
                    rec['fine_prov'] = f
                    rec['err_fine_pct'] = round((cl - f) / f * 100, 2)
            print('  %-24s CL=%.5f %s%s' % (
                rec['key'], cl, rec['status'],
                ('  fine_prov=%+.1f%%' % rec['err_fine_pct'])
                if rec.get('err_fine_pct') is not None else ''), flush=True)
        except Exception as e:
            rec.update({'status': 'FAIL', 'reason': 'exception: %s' % e,
                        'cl_tsd': None})
        results.append(rec)
    outdir = os.path.join(ROOT, 'tests3d', 'validation', 'results')
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, 'ablate_%s.json' % tag)
    json.dump({'tag': tag, 'results': results,
               'env': {'free_k': free_k, 'bc_rows': bc_rows,
                       'te_gap': te_gap, 'span': span, 'ep': ep}},
              open(out, 'w'), indent=1)
    fails = [r for r in results if r['status'] == 'FAIL']
    print('== ablate %s: %d measured, %d FAIL, %d PENDING/SKIP -> %s' % (
        tag, sum(1 for r in results
                 if r['status'] in ('MEASURED', 'PASS', 'FAIL') and
                 r.get('cl_tsd') is not None),
        len(fails), sum(1 for r in results if r['status'] in ('PENDING',
                                                               'SKIP')),
        out), flush=True)
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
