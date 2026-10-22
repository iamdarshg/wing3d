"""Direct dfvlr M-probe helper for the A6 ladder (ctest span 6.0).

published_refs flow targets only include dfvlr_r4 at M0.925. The
campaign-matched 3-probe Ms are 0.85 / 0.925 / 1.00 (dfvlr), so the
0.85 and 1.00 dfvlr probes are run directly through this helper with
the same solve path as tsd_ablate.py (run_sup, g80x32x32 ep8, span
from WING3D_TSD_SPAN default 6.0, alpha 1.25) and the same knob env
(WING3D_FREE_K / WING3D_BC_ROWS / WING3D_TE_GAP). RAM-gated exactly
like the ladder (WING3D_ABLATE_RAM_GB default 2.5).

Usage:
  python tests3d/tsd_ablate_mp.py <tag> <M> [--span 6.0]
Results: tests3d/validation/results/ablate_mp_<tag>_M<M>.json
"""
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'tests3d'))


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
    tag = sys.argv[1]
    M = float(sys.argv[2])
    span = float(os.environ.get('WING3D_TSD_SPAN', '1.6'))
    if len(sys.argv) > 3 and sys.argv[3] == '--span':
        span = float(sys.argv[4])
    free_k = os.environ.get('WING3D_FREE_K', '') == '1'
    bc_rows = os.environ.get('WING3D_BC_ROWS', '') == '1'
    te_gap = float(os.environ.get('WING3D_TE_GAP', '0.002'))
    ep = int(os.environ.get('WING3D_CTEST_EP', '8'))
    print('== mp[%s] dfvlr M=%.3f span=%.1f ep=%d free_k=%s bc_rows=%s '
          'te_gap=%.4f ==' % (tag, M, span, ep, free_k, bc_rows, te_gap),
          flush=True)
    if not wait_ram():
        print('RAM-GATED skip M=%.3f' % M, flush=True)
        return 2
    from tsd_campaign import build_profile, build_section
    from tsd_m_sup import run_sup
    prof = build_profile('dfvlr_r4_true')
    sec = build_section('dfvlr_r4_true')
    t0 = time.time()
    L, xc, ix, j0 = run_sup(M, span=span, epochs=ep, free_k=free_k,
                            profile=prof, section=sec, alpha_deg=1.25)
    cl = float(L['CL'])
    secs = round(time.time() - t0, 1)
    outdir = os.path.join(ROOT, 'tests3d', 'validation', 'results')
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, 'ablate_mp_%s_M%.3f.json' % (tag, M))
    json.dump({'tag': tag, 'geo': 'dfvlr_r4', 'M': M, 'span': span, 'ep': ep,
               'free_k': free_k, 'bc_rows': bc_rows, 'te_gap': te_gap,
               'cl_tsd': round(cl, 6), 'secs': secs},
              open(out, 'w'), indent=1)
    print('  dfvlr M=%.3f CL=%.6f (%.0fs) -> %s' % (M, cl, secs, out),
          flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
