"""U3 ablation: savgol section smoothing vs raw slopes (tsdcli STL path).

Audit 2026-10-22: savgol (w=41) section smoothing in tsdcli is
SHAKY/tuned ("parameters chosen until anchor CLs matched, rae +0.04%").
This test runs the shipped STL CLI on rae2822.stl at M0.85 with the
default (savgol) and WING3D_RAW_SECTION=1 (raw slopes) and records CL,
section cl, t/c, camber for both. Same RAM discipline as tsd_ablate
(gate WING3D_ABLATE_RAM_GB default 2.5, one solve at a time).

Usage: python tests3d/tsd_ablate_u3.py [--mach 0.85] [--epochs 8]
Results: tests3d/validation/results/ablate_u3.json
"""
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STL = os.path.join(ROOT, 'cases', 'openfoam', 'transonic', 'stl',
                   'rae2822.stl')


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
    mach = 0.85
    ep = int(os.environ.get('WING3D_CTEST_EP', '8'))
    i = 0
    while i < len(args):
        if args[i] == '--mach' and i + 1 < len(args):
            mach = float(args[i + 1]); i += 2
        elif args[i] == '--epochs' and i + 1 < len(args):
            ep = int(args[i + 1]); i += 2
        else:
            i += 1
    print('== U3 savgol ablation rae2822 M=%.3f ep=%d ==' % (mach, ep),
          flush=True)
    results = []
    for raw in (False, True):
        if not wait_ram():
            results.append({'raw_section': raw, 'status': 'FAIL',
                            'reason': 'ram_gated'})
            print('RAM-GATED skip raw_section=%s' % raw, flush=True)
            continue
        env = dict(os.environ)
        if raw:
            env['WING3D_RAW_SECTION'] = '1'
        else:
            env.pop('WING3D_RAW_SECTION', None)
        out = os.path.join(ROOT, 'tests3d', 'validation', 'results',
                           'ablate_u3_%s.json' % ('raw' if raw else 'savgol'))
        cmd = [sys.executable, '-m', 'wing3d.tsdcli', STL,
               '--mach', '%.3f' % mach, '--epochs', str(ep),
               '--json', out]
        t0 = time.time()
        r = subprocess.run(cmd, cwd=ROOT, env=env,
                           capture_output=True, text=True, timeout=7200)
        rec = {'raw_section': raw, 'rc': r.returncode,
               'secs': round(time.time() - t0, 1),
               'tail': r.stdout.strip().splitlines()[-4:] if r.stdout
               else r.stderr.strip().splitlines()[-4:]}
        if os.path.isfile(out):
            d = json.load(open(out))
            rec.update(cl=d['CL'], cd=d['CD'], section_cl=d['section_cl'],
                       thick_max=d['thick_max'], camber_max=d['camber_max'])
        results.append(rec)
        print(json.dumps(rec, indent=1), flush=True)
    res = {'mach': mach, 'epochs': ep, 'stl': STL, 'results': results}
    resdir = os.path.join(ROOT, 'tests3d', 'validation', 'results')
    os.makedirs(resdir, exist_ok=True)
    json.dump(res, open(os.path.join(resdir, 'ablate_u3.json'), 'w'),
              indent=1)
    print('saved ablate_u3.json', flush=True)


if __name__ == '__main__':
    main()
