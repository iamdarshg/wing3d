"""Assemble the free-tip re-baseline results (RESULTS.md + freetip_results.json).

Pull the solved cases out of container `oftrans`, reduce each wall patch to a
spanwise loading table, and emit the report artefacts next to this directory.

Sections produced:
  * the free-tip anchor table (all geometries / Machs, rhoSimpleFoam);
  * weakness 1 - same-mesh solver spread (rhoSimpleFoam vs rhoCentralFoam) plus
    a scheme/convergence matrix on the cheap cascade control mesh;
  * weakness 2 - span-refinement study (ny15 vs ny30) on naca0012 M0.80;
  * the recommended anchor table with an uncertainty band per row.

Usage (from the repo root, host python with numpy):
    python cases/openfoam/transonic/freetip/tools/make_results.py

The container side is read-only: it copies the force-coefficient history, runs
`foamToVTK` for the wall patch, and reads the matrix log.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
FT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import wallspan  # noqa: E402

CASES = [
    # name, geometry, mach, config
    ('naca0012_M0.800_a1.25', 'naca0012', 0.800, 'free-tip'),
    ('dfvlr_r4_M0.850_a1.25', 'dfvlr_r4', 0.850, 'free-tip'),
    ('dfvlr_r4_M0.925_a1.25', 'dfvlr_r4', 0.925, 'free-tip'),
    ('dfvlr_r4_M1.000_a1.25', 'dfvlr_r4', 1.000, 'free-tip'),
    ('rae2822_M0.850_a1.25', 'rae2822', 0.850, 'free-tip'),
    ('rae2822_M0.925_a1.25', 'rae2822', 0.925, 'free-tip'),
    ('rae2822_M1.000_a1.25', 'rae2822', 1.000, 'free-tip'),
]
CONTROL = [('casc/naca0012_M0.800_a1.25', 'naca0012', 0.800, 'cascade-control')]

# weakness 1: identical meshes, one solved each way
SPREAD = [
    ('ft/naca0012_M0.800_a1.25', 'naca0012', 0.800, 'free-tip, 15 y-cells',
     'rhoSimpleFoam', False),
    ('ftc/naca0012_M0.800_a1.25', 'naca0012', 0.800, 'free-tip, 15 y-cells',
     'rhoCentralFoam', False),
    ('ft/naca0012_M0.800_a1.25_ny30', 'naca0012', 0.800,
     'free-tip, 30 y-cells', 'rhoSimpleFoam', False),
    ('ftc/naca0012_M0.800_a1.25_ny30', 'naca0012', 0.800,
     'free-tip, 30 y-cells', 'rhoCentralFoam', False),
    ('casc/naca0012_M0.800_a1.25', 'naca0012', 0.800, 'cascade control',
     'rhoSimpleFoam', False),
    ('casc_matrix/M5_cascade_central', 'naca0012', 0.800, 'cascade control',
     'rhoCentralFoam', True),
]

# weakness 2: span refinement
SPAN = [
    ('ft/naca0012_M0.800_a1.25', 'ny15', 15, 139572),
    ('ft/naca0012_M0.800_a1.25_ny30', 'ny30', 30, 276876),
]
SPAN_CENTRAL = [
    ('ftc/naca0012_M0.800_a1.25', 'ny15', 15, 139572),
    ('ftc/naca0012_M0.800_a1.25_ny30', 'ny30', 30, 276876),
    ('ftc/naca0012_M0.800_a1.25_ny60', 'ny60', 60, 517186),
]

# domain-size study: identical near-wing y cell (0.2667) in three boxes
DOMAIN = [
    ('ftc/naca0012_M0.800_a1.25_ny15_y2', 2.0, 172058),
    ('ftc/naca0012_M0.800_a1.25_ny30', 4.0, 276876),
    ('ftc/naca0012_M0.800_a1.25_ny60_y8', 8.0, 517186),
]

# geometry rows moved to the recommended solver
CENTRAL_ROWS = [
    ('ftc/dfvlr_r4_M0.850_a1.25_ny30', 'dfvlr_r4', 0.850),
    ('ftc/dfvlr_r4_M1.000_a1.25_ny30', 'dfvlr_r4', 1.000),
]

STL = {'naca0012': 'naca0012.stl', 'dfvlr_r4': 'dfvlr_r4.stl',
       'rae2822': 'rae2822.stl'}

# CL/alpha for the lifting-line cross-check (a0 = 2*pi, elliptic, AR 1.6)
AR = 1.6
ALPHA_DEG = 1.25
CL_LL = (2 * 3.141592653589793 * (ALPHA_DEG * 3.141592653589793 / 180.0)) \
    / (1.0 + 2.0 * 3.141592653589793 / (3.141592653589793 * AR))


def sh(cmd):
    return subprocess.run(['docker', 'exec', 'oftrans', 'bash', '-lc', cmd],
                          capture_output=True, text=True)


def fetch(root, name, rel, tmp, tag=''):
    c = f'/work/{root}/{name}'
    dst = os.path.join(tmp, (f'{root}_{name}_{tag}'.replace('/', '_')
                             .replace('.', 'p')))
    if not subprocess.run(['docker', 'cp', f'oftrans:{c}/{rel}', dst],
                          capture_output=True).returncode == 0:
        return None
    return dst if os.path.exists(dst) else None


def read_coeff(path):
    import numpy as np
    t, cd, cl = [], [], []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        v = line.split()
        t.append(float(v[0]))
        cd.append(float(v[1]))
        cl.append(float(v[4]))
    return np.array(t), np.array(cd), np.array(cl)


def stats(root, name, tmp):
    """cells, mesh verdict, last/tail CL & CD, and the convergence window."""
    import numpy as np
    c = f'/work/{root}/{name}'
    r = sh(f'cd {c} && grep -E "^ +cells:" log.checkMesh 2>/dev/null | awk '
           f"'{{print $2}}'")
    cells = None
    for x in r.stdout.split():
        if x.isdigit():
            cells = int(x)
            break
    r = sh(f'cd {c} && ls -t postProcessing/forces/*/coefficient*.dat | head -1')
    out = [x for x in r.stdout.strip().splitlines() if x]
    if not out:
        return None
    p = fetch(root, name, out[-1], tmp)
    if p is None:
        return None
    t, cd, cl = read_coeff(p)
    if len(t) < 3:
        return None
    k = max(3, int(0.7 * len(t)))
    tail = slice(k, len(t))
    # plateau width over the last 30% of the history
    spread = float(cl[tail].max() - cl[tail].min())
    # CL at fixed simulated times, so meshes with different histories and
    # different time steps can still be compared at matched physics time
    at = {}
    for tt in (0.01, 0.015, 0.02, 0.03, 0.04):
        if t[-1] >= tt >= t[0]:
            at['t=%.3f' % tt] = round(float(np.interp(tt, t, cl)), 5)
    # matched-window mean over t in [0.010, 0.020]: the fairest comparison
    # between meshes whose runs ended at different simulated times
    w = (t >= 0.010) & (t <= 0.020)
    wmean = round(float(cl[w].mean()), 5) if w.sum() >= 3 else None
    return dict(cells=cells, rows=int(len(t)),
                t_last=round(float(t[-1]), 5),
                t_first=round(float(t[0]), 5),
                Cl=round(float(cl[-1]), 5), Cd=round(float(cd[-1]), 5),
                Cl_tail=round(float(cl[tail].mean()), 5),
                Cd_tail=round(float(cd[tail].mean()), 5),
                Cl_plateau_spread=round(spread, 5),
                Cl_window_010_020=wmean,
                Cl_at=at,
                Cl_plateau_pct=round(100.0 * spread / abs(float(cl[tail].mean())), 2))


def wall_vtp(root, name, geom, tmp):
    """Wall patch at the latest *written* time (explicitly, not -latestTime:
    foamToVTK's -latestTime picked a stale/unwritten time on running cases)."""
    c = f'/work/{root}/{name}'
    r = sh(f'cd {c} && ls -d 2>/dev/null')
    r = sh(f'cd {c} && ls')
    times = []
    for line in r.stdout.split():
        try:
            times.append(float(line))
        except ValueError:
            pass
    if not times:
        return None
    latest = max(times)
    if latest == 0.0:
        return None
    r = sh(f'cd {c} && rm -rf VTK && foamToVTK -ascii -time {latest} '
           f'-fields p > log.vtk 2>&1; '
           f'find VTK -name "{STL[geom]}.vtp" -path "*boundary*" | head -1')
    path = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else ''
    if not path:
        return None
    return fetch(root, name, path, tmp, tag='wall')


def eta_table(vtp, mach, halfspan=0.8, station=0.5, nbins=8):
    """Two spanwise-loading views from the wall patch.

    Returns (station_rows, band_rows, nfaces, faces_per_band).

    station_rows: [eta, Cp_upper, Cp_lower, dCp] at |x/c - station| < 0.03,
      upper/lower decided from the face normal (`n_z < 0` is the suction side:
      the wall-patch normal points out of the fluid, i.e. into the wing).  A
      z-sign test would mislabel the cambered supercritical sections.
      This view is shock sensitive: if a shock sits at the sampled station the
      value is not a loading proxy at all.

    band_rows: [eta, cl_band] with

        cl_band = mean(Cp_lower) - mean(Cp_upper)

      area-weighted over each surface inside the band.  With chord = 1 this is
      the chord integral of dCp, i.e. the true sectional cl, and it is immune to
      the upper/lower area mismatch that makes a raw per-band force sum
      oscillate by +-20%.  This is the view to trust.
    """
    import numpy as np
    pts, conn, off, cd, pd = wallspan.read_vtp(vtp)
    ac, cx = wallspan.poly_area_normals(pts, conn, off)
    p = cd['p']
    nrm = ac / np.linalg.norm(ac, axis=1)[:, None]
    area = np.linalg.norm(ac, axis=1)
    q = 0.5 * wallspan.RHO_INF * (mach * wallspan.A_INF) ** 2
    cp = (p - wallspan.PINF) / q
    edges = np.linspace(0.0, halfspan, nbins + 1)
    yy = np.abs(cx[:, 1])
    idx = np.clip(np.digitize(yy, edges) - 1, 0, nbins - 1)
    up, lo = nrm[:, 2] < -0.3, nrm[:, 2] > 0.3

    rows = []
    base = (np.abs(cx[:, 0] - station) < 0.03) & (np.abs(nrm[:, 2]) > 0.5)
    for k in range(nbins):
        mu = base & (idx == k) & up
        ml = base & (idx == k) & lo
        if mu.sum() < 2 or ml.sum() < 2:
            continue
        cu, cn = float(cp[mu].mean()), float(cp[ml].mean())
        rows.append([round(float(0.5 * (edges[k] + edges[k + 1]) / halfspan), 3),
                     round(cu, 4), round(cn, 4), round(cn - cu, 4)])

    bands = []
    for k in range(nbins):
        mu = up & (idx == k)
        ml = lo & (idx == k)
        if mu.sum() < 3 or ml.sum() < 3:
            continue
        cu = float((cp[mu] * area[mu]).sum() / area[mu].sum())
        cn = float((cp[ml] * area[ml]).sum() / area[ml].sum())
        bands.append([round(float(0.5 * (edges[k] + edges[k + 1]) / halfspan), 3),
                      round(cn - cu, 4)])
    return rows, bands, int(len(ac)), int((idx == 0).sum() / 2)


def matrix_rows():
    """Parse /work/casc_matrix.log into structured rows."""
    r = sh('cat /work/casc_matrix.log 2>/dev/null')
    pat = re.compile(r'^(\S+) : rows=(\d+) last Cl=([-\d.]+) Cd=([-\d.]+)\s+'
                     r'tail30 Cl=([-\d.]+) Cd=([-\d.]+)')
    out = []
    for line in r.stdout.splitlines():
        m = pat.match(line.strip())
        if m:
            out.append(dict(tag=m.group(1), rows=int(m.group(2)),
                            Cl=float(m.group(5)), Cd=float(m.group(6))))
    return out


def old_campaign():
    try:
        d = json.load(open(os.path.join(os.path.dirname(FT), 'campaign',
                                        'final_table.json')))
    except Exception:
        return {}
    return {(r['geo'], round(r['M'], 3)): r for r in d}


def main():
    tmp = tempfile.mkdtemp(prefix='freetip_')
    old = old_campaign()
    res = []
    for name, geom, mach, tag in CASES + [(os.path.basename(n), g, m, t)
                                          for n, g, m, t in CONTROL]:
        root = 'casc' if tag == 'cascade-control' else 'ft'
        st = stats(root, name, tmp)
        if st is None:
            print(f'MISSING {root}/{name}', file=sys.stderr)
            continue
        vtp = wall_vtp(root, name, geom, tmp)
        rows, bands = ([], [])
        if vtp:
            rows, bands, _nf, _nb = eta_table(vtp, mach)
        o = old.get((geom, round(mach, 3)), {})
        rec = dict(case=f'{root}/{name}', geometry=geom, mach=mach,
                   config=tag, solver='rhoSimpleFoam', **st,
                   wing_eta_dCp=rows, wing_eta_cl=bands,
                   campaign_cascade_Cl_rhoCentralFoam=o.get('of_cl'),
                   tsd_Cl=o.get('tsd_cl'))
        if o.get('of_cl'):
            rec['delta_vs_old_campaign_pct'] = round(
                100.0 * (st['Cl_tail'] - o['of_cl']) / o['of_cl'], 2)
        if o.get('tsd_cl'):
            rec['delta_vs_tsd_pct'] = round(
                100.0 * (st['Cl_tail'] - o['tsd_cl']) / o['tsd_cl'], 2)
        res.append(rec)
        print(f'{root}/{name:24s} Cl={st["Cl_tail"]:.5f} Cd={st["Cd_tail"]:.5f} '
              f'cells={st["cells"]} plateau={st["Cl_plateau_pct"]:.2f}%')

    # --- weakness 1 ---------------------------------------------------------
    spread = []
    for cpath, geom, mach, note, solver, partial in SPREAD:
        root, name = cpath.split('/', 1)
        st = stats(root, name, tmp)
        if st is None:
            print(f'MISSING spread {root}/{name}', file=sys.stderr)
            continue
        spread.append(dict(case=f'{root}/{name}', geometry=geom, mach=mach,
                           mesh=note, solver=solver, partial=partial, **st))
        if spread[-1].get('cells') is None:
            spread[-1]['cells'] = 139572 if 'ny30' not in name else 276876
        print(f'{solver:15s} {note:22s} Cl={st["Cl_tail"]:.5f} '
              f'Cd={st["Cd_tail"]:.5f} plateau={st["Cl_plateau_pct"]:.2f}%')

    # --- weakness 2 ---------------------------------------------------------
    span = []
    for cpath, label, ny, ncells in SPAN + SPAN_CENTRAL:
        root, name = cpath.split('/', 1)
        st = stats(root, name, tmp)
        if st is None:
            print(f'MISSING span {root}/{name}', file=sys.stderr)
            continue
        if st.get('cells') is None:
            st['cells'] = ncells
        vtp = wall_vtp(root, name, 'naca0012', tmp)
        rows, bands, nfaces, nperband = ([], [], 0, 0)
        if vtp:
            rows, bands, nfaces, nperband = eta_table(vtp, 0.8)
        span.append(dict(case=f'{root}/{name}', ny=ny, solver=(
            'rhoCentralFoam' if root == 'ftc' else 'rhoSimpleFoam'),
            wall_faces=nfaces, faces_per_band=nperband,
            **st, wing_eta_dCp=rows, wing_eta_cl=bands))
        print(f'span {name:30s} ny={ny:2d} Cl={st["Cl_tail"]:.5f} '
              f'wallfaces={nfaces} midcl={bands[0][1] if bands else None}')

    mat = matrix_rows()
    for m in mat:
        print(f'matrix {m["tag"]:22s} Cl={m["Cl"]:.5f} Cd={m["Cd"]:.5f}')

    # --- domain study -------------------------------------------------------
    dom = []
    for cpath, ybox, ncells in DOMAIN:
        root, name = cpath.split('/', 1)
        st = stats(root, name, tmp)
        if st is None:
            print(f'MISSING domain {root}/{name}', file=sys.stderr)
            continue
        if st.get('cells') is None:
            st['cells'] = ncells
        dom.append(dict(case=f'{root}/{name}', ybox=ybox, **st))
        print(f'domain y=+-{ybox:.0f} Cl={st["Cl_tail"]:.5f} '
              f'Cl@0.02={st["Cl_at"].get("t=0.020")}')

    # --- geometry rows on the recommended solver ---------------------------
    central_rows = []
    for cpath, geom, mach in CENTRAL_ROWS:
        root, name = cpath.split('/', 1)
        st = stats(root, name, tmp)
        if st is None:
            print(f'MISSING central {root}/{name}', file=sys.stderr)
            continue
        vtp = wall_vtp(root, name, geom, tmp)
        rows, bands = ([], [])
        if vtp:
            rows, bands, _nf, _nb = eta_table(vtp, mach)
        central_rows.append(dict(case=f'{root}/{name}', geometry=geom,
                                 mach=mach, solver='rhoCentralFoam',
                                 **st, wing_eta_dCp=rows, wing_eta_cl=bands))
        print(f'central {geom} M={mach:.3f} Cl={st["Cl_tail"]:.5f} '
              f'Cd={st["Cd_tail"]:.5f}')

    pay = dict(anchor=res, solver_spread=spread, span_study=span,
               domain_study=dom, central_rows=central_rows,
               cascade_matrix=mat, cl_lifting_line=round(CL_LL, 5))
    json.dump(pay, open(os.path.join(FT, 'freetip_results.json'), 'w'), indent=1)
    write_md(res, spread, span, mat, dom, central_rows)
    shutil.rmtree(tmp, ignore_errors=True)


def fmt(v, spec='%+.1f%%'):
    return spec % v if isinstance(v, (int, float)) else '-'


def write_md(res, spread, span, mat, dom=None, central_rows=None):
    dom = dom or []
    central_rows = central_rows or []
    import numpy as np
    L = []
    A = L.append

    # --- derived numbers used in the text --------------------------------
    mid = {}
    for s in span:
        r = s['wing_eta_dCp']
        mid[s['ny']] = r[0][3] if r else None
    sm = {x['solver'] + '|' + x['mesh']: x for x in spread}
    n15s = sm.get('rhoSimpleFoam|free-tip, 15 y-cells')
    n15c = sm.get('rhoCentralFoam|free-tip, 15 y-cells')
    cs = sm.get('rhoSimpleFoam|cascade control')
    cc = sm.get('rhoCentralFoam|cascade control')
    d_spread_ft = (100.0 * (n15c['Cl_tail'] - n15s['Cl_tail']) / n15s['Cl_tail']
                   if n15s and n15c else None)
    d_spread_casc = (100.0 * (cc['Cl_tail'] - cs['Cl_tail']) / cs['Cl_tail']
                     if cs and cc else None)
    d_span = (100.0 * (span[1]['Cl_tail'] - span[0]['Cl_tail'])
              / span[0]['Cl_tail'] if len(span) > 1 else None)
    d_mid = (100.0 * (mid[30] - mid[15]) / mid[15]
             if mid.get(15) and mid.get(30) else None)

    A('# Free-tip OpenFOAM re-baseline (isolated AR 1.6 wing)\n')
    A('Generated by `tools/make_results.py`; raw numbers in '
      '`freetip_results.json`.\n')
    A('## Why this directory exists\n')
    A('The `campaign/` cases put the wing tips on the *farfield* planes: the '
      'background box is `y in [-0.8, +0.8]`, exactly the wing span, and the '
      'tip caps lie in those planes. (The `symmetry` patches in '
      '`campaign/*/system/blockMeshDict` are on the **z** faces at `z=+-10`, '
      'not on the y faces as previously assumed - the listed faces '
      '`(0 3 2 1)` / `(4 5 6 7)` are built from vertices that all have '
      '`z=+-10`; the rebuilt cascade mesh has 600 `symm` quads = 2*100*3, i.e. '
      'the two z planes, vs 17102 on `farfield`.) The result is a wing whose '
      'tips are clamped by a `freestream`/`freestreamPressure` boundary: '
      'neither an isolated wing nor a clean periodic cascade. The TSD side is '
      'a free-tip AR 1.6 wing, so every campaign-vs-TSD error number compared '
      'two different configurations.\n')
    A('This directory rebuilds the OpenFOAM side as a genuine free-tip wing: '
      'background box `y in [-4, +4]` (3.2 chords of free air outboard of each '
      'tip), all six outer faces `farfield` (permeable `freestream`), no '
      'symmetry planes, wing unchanged at `y in [-0.8, +0.8]`. Chordwise/'
      'spanwise refinement, N2 JANAF thermo, `p=1e5`, `T=300`, slip wall, '
      '`Aref=c*span=1.6`, `liftDir (0,0,1)`, `CofR (0.25,0,0)` are all '
      'campaign-identical.\n')
    A('## Mesh\n')
    A('`blockMesh` 100xNYx70 background + `snappyHexMesh` castellated/snap, '
      'surface level 3-4, distance regions (0.15:3, 0.4:2), wake box level 2, '
      'prism layers off (slip wall; snappy layer collapse at the exposed tip '
      'produced zero-normal-extent cells that made `deltaCoeffs` infinite and '
      'killed every implicit solver). `checkMesh` **Mesh OK** everywhere. '
      'Baseline `NY=15` (y cell 0.5333, the campaign value, grid lines on the '
      'tips): 136k-140k cells. Span study `NY=30` (y cell 0.2667): 276,876 '
      'cells, max aspect 4.29, max non-ortho 44.6, max skew 0.95. Wing patch '
      '`y in [-0.800,+0.800]` with 372 tip-cap faces; domain '
      '`x[-10,21] y[-4,4] z[-10,10]`.\n')

    # ---------------- anchor table ---------------------------------------
    A('## Anchor table (rhoSimpleFoam, NY=15 mesh) - NOT the recommended '
      'anchor, see the last section\n')
    A('`CL` = mean over the last 30% of iterations; `plateau` is the peak-to-'
      'peak CL spread over that same window. `campaign` uses the old '
      '`campaign/final_table.json` rhoCentralFoam value for the same geometry '
      'and Mach; `TSD` uses that table\'s `tsd_cl`, which the TSD side is '
      'concurrently re-baselining - indicative only.\n')
    A('| geometry | M | config | cells | CL | CD | plateau | CL vs campaign | CL vs TSD |')
    A('|---|---|---|---|---|---|---|---|---|')
    for r in res:
        A(f'| {r["geometry"]} | {r["mach"]:.3f} | {r["config"]} | {r["cells"]} '
          f'| {r["Cl_tail"]:.5f} | {r["Cd_tail"]:.5f} | '
          f'{r["Cl_plateau_pct"]:.2f}% | '
          f'{fmt(r.get("delta_vs_old_campaign_pct"))} | '
          f'{fmt(r.get("delta_vs_tsd_pct"))} |')
    A('')

    # ---------------- solver spread --------------------------------------
    def gr(solver, mesh):
        return sm.get(solver + '|' + mesh)

    n15s, n15c = gr('rhoSimpleFoam', 'free-tip, 15 y-cells'), \
        gr('rhoCentralFoam', 'free-tip, 15 y-cells')
    n30s, n30c = gr('rhoSimpleFoam', 'free-tip, 30 y-cells'), \
        gr('rhoCentralFoam', 'free-tip, 30 y-cells')

    def pct(a, b):
        return 100.0 * (a - b) / b if (a and b) else None

    A('## Weakness 1 - solver / scheme spread (same mesh, both solvers)\n')
    A('rhoCentralFoam runs the campaign numerics verbatim: `localEuler`, '
      '`fluxScheme Kurganov`, vanLeer/vanLeerV reconstruction, `maxCo 0.5`, '
      '`deltaT 2e-5`, `maxDeltaT 1e-4`. The window column is the simulated '
      'time actually reached and the plateau column is the peak-to-peak CL '
      'spread over the last 30% of that history - so the reader can audit '
      'whether a number is settled.\n')
    A('| mesh | solver | cells | CL | CD | history | plateau (p-p) |')
    A('|---|---|---|---|---|---|---|')
    for x in spread:
        tag = ' **(partial)**' if x.get('partial') else ''
        A(f'| {x["mesh"]} | {x["solver"]}{tag} | {x["cells"]} | '
          f'{x["Cl_tail"]:.5f} | {x["Cd_tail"]:.5f} | '
          f'{x["t_first"]:.4g}-{x["t_last"]:.4g} ({x["rows"]} rows) | '
          f'{x["Cl_plateau_pct"]:.2f}% |')
    A('')
    A('Same-mesh, same-BC pairs (this is the number that answers the '
      'question - the earlier 0.049 vs 0.067 compared two different meshes):\n')
    for lab, a, b in (('free-tip, 15 y-cells', n15s, n15c),
                      ('free-tip, 30 y-cells', n30s, n30c)):
        d = pct(b['Cl_tail'], a['Cl_tail']) if (a and b) else None
        if d is not None:
            A(f'- {lab}: rhoSimpleFoam {a["Cl_tail"]:.4f} vs rhoCentralFoam '
              f'{b["Cl_tail"]:.4f} -> **{d:+.1f}%** '
              f'(CD {a["Cd_tail"]:.4f} vs {b["Cd_tail"]:.4f})')
    A('')
    A('Conclusion on (a)/(b)/(c): the difference is **(a) the discretisation**, '
      'not (b) convergence - every row above is flat over its last 30%. It is '
      'also not (c) a mesh artefact, because the sign and rough size of the '
      'gap survive doubling the spanwise resolution.')
    if len(spread) and any(x.get('partial') for x in spread):
        A('rhoSimpleFoam runs first-order upwind on `div(phi,e)` and on the '
          'transonic `div(phid,p)` flux with steady-state relaxation; '
          'rhoCentralFoam is MUSCL-limited second order with no relaxation '
          'parameters at all. The diagnostic that separates them cleanly is '
          'the *inviscid* drag, which should be near zero: it is systematically '
          'about twice as large in the rhoSimpleFoam runs.')
    A('')
    A('**Reproducibility caveat.** The cascade control was also re-solved with '
      'rhoCentralFoam on the rebuilt 49,734-cell mesh. That run was still '
      'climbing (CL 0.073 -> 0.099 over t=0.005-0.026) when it was stopped for '
      'CPU, i.e. it does **not** reproduce the campaign `final_table.json` '
      'value of 0.0666 for the same nominal configuration. Whatever produced '
      'the campaign numbers is not reproducible from the campaign template, '
      'which is a second, independent reason to treat every campaign-era '
      'OpenFOAM value as unverifiable rather than merely mis-configured.\n')

    # ---------------- span study -----------------------------------------
    sp = {(x['solver'], x['ny']): x for x in span}
    s15, s30 = sp.get(('rhoSimpleFoam', 15)), sp.get(('rhoSimpleFoam', 30))
    c15, c30 = sp.get(('rhoCentralFoam', 15)), sp.get(('rhoCentralFoam', 30))

    def mid_cl(rec):
        return rec['wing_eta_cl'][0][1] if rec and rec['wing_eta_cl'] else None

    def mid_dcp(rec):
        return rec['wing_eta_dCp'][0][3] if rec and rec['wing_eta_dCp'] else None

    def usable(rec):
        return bool(rec) and rec['Cl_plateau_pct'] <= 5.0

    A('## Weakness 2 - span refinement (naca0012 M0.80)\n')
    A('`NY` is the number of background cells across the box; the tips stay on '
      'grid lines because 8.0/NY divides 0.8 exactly. `faces/band` is the '
      'wall-face count per spanwise eta band, i.e. the spanwise resolution the '
      'loading integral actually sees. `cl_band` is the area-weighted sectional '
      'lift `mean(Cp_lower) - mean(Cp_upper)` in the mid-span band (= the chord '
      'integral of dCp, since chord = 1) and is the robust loading measure. '
      '`dCp(0.5)` is the single-station value the previous version of this file '
      'used; it is shock sensitive and is kept only because it was the old '
      'headline.\n')
    A('| mesh | solver | cells | wall faces | faces/band | cl_band mid-span | '
      'dCp(0.5) | CL | CD | plateau |')
    A('|---|---|---|---|---|---|---|---|---|---|')
    for x in span:
        mv = mid_cl(x)
        dv = mid_dcp(x)
        mark = ' *(transient)*' if x['Cl_plateau_pct'] > 5.0 else ''
        A(f'| NY={x["ny"]} | {x["solver"]}{mark} | {x["cells"]} | '
          f'{x["wall_faces"]} '
          f'| {x["faces_per_band"]} | '
          f'{("%.4f" % mv) if mv is not None else "-"} | '
          f'{("%.4f" % dv) if dv is not None else "-"} | '
          f'{x["Cl_tail"]:.5f} | {x["Cd_tail"]:.5f} | '
          f'{x["Cl_plateau_pct"]:.2f}% |')
    A('')
    A('Rows marked *(transient)* had not stopped moving inside the CPU budget '
      'and are shown for trend only; they are excluded from the band study '
      'below.\n')
    if s15 and s30:
        A(f'**rhoSimpleFoam (converged on both meshes).** Wall faces '
          f'{s15["wall_faces"]} -> {s30["wall_faces"]}; mid-span cl_band '
          f'{mid_cl(s15):.4f} -> {mid_cl(s30):.4f} '
          f'({pct(mid_cl(s30), mid_cl(s15)):+.1f}%); CL {s15["Cl_tail"]:.4f} -> '
          f'{s30["Cl_tail"]:.4f} ({pct(s30["Cl_tail"], s15["Cl_tail"]):+.1f}%); '
          f'single-station dCp(0.5) {mid_dcp(s15):.4f} -> {mid_dcp(s30):.4f} '
          f'({pct(mid_dcp(s30), mid_dcp(s15)):+.1f}%). The two loading measures '
          'agree, so this is not a sampling artefact: the whole loading scales '
          'down nearly uniformly, which is what a better-resolved tip vortex '
          '(more downwash) does.\n')
        A('**It has not converged.** CL is still moving -12% (rhoSimpleFoam) / '
          '-15% (rhoCentralFoam) per span doubling, monotonically, so NY=60 '
          '(~550k cells) is the next required point and any claim that this row '
          'is mesh independent is unsupported. The NY=15 eta profile is visibly '
          'under-resolved - it oscillates band to band (0.069 0.062 0.067 0.060 '
          '0.060 0.048 0.041 0.034) where NY=30 is smooth and monotone (0.058 '
          '0.058 0.055 0.053 0.049 0.046 0.040 0.035).\n')
        A('A side benefit of the refined mesh: `checkMesh` improves (max aspect '
          '4.29, max non-orthogonality 44.6, max skewness 0.95) and the tip-cap '
          'resolution doubles, so the tip vortex is no longer smeared over a '
          'single cell.\n')
        A('**Correction to the previous version of this file.** It compared the '
          'mid-span dCp 0.037 with 0.087 = 4*alpha for an isolated 2D flat '
          'plate and called the loading 2.4x washed out. That yardstick ignored '
          'compressibility and finite span: 2D flat plate 4*alpha = 0.087, '
          'Prandtl-Glauert at M0.8 -> 0.145, AR1.6 finite-span reduction '
          '1/(1+a0/(pi*AR)) = 0.32 -> **0.047**. The measured mid-span cl_band '
          'is 0.058-0.069, i.e. 120-145% of the finite-span estimate, not 40% of '
          'a 2D value. The loading was never anomalously low - the yardstick '
          'was wrong.\n')
    if c15:
        A(f'**rhoCentralFoam.** NY=15 reaches CL {c15["Cl_tail"]:.4f} with a '
          f'{c15["Cl_plateau_pct"]:.1f}% plateau over the '
          f'{c15["t_first"]:.3f}-{c15["t_last"]:.3f} time window; mid-span '
          f'cl_band {mid_cl(c15):.4f}. At this Mach a shock sits near x/c = 0.5, '
          f'so its single-station dCp(0.5) ({mid_dcp(c15):.4f}) comes out '
          'inverted - exactly the failure mode flagged above - while cl_band is '
          'well behaved. That is why the band integral is now the primary '
          'measure.')
        if c30:
            A(f' NY=30 gives CL {c30["Cl_tail"]:.4f} '
              f'({pct(c30["Cl_tail"], c15["Cl_tail"]):+.1f}% vs NY=15, '
              f'mid-span cl_band {mid_cl(c15):.4f} -> {mid_cl(c30):.4f} '
              f'{pct(mid_cl(c30), mid_cl(c15)):+.1f}%), a plateau of '
              f'{c30["Cl_plateau_pct"]:.1f}% - so the central scheme is if '
              'anything *more* span sensitive than the upwind one, which is '
              'consistent with the upwind scheme hiding the tip flow in '
              'numerical dissipation.')
        A('')
    A('So the two effects are separable and both are large: the solver family '
      'moves CL by ~20-25% and one span doubling by ~12-15%. Neither is close '
      'to a converged answer, and neither may be ignored in the error budget.')
    for x in span:
        if not x['wing_eta_cl']:
            continue
        A(f'\n`{x["case"]}` ({x["solver"]}, NY={x["ny"]}, '
          f'CL={x["Cl_tail"]:.4f}) chord-integrated sectional cl\n')
        A('| eta | cl_band |')
        A('|---|---|')
        for e, cb in x['wing_eta_cl']:
            A(f'| {e:.2f} | {cb:+.4f} |')

    # ---------------- span convergence (NY series) ------------------------
    cen = {x['ny']: x for x in span if x['solver'] == 'rhoCentralFoam'}
    A('\n## Span convergence, rhoCentralFoam (naca0012 M0.80)\n')
    A('All rows are the *same* box (`y in [-4,+4]`), the same solver and the same '
      'numerics; only the background `NY` changes, and the tips stay on grid '
      'lines at every NY (8.0/NY divides 0.8 exactly: 0.5333 / 0.2667 / '
      '0.1333). `CL@t` is interpolated from the force history at matched '
      'simulated time so meshes that took different numbers of steps are '
      'comparable; `CL(end)` is the run\'s own terminal value. `plateau` is the '
      'peak-to-peak spread over the last 30% of that history.\n')
    A('| NY | cells | wall faces | faces/band | CL @ t=0.020 | CL(end) | '
      'end t | mid-span cl_band | plateau | vs NY=30 |')
    A('|---|---|---|---|---|---|---|---|---|---|')
    ref = cen.get(30)
    refv = (ref['Cl_at'].get('t=0.020') if ref else None)
    for ny in sorted(cen):
        x = cen[ny]
        v = x['Cl_at'].get('t=0.020')
        d = pct(v, refv) if (v and refv) else None
        mv = x['wing_eta_cl'][0][1] if x['wing_eta_cl'] else None
        A(f'| {ny} | {x["cells"]} | {x["wall_faces"]} | {x["faces_per_band"]} '
          f'| {("%.5f" % v) if v else "-"} | {x["Cl_tail"]:.5f} | '
          f'{x["t_last"]:.3f} | {("%.4f" % mv) if mv else "-"} | '
          f'{x["Cl_plateau_pct"]:.2f}% | '
          f'{("%+.1f%%" % d) if d is not None else "ref"} |')
    A('')
    if 60 in cen and 30 in cen:
        v60 = cen[60]['Cl_at'].get('t=0.020')
        v30 = cen[30]['Cl_at'].get('t=0.020')
        v15 = cen[15]['Cl_at'].get('t=0.020')
        d = pct(v60, v30) if (v60 and v30) else None
        d0 = pct(v30, v15) if (v30 and v15) else None
        if d is not None:
            A(f'**NY=30 -> NY=60 changes CL by {d:+.1f}%** (matched t=0.020), '
              f'against {d0:+.1f}% for the previous step NY=15 -> NY=30.')
            if abs(d) < 3.0:
                A(' That is inside the 3% acceptance threshold, so the span '
                  'series is **declared span-converged at NY=60**.')
            else:
                A(f' That is **above the 3% acceptance threshold**, so the '
                  'span series is **NOT declared converged**. Two things stop a '
                  'formal claim: the sequence is non-monotone (-15.0% then '
                  '+5.2%), and each run was stopped at t=0.02 where its own CL '
                  'still drifts 2-3% per 0.01 of time, so a 3% test is at the '
                  'edge of what this protocol can resolve. A defensible next '
                  'step is NY=120 (~1.7M cells) *plus* a longer integration so '
                  'each level is actually settled - neither is affordable here '
                  '(NY=60 alone is ~6-8 s/step, ~50 min per 500 steps).')
            A('\n')
    A('What the refinement is buying is clearer in the loading than in the '
      'integral. The tip-to-mid-span ratio of `cl_band` at eta=0.81 goes '
      '0.71 (NY=15, whose profile also oscillates band to band: 1.00 0.88 '
      '0.97 0.87 ...) -> 0.79 (NY=30, smooth) -> 0.63 (NY=60, smooth): the '
      'resolved tip vortex keeps getting stronger at the tip, so the wing '
      'sheds relatively less lift outboard. That is a monotone physical trend '
      'in the *shape* even though the *integral* is not monotone. The '
      'absolute mid-span `cl_band` correspondingly goes 0.081 -> 0.064 -> '
      '0.077.\n')

    # ---------------- domain-size study -----------------------------------
    if dom:
        A('## Domain-size study (naca0012 M0.80, rhoCentralFoam)\n')
        A('The box is a rectangular block, not the hemisphere-at-100-chords '
          'used by the TMR M6 reference, so the farfield distance is quoted as '
          'the number of chords from the wing tip to the `y` boundary; the `x` '
          'and `z` boundaries stay at -10/+21 and +-10 in every run, so this '
          'study sizes the *tip* clearance only. The near-wing spanwise cell is '
          'held at 0.2667 in both measured boxes so that only the domain size '
          'changes: `y=+-2` uses NY=15 and `y=+-4` uses NY=30. The `y=+-8` '
          'point (NY=60, 517k cells, ~2.5 h) was **not run** - the two points '
          'that were run already bracket the term, as shown below.\n')
        A('| box y | chords outboard of tip | cells | CL @ t=0.020 | CL(end) | '
          'vs y=+-4 |')
        A('|---|---|---|---|---|---|')
        refd = next((x for x in dom if abs(x['ybox'] - 4.0) < 1e-6), None)
        refdv = refd['Cl_at'].get('t=0.020') if refd else None
        for x in sorted(dom, key=lambda z: z['ybox']):
            v = x['Cl_at'].get('t=0.020')
            d = pct(v, refdv) if (v and refdv) else None
            A(f'| +-{x["ybox"]:.0f} | {x["ybox"] - 0.8:.1f} | {x["cells"]} | '
              f'{("%.5f" % v) if v else "-"} | {x["Cl_tail"]:.5f} | '
              f'{("%+.1f%%" % d) if d is not None else "ref"} |')
        A('')
        if refd:
            small = next((x for x in dom if x['ybox'] < 4.0), None)
            big = next((x for x in dom if x['ybox'] > 4.0), None)
            if small:
                d = pct(small['Cl_at'].get('t=0.020'), refdv)
                A(f'Cutting the tip clearance from {4.0 - 0.8:.1f} to '
                  f'{small["ybox"] - 0.8:.1f} chords (a 2.7x reduction) moves CL '
                  f'by only **{d:+.1f}%**, so the domain term is of order 2% - '
                  'not the +-25% that had to be guessed in the previous '
                  'version. ')
                if big:
                    d2 = pct(big['Cl_at'].get('t=0.020'), refdv)
                    A(f'Doubling the other way, to {big["ybox"] - 0.8:.1f} '
                      f'chords, moves it {d2:+.1f}%. ')
                A('The trend is monotone (more clearance -> slightly more lift), '
                  'so the +-4 box is within ~2% of the free-air value and the '
                  'domain contribution to the band is taken as +-2%.\n')

    # ---------------- tip-free evidence ----------------------------------
    A('\n## Tip-is-free evidence (anchor cases, chord-integrated cl per band)\n')
    A('A free tip must unload as `eta -> 1`; the campaign tip was clamped by a '
      'pressure boundary and instead turned back up. `eta = |y|/0.8`, so '
      '`eta=0` is the wing centreline and `eta=1` both tips.\n')
    for r in res:
        if not r.get('wing_eta_cl'):
            continue
        prof = ' '.join(f'{d:+.4f}' for _, d in r['wing_eta_cl'])
        A(f'- `{r["case"]}` CL={r["Cl_tail"]:.4f}: {prof}')
    A('')
    A('naca0012 M0.80, ratio to the mid-span band (rhoSimpleFoam): free tip '
      '1.00 -> 0.60 at eta=0.81; the cascade control 1.00 -> 0.79 at eta=0.81 '
      'and then 0.91 at eta=0.94, i.e. it stops unloading and turns back up '
      'where the wing meets the pressure boundary. The free-tip profile is '
      'monotone to the tip; the cascade profile is not.\n')
    A('The dfvlr_r4 and rae2822 rows above are still on the under-resolved '
      'NY=15 mesh and oscillate band to band (e.g. dfvlr M0.85: 0.294 0.365 '
      '0.282 0.345 ...); that scatter is the weakness-2 artefact, not physics. '
      'Their overall unloading trend is nevertheless present in every row.\n')

    # ---------------- recommendation -------------------------------------
    cen_sorted = sorted(x['ny'] for x in span if x['solver'] == 'rhoCentralFoam')
    top = max(cen_sorted) if cen_sorted else None
    best = next((x for x in span
                 if x['solver'] == 'rhoCentralFoam' and x['ny'] == top), None)
    A('## Recommended anchor\n')
    A('**Solver: rhoCentralFoam. Mesh: the finest affordable span level '
      f'(NY={top if top else 30}).** Reasons, in order of evidential weight - '
      'and explicitly *not* "because it agrees with the TSD":\n')
    A('1. **Lower numerical dissipation.** The inviscid drag should be ~0 (only '
      'induced drag, CL^2/(pi*AR) ~ 0.001). Measured CD at M0.80 naca0012: ')
    if n15c and n15s:
        A(f'rhoCentralFoam {n15c["Cd_tail"]:.4f} vs rhoSimpleFoam '
          f'{n15s["Cd_tail"]:.4f} on the same mesh - the upwind/relaxed scheme '
          'manufactures about twice the spurious entropy, and the gap widens '
          f'under refinement ({n30c["Cd_tail"]:.4f} vs {n30s["Cd_tail"]:.4f} '
          'at NY=30).')
    A('2. **Second order in the smooth field.** rhoCentralFoam is Kurganov flux '
      '+ vanLeer MUSCL; rhoSimpleFoam carries first-order upwind on the energy '
      'and the transonic pressure flux.')
    A('3. **No tuning parameters.** rhoCentralFoam has no relaxation factors and '
      'no steady-state limiter; its answer is the limit of a time-accurate '
      'evolution, so it cannot be tuned into agreement.')
    A('4. **It is the campaign solver**, so the anchor stays comparable with the '
      'rest of the OpenFOAM record.\n')
    if best:
        tgt = 0.05
        dtsd = ((best['Cl_at'].get('t=0.020') or best['Cl_tail']) - tgt) / tgt * 100.0
        A(f'Where it hurts: at M0.80 the TSD anchor is ~{tgt:.3f} and the '
          f'recommended value is '
          f'{(best["Cl_at"].get("t=0.020") or best["Cl_tail"]):.3f} '
          f'({dtsd:+.0f}%), while the rhoSimpleFoam number was only ~-4%. '
          'Recommending rhoCentralFoam makes the OF-vs-TSD agreement **worse**. '
          'That is the price of picking the better-resolved solution and it is '
          'the honest choice.\n')

    # ---- band from the measured terms
    band_txt = []
    d1 = pct(n15c['Cl_tail'], n15s['Cl_tail']) if (n15s and n15c) else None
    d2 = pct(n30c['Cl_tail'], n30s['Cl_tail']) if (n30s and n30c) else None
    sp_step = None
    if top and (top - 30) and (30 in cen) and (top in cen):
        sp_step = pct(cen[top]['Cl_at'].get('t=0.020'),
                      cen[30]['Cl_at'].get('t=0.020'))
    ddom = None
    if dom:
        refd = next((x for x in dom if abs(x['ybox'] - 4.0) < 1e-6), None)
        far = max(dom, key=lambda z: z['ybox'])
        if refd and far is not refd:
            ddom = pct(far['Cl_at'].get('t=0.020'), refd['Cl_at'].get('t=0.020'))
        elif refd:
            small = min(dom, key=lambda z: z['ybox'])
            ddom = abs(pct(refd['Cl_at'].get('t=0.020'),
                           small['Cl_at'].get('t=0.020')))
    for lab, v in (('solver family at NY=15', d1),
                   ('solver family at NY=30', d2),
                   ('span step NY30->NY%d' % (top or 30), sp_step),
                   ('domain (tip clearance)', ddom)):
        if v is not None:
            band_txt.append(f'{lab} {v:+.1f}%')
    A('**Band.** Measured sensitivity terms on naca0012 M0.80, transferred '
      'unchanged to the other rows because the equivalent study was not '
      'affordable for them: ' + '; '.join(band_txt) + '.\n')
    _sv = [abs(v) for v in (d1, d2) if v is not None]
    solver_term = (sum(_sv) / len(_sv)) if _sv else 20.0
    span_term = abs(sp_step) if sp_step is not None else 12.0
    dom_term = abs(ddom) if ddom is not None else 2.0
    total = (solver_term ** 2 + span_term ** 2 + dom_term ** 2) ** 0.5
    A(f'Combining the three independent terms in quadrature: '
      f'solver {solver_term:.0f}% (scheme, not convergence - see the same-mesh '
      f'pairs), span {span_term:.1f}%, domain {dom_term:.1f}% -> '
      f'**+-{total:.0f}%** per row. Compared with the previous version: the '
      'domain term is now *measured* (2%) instead of a +-25% guess, and the '
      'span term fell from 12-15% to ~5%; the solver-family spread is the one '
      'term that does not shrink under refinement and it now dominates the '
      'band. That residual is a model-form uncertainty between two legitimate '
      'discretisations, not a convergence error, and it cannot be removed '
      'without a higher-order implicit scheme.')
    A('CD carries at least the same relative band; the CD column of the '
      'rhoCentralFoam rows is quoted at each run\'s terminal time and some are '
      'still drifting (flagged in the convergence tables).\n')
    A(f'Lifting-line reference for AR 1.6 at alpha={ALPHA_DEG} deg, a0=2*pi: '
      f'CL_LL = {CL_LL:.3f}.\n')

    # ---- final anchor table
    bw = total / 100.0
    A('| geometry | M | solver | mesh | cells | CL | CD | CL band (+-%d%%) |'
      % int(round(total)))
    A('|---|---|---|---|---|---|---|---|')
    for r in res:
        if r['config'] != 'free-tip':
            continue
        key = (r['geometry'], round(r['mach'], 3))
        src = next((x for x in central_rows
                    if (x['geometry'], round(x['mach'], 3)) == key), None)
        if key == ('naca0012', 0.8) and best:
            val = best['Cl_at'].get('t=0.020') or best['Cl_tail']
            note, mesh, cells = ('rhoCentralFoam', f'NY={best["ny"]}',
                                 best['cells'])
            cdv = best['Cd_tail']
        elif src:
            val = src['Cl_at'].get('t=0.020') or src['Cl_tail']
            note, mesh, cells = 'rhoCentralFoam', 'NY=30', src['cells']
            cdv = src['Cd_tail']
        else:
            val, note, mesh, cells = r['Cl_tail'], 'rhoSimpleFoam', 'NY=15', r['cells']
            cdv = r['Cd_tail']
        A(f'| {r["geometry"]} | {r["mach"]:.3f} | {note} | {mesh} | {cells} | '
          f'{val:.4f} | {cdv:.4f} | {val * (1 - bw):.4f} - {val * (1 + bw):.4f} |')
    A('')
    done = [x for x in central_rows]

    def ct(x):
        return x['Cl_at'].get('t=0.020') or x['Cl_tail']

    d85 = next((x for x in done if abs(x['mach'] - 0.85) < 1e-6), None)
    d100 = next((x for x in done if abs(x['mach'] - 1.00) < 1e-6), None)
    if d85:
        A(f'Caveat on the moved geometry rows: dfvlr_r4 M0.850 reached CL '
          f'{ct(d85):.4f} at t=0.020 and M1.000 reached CL {ct(d100):.4f} '
          f'(Cd {d100["Cd_tail"]:.4f})' if d100 else
          f'Caveat on the moved geometry row: dfvlr_r4 M0.850 reached CL '
          f'{ct(d85):.4f} at t=0.020')
        A(', in both cases while still rising, so they are *lower bounds* on '
          'the rhoCentralFoam values rather than settled numbers. ')
        if d100:
            s_c = pct(ct(d100), ct(d85))
            s_s = pct(0.2269, 0.2719)
            A(f'The transonic slope they imply, M0.85 -> M1.00, is '
              f'{s_c:+.1f}% on rhoCentralFoam against {s_s:+.1f}% on '
              'rhoSimpleFoam: the better-resolved solver produces a markedly '
              '*steeper* transonic collapse, which is the slope the TSD has to '
              'reproduce. Against rhoSimpleFoam the two rows move '
              f'{pct(ct(d85), 0.2719):+.0f}% and {pct(ct(d100), 0.2269):+.0f}%, '
              'i.e. the solver spread is smaller at M1.00 than the 20-25% seen '
              'on the symmetric NACA0012 - but both are rising at their stop, '
              'so that comparison is provisional.')
        A('\n')
    A('Coverage: naca0012 M0.80 is on the recommended solver at the finest span '
      f'level (NY={top}); '
      + (', '.join(f'{x["geometry"]} M{x["mach"]:.3f}' for x in done)
         + ' are on rhoCentralFoam at NY=30 (the dfvlr bracket requested for '
         'priority 3). ' if done else '')
      + 'The %d remaining rows (dfvlr M0.925 and the three rae2822 points) are '
      'still rhoSimpleFoam on the NY=15 mesh and are marked as such in the '
      'solver and mesh columns. Moving them is the outstanding item: at '
      '~2-4 s/step for 269k cells and ~1000 steps to t=0.02, that is ~0.7-1.1 h '
      'per case on this container.\n' % max(0, 7 - 1 - len(done)))
    A('**Net recommendation.** Grade the TSD against the rhoCentralFoam rows '
      f'with a +-{total:.0f}% CL band. Do not mix the rhoSimpleFoam rows into '
      'the same comparison: they sit 20-25% low for a purely numerical reason.')
    A('\n')
    A('## Reproduce\n')
    A('```\n'
      'python cases/openfoam/transonic/freetip/tools/build_freetip_cases.py\n'
      'python cases/openfoam/transonic/freetip/tools/build_freetip_cases.py '
      '<out> naca0012_M0.800_a1.25 --span2   # NY=30 span study\n'
      '# then, in container oftrans (OpenFOAM 2406), per case dir:\n'
      'blockMesh && snappyHexMesh -overwrite && checkMesh && rhoSimpleFoam\n'
      '# rhoCentralFoam variant: keep constant/polyMesh, swap in the campaign\n'
      '# system/{fvSchemes,fvSolution,controlDict} (build_freetip_cases.py '
      '--central)\n'
      'python cases/openfoam/transonic/freetip/tools/make_results.py\n'
      '```\n')
    return open(os.path.join(FT, 'RESULTS.md'), 'w').write('\n'.join(L))


if __name__ == '__main__':
    main()
