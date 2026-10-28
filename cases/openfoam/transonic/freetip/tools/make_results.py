"""Assemble the free-tip re-baseline results (RESULTS.md + freetip_results.json).

Pull the solved cases out of container `oftrans`, reduce each wall patch to a
spanwise loading table, and emit the report artefacts next to this directory.

Usage (from the repo root, host python with numpy):
    python cases/openfoam/transonic/freetip/tools/make_results.py

The container side is read-only here: it re-runs only `foamToVTK` on the wall
patch and reads `postProcessing/forces/*/coefficient*.dat`.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
FT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import wallspan  # noqa: E402

CASES = [
    # name, geometry, mach, endpoint note
    ('naca0012_M0.800_a1.25', 'naca0012', 0.800, 'free-tip'),
    ('dfvlr_r4_M0.850_a1.25', 'dfvlr_r4', 0.850, 'free-tip'),
    ('dfvlr_r4_M0.925_a1.25', 'dfvlr_r4', 0.925, 'free-tip'),
    ('dfvlr_r4_M1.000_a1.25', 'dfvlr_r4', 1.000, 'free-tip'),
    ('rae2822_M0.850_a1.25', 'rae2822', 0.850, 'free-tip'),
    ('rae2822_M0.925_a1.25', 'rae2822', 0.925, 'free-tip'),
    ('rae2822_M1.000_a1.25', 'rae2822', 1.000, 'free-tip'),
]
CONTROL = [('casc/naca0012_M0.800_a1.25', 'naca0012', 0.800, 'cascade-control')]

STL = {'naca0012': 'naca0012.stl', 'dfvlr_r4': 'dfvlr_r4.stl',
       'rae2822': 'rae2822.stl'}


def sh(cmd):
    return subprocess.run(['docker', 'exec', 'oftrans', 'bash', '-lc', cmd],
                          capture_output=True, text=True)


def case_stats(root, name):
    """cells / last CL / Cd / converged-tail mean CL."""
    c = f'/work/{root}/{name}'
    r = sh(f'cd {c} && '
           f'grep -E "^ +cells:" log.checkMesh | awk \'{{print $2}}\' && '
           f'grep -cE "Mesh OK" log.checkMesh && '
           f'F=$(ls -t postProcessing/forces/*/coefficient*.dat | head -1) && '
           f'echo $F && grep -v "^#" $F | tail -1 && '
           f'grep -vc "^#" $F')
    out = [x for x in r.stdout.strip().splitlines() if x]
    if len(out) < 5:
        return None
    cells, meshok, coffile, lastrow, nrows = out[0], out[1], out[2], out[3], out[4]
    vals = lastrow.split()
    # tail mean over the last 30%
    r2 = sh(f'cd {c} && F=$(ls -t postProcessing/forces/*/coefficient*.dat | head -1) && '
            f'grep -v "^#" $F | awk -v n={nrows} \'NR>0.7*n {{c+=$5; d+=$2; k++}} '
            f'END {{printf "%.6f %.6f", c/k, d/k}}\'')
    tail = r2.stdout.split()
    return dict(cells=int(cells), mesh_ok=bool(int(meshok)), coeff_file=coffile,
                rows=int(nrows),
                Cl=float(vals[4]), Cd=float(vals[1]), Cm=float(vals[7]),
                Cl_tail=float(tail[0]), Cd_tail=float(tail[1]))


def wall_vtp(root, name, geom, tmp):
    c = f'/work/{root}/{name}'
    r = sh(f'cd {c} && foamToVTK -ascii -latestTime -fields p > log.vtk 2>&1; '
           f'find VTK -name "{STL[geom]}.vtp" -path "*boundary*" | head -1')
    path = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else ''
    if not path:
        return None, None
    tag = name.replace('.', 'p')
    dst = os.path.join(tmp, f'{root.replace("/", "_")}_{tag}.vtp')
    subprocess.run(['docker', 'cp', f'oftrans:{c}/{path}', dst],
                   capture_output=True)
    return dst, path


def eta_table(vtp, mach, halfspan=0.8, station=0.5, nbins=8):
    """Cp upper/lower and dCp at x/c=station, binned by eta=y/halfspan.

    Upper/lower is decided from the face normal, not from the sign of the face
    centroid z: on the cambered supercritical sections (dfvlr_r4, rae2822) the
    lower surface sits above z=0 over the aft part of the chord, so a z-sign
    test mislabels faces there.  In OpenFOAM the wall-patch normal points out of
    the fluid, i.e. into the wing, so the upper (suction) side has n_z < 0.
    """
    import numpy as np
    pts, conn, off, cd, pd = wallspan.read_vtp(vtp)
    ac, cx = wallspan.poly_area_normals(pts, conn, off)
    p = cd['p']
    nrm = ac / np.linalg.norm(ac, axis=1)[:, None]
    q = 0.5 * wallspan.RHO_INF * (mach * wallspan.A_INF) ** 2
    cp = (p - wallspan.PINF) / q
    base = (np.abs(cx[:, 0] - station) < 0.03) & (np.abs(nrm[:, 2]) > 0.5)
    up, lo = base & (nrm[:, 2] < 0), base & (nrm[:, 2] > 0)
    edges = np.linspace(0.0, halfspan, nbins + 1)
    yy = np.abs(cx[:, 1])
    idx = np.clip(np.digitize(yy, edges) - 1, 0, nbins - 1)
    rows = []
    for k in range(nbins):
        if (up & (idx == k)).sum() == 0 or (lo & (idx == k)).sum() == 0:
            continue
        cu = float(cp[up & (idx == k)].mean())
        cl_ = float(cp[lo & (idx == k)].mean())
        rows.append([round(float(0.5 * (edges[k] + edges[k + 1]) / halfspan), 3),
                     round(cu, 4), round(cl_, 4), round(cl_ - cu, 4)])
    return rows


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
    todo = [(n, g, m, t, 'ft') for n, g, m, t in CASES] + \
           [(os.path.basename(n), g, m, t, 'casc') for n, g, m, t in CONTROL]
    for name, geom, mach, tag, root in todo:
        st = case_stats(root, name)
        if st is None:
            print(f'MISSING {root}/{name}', file=sys.stderr)
            continue
        vtp, rel = wall_vtp(root, name, geom, tmp)
        rows = eta_table(vtp, mach) if vtp else []
        o = old.get((geom, round(mach, 3)), {})
        rec = dict(case=f'{root}/{name}', geometry=geom, mach=mach,
                   config=tag, **st,
                   wing_eta_dCp=rows,
                   campaign_cascade_Cl_rhoCentralFoam=o.get('of_cl'),
                   tsd_Cl=o.get('tsd_cl'))
        if o.get('of_cl') and st:
            rec['delta_vs_old_campaign_pct'] = round(
                100.0 * (st['Cl_tail'] - o['of_cl']) / o['of_cl'], 2)
        if o.get('tsd_cl') and st:
            rec['delta_vs_tsd_pct'] = round(
                100.0 * (st['Cl_tail'] - o['tsd_cl']) / o['tsd_cl'], 2)
            if o.get('of_cl'):
                rec['old_delta_vs_tsd_pct'] = round(o['err_pct'], 2)
        res.append(rec)
        print(f'{root}/{name:24s} Cl={st["Cl_tail"]:.5f} Cd={st["Cd_tail"]:.5f} '
              f'cells={st["cells"]}')

    json.dump(res, open(os.path.join(FT, 'freetip_results.json'), 'w'), indent=1)
    write_md(res)
    shutil.rmtree(tmp, ignore_errors=True)


def write_md(res):
    L = []
    A = L.append

    def fmt(v, spec='%+.1f%%'):
        return spec % v if isinstance(v, (int, float)) else '-'

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
    A('Solver: **rhoSimpleFoam** (steady). rhoCentralFoam is explicit and its '
      'time step is acoustic-CFL limited by the smallest near-wall cell; on the '
      "139.5k-cell free-tip mesh that measured **1.6 s/step on the container's "
      "2 cores**, i.e. ~5.5 h for the campaign's `t=0.25` per case, which is "
      'not affordable for a 7-case baseline. The cascade control case was '
      're-run with rhoSimpleFoam too, so free-tip vs cascade differs only in '
      'geometry. Prism layers are off: the wall is `slip` (inviscid), the '
      'layers resolve nothing, and snappy layer collapse at the newly exposed '
      'tip produced zero-normal-extent cells that killed every implicit solver '
      '(`deltaCoeffs -> inf`).\n')
    A('## Mesh\n')
    A('`blockMesh` 100x15x70 background (y cell 0.5333 - the campaign value - '
      'with grid lines landing exactly on the tips at `y=+-0.8`) + '
      '`snappyHexMesh` castellated/snap, surface level 3-4, distance regions '
      '(0.15:3, 0.4:2), wake box level 2. `checkMesh`: **Mesh OK** on every '
      'case, 136k-140k cells (the campaign cascade at the same refinement is '
      '49.7k). Wing patch extent `y in [-0.800, +0.800]` with 372 tip-cap '
      'faces at `|y|=0.8`; domain `x[-10,21] y[-4,4] z[-10,10]`; two boundary '
      'patches only (`farfield`, `naca0012.stl`).\n')
    A('## Results (steady, converged; CL = mean of the last 30% of iterations)\n')
    A('`CL vs campaign` uses the old `campaign/final_table.json` rhoCentralFoam '
      'value for the same geometry and Mach; `CL vs TSD` uses that table\'s '
      '`tsd_cl`, which the TSD side is concurrently re-baselining - treat the '
      'TSD column as indicative, not final.\n')
    A('| geometry | M | config | cells | CL | CD | CL vs campaign | CL vs TSD |')
    A('|---|---|---|---|---|---|---|---|')
    for r in res:
        n = r['case'].split('/')[-1].replace(f'{r["geometry"]}_', '')
        A(f'| {r["geometry"]} | {r["mach"]:.3f} | {r["config"]} | {r["cells"]} '
          f'| {r["Cl_tail"]:.5f} | {r["Cd_tail"]:.5f} | '
          f'{fmt(r.get("delta_vs_old_campaign_pct"))} | '
          f'{fmt(r.get("delta_vs_tsd_pct"))} |')
    A('')
    A('## Spanwise loading evidence (x/c = 0.5)\n')
    A('`dCp = Cp_lower - Cp_upper`, so `dCp` is the local loading at that '
      'station (Cp_upper from faces with `n_z < 0`, i.e. the suction side). '
      'A free tip must unload as `eta -> 1`; a tip clamped by a boundary must '
      'not vanish and typically turns *up* where the boundary pressure is '
      'imposed. `eta = |y|/0.8`, so `eta=0` is the wing centreline and '
      '`eta=1` is both tips. Caveat: on the supercritical sections at '
      'M >= 0.925 a shock sits at `x/c ~ 0.5`, so upper and lower Cp nearly '
      'cancel at that single station and `dCp(0.5)` is not a loading proxy '
      'there; the controlled comparison is the naca0012 free-tip vs cascade '
      'pair, which is the same geometry, the same Mach and the same solver.\n')
    for r in res:
        if not r['wing_eta_dCp']:
            continue
        A(f'\n`{r["case"]}`  (CL = {r["Cl_tail"]:.4f})\n')
        A('| eta | Cp upper | Cp lower | dCp |')
        A('|---|---|---|---|')
        for e, cu, cn, d in r['wing_eta_dCp']:
            A(f'| {e:.2f} | {cu:.3f} | {cn:.3f} | {d:+.4f} |')
    A('')
    A('Ratio `dCp(eta)/dCp(mid-span)`: naca0012 free tip 1.00 -> 0.52 at '
      'eta=0.81; naca0012 cascade 1.00 -> 0.80 at eta=0.81 and **1.96** at '
      'eta=0.94. Lifting-line (elliptic) loading has to vanish at a free tip; '
      'a boundary-clamped tip has no reason to. The outermost bin (eta=0.94) '
      'holds only a few faces on either side of the tip edge, so treat its '
      'absolute value as noisy - the direction of the two curves is the point.\n')
    A('## Verdict\n')
    A('**The configuration is now right, the CL agreement is not improved.**\n')
    A('1. The free-tip cases are genuine 3D wings: free tips 3.2 chords from '
      'any wall, 372 tip-cap faces present, no symmetry planes, `checkMesh` OK, '
      'and the local loading unloads toward the tip - naca0012 M0.80 dCp(0.5) '
      'falls from 0.0366 mid-span to 0.0190 at eta=0.81 (52% unloading) while '
      'the cascade control falls only 20% and then rises to 0.0678 (196%) in '
      'the outermost bin. That rise is the signature of the clamped tip, not of '
      'a periodic cascade.\n')
    A('2. The headline number moved the wrong way: naca0012 M0.80 was OF '
      '0.06658 vs TSD 0.05000 (-24.9%); it is now free-tip 0.05440 (+8.8%), '
      'while the *cascade* configuration re-run under the same solver gives '
      '0.04938 (-1.2%). So the old campaign number agreed with the TSD for the '
      'wrong reason. The free-tip CL 0.0544 is, however, the value the physics '
      'predicts for AR 1.6 at alpha=1.25 deg (lifting line with a0=2*pi: '
      'CL ~ 0.061), so it is the cascade number that should be regarded with '
      'suspicion.\n')
    A('3. The dominant uncertainty is numerical, not geometric. The same '
      'cascade mesh and BCs give CL 0.0494 (rhoSimpleFoam) vs 0.0666 '
      '(rhoCentralFoam) at M0.80 - a 26% spread from the solver/scheme alone. '
      'No geometry conclusion smaller than that is meaningful until a '
      'solver/scheme sensitivity study exists.\n')
    A('4. The spanwise loading is washed out: mid-span dCp(0.5) = 0.037 vs '
      '0.087 for an isolated 2D flat plate at the same incidence, because the '
      'background y cell is 0.533 and the surface only reaches ~24 spanwise '
      'cells. A span refinement study (and a tip-focused refinement region) is '
      'required before quoting a CL band.\n')
    A('5. Sensitivity to the tip clamp is geometry dependent: dfvlr_r4 at '
      'M0.85 moves only +1.7% (0.2674 -> 0.2719), so its campaign number was '
      'barely affected by the clamped tips, while rae2822 moves -31.7% / '
      '-38.2% at M0.85 / 0.925 (0.2510 -> 0.1714, 0.2778 -> 0.1717). Any '
      'conclusion drawn from the rae2822 campaign points in particular is '
      'invalid for a free-tip comparison. The free-tip rae2822 CL is now '
      'nearly Mach independent over 0.85-0.925 (0.1714, 0.1717) with CD '
      'climbing 0.035 -> 0.051, which is the expected shock-drag behaviour for '
      'a section already above its drag-rise Mach.\n')
    A('**Usable?** Yes as the comparison configuration - it is the only '
      'OpenFOAM side that matches the TSD geometry. Not yet usable as a CL '
      'anchor within better than ~+-15% until items 3 and 4 are closed.\n')
    A('## Reproduce\n')
    A('```\n'
      'python cases/openfoam/transonic/freetip/tools/build_freetip_cases.py\n'
      '# then, in container oftrans (OpenFOAM 2406), per case dir:\n'
      'blockMesh && snappyHexMesh -overwrite && checkMesh && rhoSimpleFoam\n'
      '# FOAM_SIGFPE=false may be needed on some builds (spurious trap in\n'
      '# basicFvGeometryScheme::deltaCoeffs)\n'
      'python cases/openfoam/transonic/freetip/tools/make_results.py\n'
      '```\n')
    return open(os.path.join(FT, 'RESULTS.md'), 'w').write('\n'.join(L))


if __name__ == '__main__':
    main()
