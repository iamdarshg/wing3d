"""of_mesh -- permanent body-fitted mesh pathway for the viscous anchor leg.

Orchestrates, on the host, driving OpenFOAM in the `oftrans` container:

    input.nml  ->  hcf_wing (NASA TMR structured O-grid, y+ target)
               ->  plot3dToFoam -noBlank
               ->  autoPatch 30 -overwrite
               ->  rename patches by geometry (wing / symmetry / farfield)
               ->  write a rhoSimpleFoam + Spalart-Allmaras case skeleton

See tools/of_mesh/README.md for why this generator and the two conversion
gotchas (formatted .p3d, -noBlank).

Usage
-----
  # NACA00XX section, AR from the semi-span
  python tests3d/of_mesh.py --airfoil naca0012 --semispan 1.6 \\
      --nspan 32 --nradial 64 --yplus 1 --re 6.0e6 \\
      --mach 0.85 --alpha 1.25 --out /work/cases/naca_ar16

  # one of our campaign sections (written out as a discrete .dat first)
  python tests3d/of_mesh.py --airfoil dfvlr_r4 --semispan 1.6 ...

Notes
-----
* semi-span b is in ROOT-CHORD units (root chord is normalised to 1).
* the grid is a half wing with a symmetry plane at the root (y=0), which is
  the TMR convention.  For a full free-tip wing use --full-wing later (not
  implemented; the half-model is what the M6/NACA0012 references use).
* y+ is a target only; verify with `checkMesh`/`yPlus` after a solve.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

CONTAINER = 'oftrans'
GEN_DIR = '/work/tmr/wing_release/wing_release_072319'
GEN_BIN = GEN_DIR + '/hcf_wing'


def cp(*parts):
    """POSIX join for container paths."""
    return '/'.join(str(x).rstrip('/') for x in parts)


# ---------------------------------------------------------------- container
def dexec(cmd, check=False):
    r = subprocess.run(['docker', 'exec', CONTAINER, 'bash', '-lc', cmd],
                       capture_output=True, text=True)
    if check and r.returncode != 0:
        raise RuntimeError('docker exec failed: %s\n%s\n%s'
                           % (cmd, r.stdout, r.stderr))
    return r


def push(local, remote):
    subprocess.run(['docker', 'cp', local, '%s:%s' % (CONTAINER, remote)],
                   check=True)


def pull(remote, local):
    subprocess.run(['docker', 'cp', '%s:%s' % (CONTAINER, remote), local],
                   check=True)


def write_remote(remote, text):
    fd, tmp = tempfile.mkstemp(suffix='.txt')
    with os.fdopen(fd, 'w') as f:
        f.write(text)
    push(tmp, remote)
    os.unlink(tmp)


def read_remote(remote):
    fd, tmp = tempfile.mkstemp(suffix='.txt')
    os.close(fd)
    pull(remote, tmp)
    with open(tmp, errors='ignore') as f:
        s = f.read()
    os.unlink(tmp)
    return s


# ------------------------------------------------------------------- inputs
def section_to_dat(geo, path):
    """Write one of our campaign sections as a Selig/PLOT3D-style .dat."""
    sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
    from gen_transonic_stl import (naca0012_xy, rae2822_xy,
                                   biconvex_xy, dfvlr_r4_xy)
    import numpy as np
    if geo == 'naca0012':
        x, yt = naca0012_xy(400)
        xu, yu, xl, yl = x, yt, x, -yt
    elif geo == 'biconvex':
        x, yb = biconvex_xy(400)
        xu, yu, xl, yl = x, yb, x, -yb
    elif geo == 'rae2822':
        xu, yu, xl, yl = rae2822_xy()
    elif geo == 'dfvlr_r4':
        xu, yu, xl, yl = dfvlr_r4_xy()
    else:
        raise SystemExit('unknown section %s' % geo)
    xu, yu, xl, yl = map(np.asarray, (xu, yu, xl, yl))
    if xl[0] > xl[-1]:
        xl, yl = xl[::-1], yl[::-1]
    if xu[0] > xu[-1]:
        xu, yu = xu[::-1], yu[::-1]
    lines = []
    # hcf_wing wants a single closed loop: upper TE->LE then lower LE->TE
    for x, y in zip(xu[::-1], yu[::-1]):
        lines.append('%.8f %.8f' % (x, y))
    for x, y in zip(xl, yl):
        lines.append('%.8f %.8f' % (x, y))
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return path


def is_naca(name):
    """'naca0012' -> True (generator builds NACA00XX from its 2-digit thickness)."""
    return name.startswith('naca') and name[-2:].isdigit() and len(name) <= 8


def make_nml(a):
    """Build the hcf_wing &input_parameters namelist."""
    if is_naca(a.airfoil):
        af_type, af_line = 1, 'nacaxx = %s' % a.airfoil[4:]
    else:
        af_type, af_line = 0, 'airfoil_data_file = "%s"' % a.dat_name
    return """&input_parameters
          target_reynolds_number = %.4e
                   target_y_plus = %.2f
                      igrid_type = 5
                    airfoil_type = %d
               %s
                        tip_type = 'round'
                     use_s_blunt = T
       tip_corner_normal_deflect = 0.0
                              tr = %.6f
                            beta = %.6f
                         R_outer = %.1f
           nnodes_cylinder_input = %d
                               b = %.6f
                           nr_gs = 8
                             nre = %d
                       wing_side = 2
                       root_le_x = 0.0
    stretching_tanh_towards_lete = 4.5
           hc_region_fixed_outer = 25.0
          stretching_towards_tip = -1
             generate_ugrid_file = T
          ugrid_file_unformatted = F
           generate_p3d_ufmt_file = T
            generate_su2grid_file = T
                 generate_k_file = T
               k_file_unformatted = F
           generate_line_file_nc = T
           generate_line_file_cc = T
             generate_tec_file_b = T
             generate_tec_file_v = F
         sf_grow_radial_outer_c1 = 6.0
         sf_grow_radial_outer_c2 = 3.0
         sf_grow_towards_tip_c1  = 5.0
         sf_grow_towards_tip_c2  = 2.5
               use_s_blunt_match = T
             use_s_blunt_ssf   = 3.0
             use_s_blunt_ssfyz = 3.0
             check_non_planarity = F
/
""" % (a.re, a.yplus, af_type, af_line, a.taper, a.sweep, a.r_outer,
       a.nspan, a.semispan, a.nradial)


# ------------------------------------------------------------------ patches
def patch_bboxes(points_text, faces_text, boundary_text):
    import numpy as np
    pts, started = [], False
    for line in points_text.splitlines():
        if not started:
            if line.strip() == '(':
                started = True
            continue
        m = re.match(r'\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)',
                     line)
        if m:
            pts.append([float(m.group(1)), float(m.group(2)),
                        float(m.group(3))])
    P = np.array(pts)
    F, started = [], False
    for line in faces_text.splitlines():
        if not started:
            if line.strip() == '(':
                started = True
            continue
        m = re.match(r'\s*(\d+)\(([\d\s]+)\)', line)
        if m:
            F.append([int(t) for t in m.group(2).split()])
    out = {}
    for name, body in re.findall(r'(\w+)\s*\{([^}]*)\}', boundary_text):
        nF = re.search(r'nFaces\s+(\d+)', body)
        sF = re.search(r'startFace\s+(\d+)', body)
        if not (nF and sF) or int(nF.group(1)) == 0:
            continue
        n0 = int(sF.group(1))
        n1 = n0 + int(nF.group(1))
        idx = set()
        for k in range(n0, min(n1, len(F))):
            idx.update(F[k])
        A = P[sorted(idx)]
        out[name] = (A[:, 0].min(), A[:, 0].max(), A[:, 1].min(),
                     A[:, 1].max(), A[:, 2].min(), A[:, 2].max())
    return out


def classify(bbox, semispan):
    x0, x1, y0, y1, z0, z1 = bbox
    if abs(y0) < 1e-9 and abs(y1) < 1e-9:
        return 'symmetry'
    if (z1 - z0) < 0.15 and x1 < 3.0 * semispan and y1 < 1.5 * semispan:
        return 'wing'
    return 'farfield'


def patch_bboxes_skip_empty(boxes, boundary_text):
    """Also return the names of zero-face patches (defaultFaces etc.)."""
    empty = []
    for name, body in re.findall(r'(\w+)\s*\{([^}]*)\}', boundary_text):
        nF = re.search(r'nFaces\s+(\d+)', body)
        if nF and int(nF.group(1)) == 0:
            empty.append(name)
    return empty


def rewrite_boundary(boundary_text, mapping):
    """Rename auto* patches to wing/symmetry/farfield with proper types."""
    def ptype(name):
        if name.startswith('wing'):
            return 'wall'
        if name.startswith('symmetry'):
            return 'symmetry'
        return 'patch'
    txt = boundary_text
    for old, new in mapping.items():
        txt = re.sub(r'(?m)^    %s\s*$' % re.escape(old), '    ' + new, txt)
        txt = re.sub(r'(%s\s*\{\s*\n\s*type\s+)\w+' % re.escape(new),
                     lambda m: m.group(1) + ptype(new), txt)
    return txt


# -------------------------------------------------------------------- cases
def field_files(a, patches):
    """0/ fields for rhoSimpleFoam + SpalartAllmaras."""
    R = 296.803          # N2
    p, T = 1.0e5, 300.0
    rho = p / (R * T)
    import math
    gam = 1.4
    aa = math.sqrt(gam * R * T)
    U = a.mach * aa
    al = math.radians(a.alpha)
    ux, uz = U * math.cos(al), U * math.sin(al)
    nut = 3.0 * 1.5e-5   # ~3 nu
    far = ' '.join(sorted(patches['farfield']))
    sym = ' '.join(sorted(patches['symmetry']))
    wall = ' '.join(sorted(patches['wing']))

    def hdr(cls, obj):
        return ('FoamFile { version 2.0; format ascii; class %s; object %s; }\n'
                % (cls, obj))

    U_f = hdr('volVectorField', 'U') + (
        'dimensions [0 1 -1 0 0 0 0];\n'
        'internalField uniform (%f 0 %f);\n'
        'boundaryField {\n' % (ux, uz)
        + ''.join(' %s { type freestream; freestreamValue uniform (%f 0 %f); }\n'
                  % (b, ux, uz) for b in patches['farfield'])
        + ''.join(' %s { type symmetry; }\n' % b for b in patches['symmetry'])
        + ''.join(' %s { type noSlip; }\n' % b for b in patches['wing'])
        + '}\n')
    p_f = hdr('volScalarField', 'p') + (
        'dimensions [1 -1 -2 0 0 0 0];\n'
        'internalField uniform %f;\n'
        'boundaryField {\n' % p
        + ''.join(' %s { type freestreamPressure; freestreamValue uniform %f; }\n'
                  % (b, p) for b in patches['farfield'])
        + ''.join(' %s { type symmetry; }\n' % b for b in patches['symmetry'])
        + ''.join(' %s { type zeroGradient; }\n' % b for b in patches['wing'])
        + '}\n')
    T_f = hdr('volScalarField', 'T') + (
        'dimensions [0 0 0 1 0 0 0];\n'
        'internalField uniform %f;\n'
        'boundaryField {\n' % T
        + ''.join(' %s { type freestream; freestreamValue uniform %f; }\n'
                  % (b, T) for b in patches['farfield'])
        + ''.join(' %s { type symmetry; }\n' % b for b in patches['symmetry'])
        + ''.join(' %s { type zeroGradient; }\n' % b for b in patches['wing'])
        + '}\n')
    nut_f = hdr('volScalarField', 'nut') + (
        'dimensions [0 2 -1 0 0 0 0];\n'
        'internalField uniform %e;\n'
        'boundaryField {\n' % nut
        + ''.join(' %s { type calculated; value uniform %e; }\n' % (b, nut)
                  for b in patches['farfield'])
        + ''.join(' %s { type symmetry; }\n' % b for b in patches['symmetry'])
        + ''.join(' %s { type fixedValue; value uniform 0; }\n' % b
                  for b in patches['wing'])
        + '}\n')
    nuT_f = hdr('volScalarField', 'nuTilda') + (
        'dimensions [0 2 -1 0 0 0 0];\n'
        'internalField uniform %e;\n'
        'boundaryField {\n' % (4.0 * nut)
        + ''.join(' %s { type inletOutlet; inletValue uniform %e; value uniform %e; }\n'
                  % (b, 4.0 * nut, 4.0 * nut) for b in patches['farfield'])
        + ''.join(' %s { type symmetry; }\n' % b for b in patches['symmetry'])
        + ''.join(' %s { type fixedValue; value uniform 0; }\n' % b
                  for b in patches['wing'])
        + '}\n')
    alpha_t = hdr('volScalarField', 'alphat') + (
        'dimensions [1 -1 -1 0 0 0 0];\n'
        'internalField uniform 0;\n'
        'boundaryField {\n'
        + ''.join(' %s { type calculated; value uniform 0; }\n' % b
                  for b in patches['farfield'])
        + ''.join(' %s { type symmetry; }\n' % b for b in patches['symmetry'])
        + ''.join(' %s { type compressible::alphatWallFunction; Prt 0.85; value uniform 0; }\n'
                  % b for b in patches['wing'])
        + '}\n')
    return {'U': U_f, 'p': p_f, 'T': T_f, 'nut': nut_f,
            'nuTilda': nuT_f, 'alphat': alpha_t}


def system_files(a):
    hdr = lambda o: ('FoamFile { version 2.0; format ascii; class dictionary; '
                     'object %s; }\n' % o)
    control = hdr('controlDict') + """application     rhoSimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         3000;
deltaT          1;
writeControl    timeStep;
writeInterval   500;
purgeWrite      0;
writeFormat     ascii;
writePrecision  10;
runTimeModifiable true;
functions
{
    forces
    {
        type            forces;
        libs            ("libforces.so");
        patches         (".*");
        rho             rhoInf;
        rhoInf          1.0;
        CofR            (0 0 0);
        log             true;
    }
}
"""
    fvSchemes = hdr('fvSchemes') + """ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; grad(p) Gauss linear; }
divSchemes {
 default none;
 div(phi,U) bounded Gauss linearUpwind grad(U);
 div(phi,K) bounded Gauss linearUpwind grad(U);
 div(phi,Ekp) bounded Gauss linearUpwind grad(U);
 div(phi,e) bounded Gauss upwind;
 div(phi,T) bounded Gauss upwind;
 div(phi,nuTilda) bounded Gauss upwind;
 div((nuEff*dev2(T(grad(U))))) Gauss linear;
 div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
fluxRequired { default no; p; }
wallDist { method meshWave; }
"""
    fvSolution = hdr('fvSolution') + """solvers
{
 p { solver GAMG; tolerance 1e-7; relTol 0.05; smoother DICGaussSeidel;
     nCellsInCoarsestLevel 1000; }
 U { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-8; relTol 0.05; }
 e { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-8; relTol 0.05; }
 T { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-8; relTol 0.05; }
 nuTilda { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-8; relTol 0.05; }
}
rhoSimple { momentumPredictor yes; nNonOrthogonalCorrectors 2; consistent no;
    rhoMin 0.05; rhoMax 8; }
relaxationFactors {
 fields { p 0.2; rho 0.3; }
 equations { U 0.5; e 0.7; nuTilda 0.7; T 0.7; }
}
"""
    mom = hdr('turbulenceProperties') + """simulationType  RAS;
RAS
{
    model           SpalartAllmaras;
    turbulence      on;
    printCoeffs     on;
}
"""
    therm = hdr('thermophysicalProperties') + """thermoType
{
    type            hePsiThermo;
    mixture         pureMixture;
    transport       sutherland;
    thermo          janaf;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleInternalEnergy;
}
mixture
{
    specie
    {
        molWeight       28.01348;
    }
    thermodynamics
    {
        Tlow            200;
        Thigh           6000;
        Tcommon         1000;
        highCpCoeffs    (3.6744 0.0005285 -4.8e-08 2.2e-12 0 -1104 0);
        lowCpCoeffs     (3.5310 -0.0001236 2.65e-07 -1.4e-11 0 -1046 0);
    }
    transport
    {
        As              1.458e-06;
        Ts              110.4;
    }
}
"""
    return {'controlDict': control, 'fvSchemes': fvSchemes,
            'fvSolution': fvSolution, 'turbulenceProperties': mom,
            'thermophysicalProperties': therm}


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--airfoil', default='naca0012')
    ap.add_argument('--semispan', type=float, default=1.6)
    ap.add_argument('--taper', type=float, default=1.0)
    ap.add_argument('--sweep', type=float, default=0.0)
    ap.add_argument('--nspan', type=int, default=32)
    ap.add_argument('--nradial', type=int, default=64)
    ap.add_argument('--r-outer', type=float, default=25.0)
    ap.add_argument('--yplus', type=float, default=1.0)
    ap.add_argument('--re', type=float, default=6.0e6)
    ap.add_argument('--mach', type=float, default=0.85)
    ap.add_argument('--alpha', type=float, default=1.25)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    a.dat_name = 'section.dat'

    # 0. sanity
    if dexec('test -x %s && echo ok' % GEN_BIN).stdout.strip() != 'ok':
        raise SystemExit('hcf_wing not found at %s -- see README' % GEN_BIN)

    # 1. namelist (+ discrete section file)
    write_remote(cp(GEN_DIR, 'input.nml'), make_nml(a))
    if not is_naca(a.airfoil):
        fd, tmp = tempfile.mkstemp(suffix='.dat')
        os.close(fd)
        section_to_dat(a.airfoil, tmp)
        push(tmp, cp(GEN_DIR, a.dat_name))
        os.unlink(tmp)

    # 2. generate
    r = dexec('cd %s && ./hcf_wing' % GEN_DIR)
    if 'successfully completed' not in r.stdout:
        print(r.stdout[-2000:], r.stderr[-2000:])
        raise SystemExit('hcf_wing failed')
    print('grid generated')

    # 3. convert
    out = a.out.rstrip('/')
    dexec('rm -rf %s && mkdir -p %s/system %s/constant' % (out, out, out),
          check=True)
    dexec('cd %s && cp %s/*.p3d . && ls *.p3d' % (out, GEN_DIR), check=True)
    write_remote(cp(out, 'system/controlDict'),
                 system_files(a)['controlDict'])
    r = dexec('source /usr/lib/openfoam/openfoam2406/etc/bashrc && cd %s && '
              'plot3dToFoam -noBlank *.p3d' % out)
    if 'Writing polyMesh' not in r.stdout:
        print(r.stdout[-2000:], r.stderr[-1000:])
        raise SystemExit('plot3dToFoam failed')
    dexec('source /usr/lib/openfoam/openfoam2406/etc/bashrc && cd %s && '
          'autoPatch 30 -overwrite' % out, check=True)
    print('mesh converted')

    # 4. classify + rename patches
    pts = read_remote(cp(out, 'constant/polyMesh/points'))
    fcs = read_remote(cp(out, 'constant/polyMesh/faces'))
    bnd = read_remote(cp(out, 'constant/polyMesh/boundary'))
    boxes = patch_bboxes(pts, fcs, bnd)
    mapping, groups = {}, {'wing': [], 'symmetry': [], 'farfield': []}
    for name, bb in boxes.items():
        kind = classify(bb, a.semispan)
        mapping[name] = kind
        groups[kind].append(name)
        print('  %-10s -> %-9s bbox x[%.2f,%.2f] y[%.2f,%.2f] z[%.3f,%.3f]'
              % (name, kind, *bb))
    # zero-face patches still need a BC entry: fold them into farfield
    for name in patch_bboxes_skip_empty(boxes, bnd):
        mapping[name] = 'farfield'
        groups['farfield'].append(name)
        print('  %-10s -> farfield (empty patch, folded)' % name)
    if not groups['wing']:
        raise SystemExit('could not identify the wing patch')
    # avoid name collisions if two patches map to the same kind
    final = {}
    for kind, names in groups.items():
        for i, n in enumerate(names):
            final[n] = kind if i == 0 else '%s%d' % (kind, i)
    # field BCs must use the FINAL (renamed) patch names
    groups = {k: [final[n] for n in v] for k, v in groups.items()}
    write_remote(cp(out, 'constant/polyMesh/boundary'),
                 rewrite_boundary(bnd, final))

    # 5. case skeleton
    sysf = system_files(a)
    for k, v in sysf.items():
        if k == 'controlDict':
            continue
        write_remote(cp(out, 'system', k), v)
    for k in ('turbulenceProperties', 'thermophysicalProperties'):
        write_remote(cp(out, 'constant', k), sysf[k])
    fields = field_files(a, groups)
    dexec('mkdir -p %s/0' % out, check=True)
    for k, v in fields.items():
        write_remote(cp(out, '0', k), v)

    print('\ncase ready: %s  (rhoSimpleFoam + SpalartAllmaras)' % out)
    print('run: docker exec %s bash -lc "source '
          '/usr/lib/openfoam/openfoam2406/etc/bashrc && cd %s && '
          'rhoSimpleFoam > log.run 2>&1 &"' % (CONTAINER, out))


if __name__ == '__main__':
    main()