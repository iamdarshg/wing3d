"""Build free-tip (isolated, AR 1.6) OpenFOAM cases for the TSD re-baseline.

Why: the `campaign/` cases put `symmetry` planes on both y-faces of a box whose
y extent is exactly the wing span (1.6), i.e. the wing spans wall-to-wall.  That
is a *cascade* with period 1.6 (solidity 0.625), not an isolated wing, and it
suppresses lift ~5x (NACA0012 M0.80: Cl 0.067 vs 0.352 isolated-2D).  The TSD
side is a free-tip AR 1.6 wing, so every campaign-vs-TSD error number compares
two different configurations.

This builder keeps every numerics/BC/thermo convention of the campaign case
byte-identical and changes ONLY the geometry of the box:

  * background box y in [-4.0, +4.0] instead of [-0.8, +0.8]
  * all six outer faces are `farfield` (permeable `freestream`), no `symmetry`
  * wing surface y in [-0.8, +0.8] (unchanged STL), tips clear of every wall
    (3.2 chords of free air outboard of each tip)
  * y-cell size kept at 2.4/9 = 0.5333 = the campaign value, with the y
    gridlines landing exactly on the tips (+-0.8) so snappy snaps a clean tip

Everything else (rhoCentralFoam, localEuler, Kurganov, vanLeer, N2 janaf,
p=1e5, T=300, slip wall, Aref=c*span=1.6, liftDir (0,0,1), CofR (0.25,0,0))
is copied from cases/openfoam/transonic/campaign/<case>/.

Run:  python cases/openfoam/transonic/freetip/tools/build_freetip_cases.py [outdir] [names...]
"""
import math
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FT = os.path.dirname(HERE)                       # cases/openfoam/transonic/freetip
TRANS = os.path.dirname(FT)                      # cases/openfoam/transonic
CAMPAIGN = os.path.join(TRANS, 'campaign')
STLDIR = os.path.join(TRANS, 'stl')

FOAM = 'FoamFile { version 2.0; format ascii;'

# --- geometry / gas constants, campaign-identical ---------------------------
TINF = 300.0
PINF = 1.0e5
MOLW = 28.01348                                  # N2, campaign value
RU = 8314.462618
RS = RU / MOLW
GAMMA = 1.4
A_INF = math.sqrt(GAMMA * RS * TINF)             # ~353.06 m/s
RHO_INF = PINF / (RS * TINF)                     # ~1.1231 kg/m3
ALPHA = 1.25
SPAN = 1.6                                       # wing y in [-0.8, +0.8]
AREF = 1.6                                       # c * span
YBOX = 4.0                                       # domain y in [-4, +4]
NY = 15                                          # 15 * 0.5333 = 8.0, tips on gridlines

# --- background box, campaign-identical in x/z ------------------------------
X0, X1 = -10.0, 21.0
Z0, Z1 = -10.0, 10.0
NX, NZ = 100, 70

# --- snappy, campaign-identical levels/regions ------------------------------
MAX_LOCAL = 140000
MAX_GLOBAL = 340000
# This is an inviscid (slip-wall) calculation: prism layers resolve nothing
# here and snappy's layer collapse at the now-exposed wing tip produced
# zero-normal-extent cells (deltaCoeffs -> inf) that killed every implicit
# solver.  Layers off.
NSURFACELAYERS = 0


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write(text)


def blockmesh():
    return f"""{FOAM} class dictionary; object blockMeshDict; }}
convertToMeters 1;
vertices (
 ({X0} {-YBOX} {Z0}) ({X1} {-YBOX} {Z0}) ({X1} {YBOX} {Z0}) ({X0} {YBOX} {Z0})
 ({X0} {-YBOX} {Z1}) ({X1} {-YBOX} {Z1}) ({X1} {YBOX} {Z1}) ({X0} {YBOX} {Z1}));
blocks (hex (0 1 2 3 4 5 6 7) ({NX} {NY} {NZ}) simpleGrading (1 1 1));
edges ();
boundary (
 farfield {{ type patch; faces ((0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }});
mergePatchPairs ();
"""


def snappy(stl):
    layers = 'true' if NSURFACELAYERS else 'false'
    return f"""{FOAM} class dictionary; object snappyHexMeshDict; }}
castellatedMesh true; snap true; addLayers {layers};
geometry {{ {stl} {{ type triSurfaceMesh; file "{stl}"; }}
 wakebox {{ type searchableBox; min (0 -0.8 -0.6); max (4 0.8 0.6); }} }}
castellatedMeshControls {{
 maxLocalCells {MAX_LOCAL}; maxGlobalCells {MAX_GLOBAL}; minRefinementCells 10;
 maxLoadUnbalance 0.1; nCellsBetweenLevels 3;
 features ();
 refinementSurfaces {{ {stl} {{ level (3 4); patchInfo {{ type wall; }} }} }}
 resolveFeatureAngle 25;
 refinementRegions {{
  {stl} {{ mode distance; levels ((0.15 3) (0.4 2)); }}
  wakebox {{ mode inside; levels ((1E15 2)); }} }}
 locationInMesh (5.5 0 0);
 allowFreeStandingZoneFaces true;
}}
snapControls {{ nSmoothPatch 3; tolerance 2.0; nSolveIter 30; nRelaxIter 5; }}
addLayersControls {{
 relativeSizes false;
 layers {{ "{stl}" {{ nSurfaceLayers {NSURFACELAYERS}; }} }}
 expansionRatio 1.2; finalLayerThickness 2e-4; minThickness 1e-5;
 nGrow 0; featureAngle 60; slipFeatureAngle 30; nRelaxIter 5;
 nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10;
 maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3;
 minMedianAxisAngle 90; minMedialAxisAngle 90; nBufferCellsNoExtrude 0;
 nLayerIter 50; }}
meshQualityControls {{ #include "meshQualityDict" nSmoothScale 4; errorReduction 0.75; }}
writeFlags (scalarLevels layerSets);
mergeTolerance 1e-6;
"""


def control(stl, M, end_t=0.25):
    magU = M * A_INF
    a = math.radians(ALPHA)
    return f"""{FOAM} class dictionary; object controlDict; }}
application rhoCentralFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {end_t};
deltaT 2e-5; writeControl adjustable; writeInterval 0.05;
purgeWrite 2; writeFormat ascii; writePrecision 6; writeCompression off;
timeFormat general; timePrecision 6; runTimeModifiable true;
adjustTimeStep yes; maxCo 0.5; maxDeltaT 1e-4;
functions {{
 forces {{
  type forceCoeffs; libs ("libforces.so");
  writeControl timeStep; writeInterval 10;
  patches ("{stl}");
  rho rhoInf; rhoInf {RHO_INF:.6f}; CofR (0.25 0 0);
  liftDir (0 0 1); dragDir (1 0 0); pitchAxis (0 1 0);
  magUInf {magU:.6f}; lRef 1.0; Aref {AREF};
 }}
 residuals {{ type residuals; libs ("libutilityFunctionObjects.so");
  writeControl timeStep; writeInterval 10; fields (p U T); }}
}}
"""


# --------------------------------------------------------------------------
# steady (rhoSimpleFoam) variants.  Chosen as the primary solver for the
# free-tip re-baseline: rhoCentralFoam is explicit, so its time step is
# acoustic-CFL limited by the smallest near-wall cell; on the 139.5k-cell
# free-tip mesh it costs ~1.6 s/step on the container's 2 cores, i.e. ~5.5 h
# for the campaign's t=0.25.  Everything else (thermo, BCs, Aref, liftDir) is
# campaign-identical -> the steady solve is a numerics-only difference, and it
# is validated against the campaign cascade level in RESULTS.md.
# --------------------------------------------------------------------------

def fvschemes_simple():
    return f"""{FOAM} class dictionary; object fvSchemes; }}
ddtSchemes {{ default steadyState; }}
gradSchemes {{ default Gauss linear; grad(p) Gauss linear; }}
divSchemes {{
 default Gauss linear;
 div(phi,U)      Gauss linearUpwind grad(U);
 div(phi,T)      Gauss upwind;
 div(phi,e)      Gauss upwind;
 div(phi,h)      Gauss upwind;
 div(phi,K)      Gauss upwind;
 div(phi,k)      Gauss upwind;
 div(phi,omega)  Gauss upwind;
 div(phi,nuTilda) Gauss upwind;
 div(phid,p)     Gauss upwind;
 div((muEff*dev2(T(grad(U))))) Gauss linear;
 div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}}
laplacianSchemes {{ default Gauss linear corrected; }}
interpolationSchemes {{ default linear; }}
snGradSchemes {{ default corrected; }}
fluxRequired {{ default no; p; }}
"""


def fvsolution_simple():
    body = (
        " momentumPredictor no;\n"
        " nNonOrthogonalCorrectors 1;\n"
        " consistent yes;\n"
        " transonic yes;\n"
    )
    head = (
        FOAM + " class dictionary; object fvSolution; }\n"
        "solvers\n"
        "{\n"
        ' "rho.*" { solver diagonal; }\n'
        " p { solver GAMG; smoother GaussSeidel; tolerance 1e-8; relTol 0.01;\n"
        "     nCellsInCoarsestLevel 1000; }\n"
        ' "(U|e|h|k|omega)" { solver smoothSolver; smoother GaussSeidel; nSweeps 2;\n'
        "     tolerance 1e-9; relTol 0.01; }\n"
        "}\n"
        "// rhoSimpleFoam in 2406 reads the SIMPLE dict (solutionControl); the\n"
        "// rhoSimple alias is given identical content so either name is harmless.\n"
        "SIMPLE\n{\n" + body + "}\n"
        "rhoSimple\n{\n" + body + "}\n"
        "relaxationFactors\n"
        "{\n"
        " fields { p 0.3; rho 0.3; }\n"
        " equations { U 0.6; e 0.6; }\n"
        "}\n"
    )
    return head


def control_simple(stl, M, end_iter=6000):
    magU = M * A_INF
    return f"""{FOAM} class dictionary; object controlDict; }}
application rhoSimpleFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {end_iter};
deltaT 1; writeControl timeStep; writeInterval 500; purgeWrite 0;
writeFormat ascii; writePrecision 8; writeCompression off;
timeFormat general; timePrecision 8; runTimeModifiable true;
functions {{
 forces {{
  type forceCoeffs; libs ("libforces.so");
  writeControl timeStep; writeInterval 1;
  patches ("{stl}");
  rho rhoInf; rhoInf {RHO_INF:.6f}; CofR (0.25 0 0);
  liftDir (0 0 1); dragDir (1 0 0); pitchAxis (0 1 0);
  magUInf {magU:.6f}; lRef 1.0; Aref {AREF};
 }}
 residuals {{ type residuals; libs ("libutilityFunctionObjects.so");
  writeControl timeStep; writeInterval 1; fields (p U e); }}
}}
"""



def field_p(stl):
    return f"""{FOAM} class volScalarField; object p; }}
dimensions [1 -1 -2 0 0 0 0];
internalField uniform {PINF};
boundaryField {{
 farfield {{ type freestreamPressure; freestreamValue uniform {PINF}; }}
 "{stl}" {{ type zeroGradient; }}
}}
"""


def field_u(stl, M):
    magU = M * A_INF
    a = math.radians(ALPHA)
    ux, uz = magU * math.cos(a), magU * math.sin(a)
    return f"""{FOAM} class volVectorField; object U; }}
dimensions [0 1 -1 0 0 0 0];
internalField uniform ({ux:.6f} 0 {uz:.6f});
boundaryField {{
 farfield {{ type freestream; freestreamValue uniform ({ux:.6f} 0 {uz:.6f}); }}
 "{stl}" {{ type slip; }}
}}
"""


def field_t(stl):
    return f"""{FOAM} class volScalarField; object T; }}
dimensions [0 0 0 1 0 0 0];
internalField uniform {TINF};
boundaryField {{
 farfield {{ type freestream; freestreamValue uniform {TINF}; }}
 "{stl}" {{ type zeroGradient; }}
}}
"""


CASES = [
    ('naca0012_M0.800_a1.25', 'naca0012.stl', 0.800),
    ('dfvlr_r4_M0.850_a1.25', 'dfvlr_r4.stl', 0.850),
    ('dfvlr_r4_M0.925_a1.25', 'dfvlr_r4.stl', 0.925),
    ('dfvlr_r4_M1.000_a1.25', 'dfvlr_r4.stl', 1.000),
    ('rae2822_M0.850_a1.25', 'rae2822.stl', 0.850),
    ('rae2822_M0.925_a1.25', 'rae2822.stl', 0.925),
    ('rae2822_M1.000_a1.25', 'rae2822.stl', 1.000),
]

# Control case: the campaign (wing-tip-clamped) geometry re-run with exactly the
# same solver/scheme settings as the free-tip cases, so "free tip vs campaign"
# differs only in the box.  Only system/{fvSchemes,fvSolution,controlDict} are
# replaced; blockMeshDict / snappyHexMeshDict / 0/* / constant/* are the
# campaign template verbatim (that is the configuration under test).
CONTROL = ('naca0012_M0.800_a1.25', 'naca0012.stl', 0.800)


def build_control(out):
    name, stl, M = CONTROL
    src = os.path.join(CAMPAIGN, name)
    d = os.path.join(out, 'cascade_control_' + name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    shutil.copytree(src, d)
    write(os.path.join(d, 'system', 'fvSchemes'), fvschemes_simple())
    write(os.path.join(d, 'system', 'fvSolution'), fvsolution_simple())
    write(os.path.join(d, 'system', 'controlDict'), control_simple(stl, M, 2500))
    print(f'cascade control: {d} (campaign blockMesh/snappy/0/constant verbatim)')


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else FT
    names = [a for a in sys.argv[2:] if not a.startswith('-')]
    central = '--central' in sys.argv
    for name, stl, M in CASES:
        if names and name not in names:
            continue
        d = os.path.join(out, name)
        # numerics/thermo copied verbatim from the campaign case (same family)
        src = os.path.join(CAMPAIGN, name)
        for f in ('system/meshQualityDict', 'constant/thermophysicalProperties',
                  'constant/turbulenceProperties'):
            s = os.path.join(src, f)
            if not os.path.exists(s):
                s = os.path.join(CAMPAIGN, 'naca0012_M0.800_a1.25', f)
            os.makedirs(os.path.dirname(os.path.join(d, f)), exist_ok=True)
            shutil.copy2(s, os.path.join(d, f))
        os.makedirs(os.path.join(d, 'constant', 'triSurface'), exist_ok=True)
        shutil.copy2(os.path.join(STLDIR, stl),
                     os.path.join(d, 'constant', 'triSurface', stl))
        write(os.path.join(d, 'system', 'blockMeshDict'), blockmesh())
        write(os.path.join(d, 'system', 'snappyHexMeshDict'), snappy(stl))
        if central:
            for f in ('system/fvSchemes', 'system/fvSolution'):
                s = os.path.join(src, f)
                if not os.path.exists(s):
                    s = os.path.join(CAMPAIGN, 'naca0012_M0.800_a1.25', f)
                shutil.copy2(s, os.path.join(d, f))
            write(os.path.join(d, 'system', 'controlDict'), control(stl, M))
        else:
            write(os.path.join(d, 'system', 'fvSchemes'), fvschemes_simple())
            write(os.path.join(d, 'system', 'fvSolution'), fvsolution_simple())
            write(os.path.join(d, 'system', 'controlDict'), control_simple(stl, M))
        write(os.path.join(d, '0', 'p'), field_p(stl))
        write(os.path.join(d, '0', 'U'), field_u(stl, M))
        write(os.path.join(d, '0', 'T'), field_t(stl))
        print(f'{name:24s} M={M:.3f} magUInf={M * A_INF:.3f} '
              f'rho={RHO_INF:.5f} a={A_INF:.4f} y=+-{YBOX} solver='
              f'{"rhoCentralFoam" if central else "rhoSimpleFoam"}')
    build_control(out)


if __name__ == '__main__':
    main()
