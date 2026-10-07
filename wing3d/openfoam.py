"""OpenFOAM case generation + result parsing for wing3d validation.

Generates simpleFoam + SpalartAllmaras cases around wing3d-exported STLs
(blockMesh background + snappyHexMesh). Run inside the container:
  docker run --rm -v <case>:/case -w /case opencfd/openfoam-default:2406 \
    bash -c "source .../bashrc && ./Allrun"
"""
import os
import numpy as np


def write_case(root, stl_name='wing.stl', chord=1.0, span=6.0, alpha_deg=4.0,
               Umag=30.0, nu=1e-5, level=(3, 4), n_layers=3):
    os.makedirs(os.path.join(root, 'system'), exist_ok=True)
    os.makedirs(os.path.join(root, '0'), exist_ok=True)
    os.makedirs(os.path.join(root, 'constant', 'triSurface'), exist_ok=True)
    a = np.radians(alpha_deg)
    ux, uz = Umag * np.cos(a), Umag * np.sin(a)
    # domain (tunnel)
    x0, x1 = -5 * chord, 10 * chord
    y0, y1 = -span / 2 - 4 * chord, span / 2 + 4 * chord
    z0, z1 = -5 * chord, 5 * chord
    nx, ny, nz = 40, 30, 30
    with open(os.path.join(root, 'system', 'blockMeshDict'), 'w') as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; object blockMeshDict; }}
convertToMeters 1;
vertices (
 ({x0} {y0} {z0}) ({x1} {y0} {z0}) ({x1} {y1} {z0}) ({x0} {y1} {z0})
 ({x0} {y0} {z1}) ({x1} {y0} {z1}) ({x1} {y1} {z1}) ({x0} {y1} {z1}));
blocks (hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1));
edges ();
boundary (
 inlet {{ type patch; faces ((0 4 7 3)); }}
 outlet {{ type patch; faces ((1 2 6 5)); }}
 walls {{ type wall; faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }});
""")
    with open(os.path.join(root, 'system', 'snappyHexMeshDict'), 'w') as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; object snappyHexMeshDict; }}
castellatedMesh true; snap true; addLayers true;
geometry {{ {stl_name} {{ type triSurfaceMesh; file "{stl_name}"; }} }}
castellatedMeshControls {{
 maxLocalCells 500000; maxGlobalCells 800000; minRefinementCells 10;
 maxLoadUnbalance 0.1; nCellsBetweenLevels 3;
 features ();
 refinementSurfaces {{ {stl_name} {{ level ({level[0]} {level[1]}); patchInfo {{ type wall; }} }} }}
 resolveFeatureAngle 25;
 refinementRegions {{ {stl_name} {{ mode distance; levels ((0.5 3) (1.0 2)); }} }}
 locationInMesh ({(x0 + x1) / 2} {(y0 + y1) / 2 + span} {(z0 + z1) / 2});
 allowFreeStandingZoneFaces true;
}}
snapControls {{ nSmoothPatch 3; tolerance 2.0; nSolveIter 30; nRelaxIter 5; }}
addLayersControls {{
 relativeSizes true; layers {{ "{stl_name}" {{ nSurfaceLayers {n_layers}; }} }}
 expansionRatio 1.2; finalLayerThickness 0.5; minThickness 0.1;
 nGrow 0; featureAngle 60; slipFeatureAngle 30;
 nRelaxIter 3; nSmoothSurfaceNormals 1; nSmoothNormals 3;
 nSmoothThickness 10; maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3;
 minMedianAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 50;
}}
meshQualityControls {{
 #include "meshQualityDict"
 nSmoothScale 4; errorReduction 0.75;
}}
writeFlags (scalarLevels layerSets);
mergeTolerance 1e-6;
""")
    with open(os.path.join(root, 'system', 'controlDict'), 'w') as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; object controlDict; }}
application simpleFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 800;
deltaT 1; writeControl timeStep; writeInterval 800;
purgeWrite 0; writeFormat ascii; writePrecision 6; writeCompression off;
timeFormat general; timePrecision 6; runTimeModifiable true;
functions {{
 forces {{
  type forceCoeffs; libs ("libforces.so");
  writeControl timeStep; writeInterval 50;
  patches ("{stl_name.replace('.stl', '')}");
  rho rhoInf; rhoInf 1.0; CofR (0.25 0 0);
  liftDir (0 0 1); dragDir (1 0 0); pitchAxis (0 1 0);
  magUInf {Umag}; lRef {chord}; Aref {chord * span};
 }}
 residuals {{ type residuals; libs ("libutilityFunctionObjects.so");
  writeControl timeStep; writeInterval 50; fields (p U k omega nuTilda); }}
}}
""")
    with open(os.path.join(root, 'system', 'fvSchemes'), 'w') as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; object fvSchemes; }
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; grad(p) Gauss linear; grad(U) Gauss linear; }
divSchemes { default none; div(phi,U) Gauss linear; div(phi,nuTilda) Gauss linear; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    with open(os.path.join(root, 'system', 'fvSolution'), 'w') as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }
solvers {
 p { solver GAMG; tolerance 1e-6; relTol 0.05; smoother GaussSeidel; }
 "(U|nuTilda)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-5; relTol 0.1; }
}
SIMPLE { nNonOrthogonalCorrectors 0; residualControl { p 1e-4; U 1e-4; "(nuTilda)" 1e-4; } }
relaxationFactors { equations { U 0.7; nuTilda 0.7; } }
""")
    with open(os.path.join(root, 'constant', 'transportProperties'), 'w') as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; object transportProperties; }}
transportModel Newtonian; nu [0 2 -1 0 0 0 0] {nu};
""")
    with open(os.path.join(root, 'constant', 'turbulenceProperties'), 'w') as f:
        f.write("""FoamFile { version 2.0; format ascii; class dictionary; object turbulenceProperties; }
simulationType RAS;
RAS { RASModel SpalartAllmaras; turbulence on; printCoeffs on; }
""")
    with open(os.path.join(root, '0', 'U'), 'w') as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class internalField; object U; }}
dimensions [0 1 -1 0 0 0 0];
internalField uniform ({ux} 0 {uz});
boundaryField {{
 inlet {{ type fixedValue; value uniform ({ux} 0 {uz}); }}
 outlet {{ type inletOutlet; inletValue uniform ({ux} 0 {uz}); value uniform ({ux} 0 {uz}); }}
 walls {{ type slip; }}
 "{stl_name.replace('.stl', '')}" {{ type noSlip; }}
}}
""")
    with open(os.path.join(root, '0', 'p'), 'w') as f:
        f.write("""FoamFile { version 2.0; format ascii; class internalField; object p; }
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField {
 inlet { type zeroGradient; }
 outlet { type fixedValue; value uniform 0; }
 walls { type zeroGradient; }
 "%s" { type zeroGradient; }
}
""" % stl_name.replace('.stl', ''))
    nut = 3 * nu
    with open(os.path.join(root, '0', 'nut'), 'w') as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class internalField; object nut; }}
dimensions [0 2 -1 0 0 0 0];
internalField uniform {nut};
boundaryField {{
 inlet {{ type calculated; value uniform {nut}; }}
 outlet {{ type calculated; value uniform {nut}; }}
 walls {{ type nutUSpaldingWallFunction; value uniform {nut}; }}
 "{stl_name.replace('.stl', '')}" {{ type nutUSpaldingWallFunction; value uniform {nut}; }}
}}
""")
    with open(os.path.join(root, '0', 'nuTilda'), 'w') as f:
        f.write(f"""FoamFile {{ version 2.0; format ascii; class internalField; object nuTilda; }}
dimensions [0 2 -1 0 0 0 0];
internalField uniform {3 * nu};
boundaryField {{
 inlet {{ type fixedValue; value uniform {3 * nu}; }}
 outlet {{ type zeroGradient; }}
 walls {{ type zeroGradient; }}
 "{stl_name.replace('.stl', '')}" {{ type fixedValue; value uniform 0; }}
}}
""")
    with open(os.path.join(root, 'Allrun'), 'w') as f:
        f.write("""#!/bin/bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc
set -e
blockMesh > log.blockMesh 2>&1
snappyHexMesh -overwrite > log.snappy 2>&1
potentialFoam -writep > log.potential 2>&1 || true
simpleFoam > log.simpleFoam 2>&1
""")
    print('wrote case', root)


def parse_force_coeffs(root):
    """Parse latest forceCoeffs output -> dict(CL, CD)."""
    import glob
    import re
    files = sorted(glob.glob(os.path.join(root, 'postProcessing', 'forces',
                                           '*', 'forceCoeffs.dat')))
    if not files:
        return None
    qs = []
    with open(files[-1]) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            try:
                vals = [float(x.strip('()')) for x in parts[1:4]]
            except Exception:
                continue
            qs.append(vals)
    if not qs:
        return None
    q = np.array(qs[-50:])
    # forceCoeffs.dat columns: time Cm Cd Cl (moment, drag, lift)
    return {'CL': float(q[:, 2].mean()), 'CD': float(q[:, 1].mean()),
            'CL_std': float(q[:, 2].std()), 'CD_std': float(q[:, 1].std())}


def parse_log_forces(logpath, last=100):
    """Parse 'forceCoeffs forces write:' blocks from a solver log.

    Returns dict with history lists + final means.
    """
    import re
    hist = {'CD': [], 'CDp': [], 'CDv': [], 'CL': [], 'CLp': [], 'CLv': []}
    cur = {}
    with open(logpath, errors='ignore') as fh:
        for line in fh:
            m = re.match(r'\s*Cd:\s*(\S+)\s+(\S+)\s+(\S+)', line)
            if m:
                cur['CD'], cur['CDp'], cur['CDv'] = map(float, m.groups())
                continue
            m = re.match(r'\s*Cl:\s*(\S+)\s+(\S+)\s+(\S+)', line)
            if m:
                cur['CL'], cur['CLp'], cur['CLv'] = map(float, m.groups())
                if len(cur) == 6:
                    for k in hist:
                        hist[k].append(cur[k])
                    cur = {}
    out = {}
    for k, v in hist.items():
        a = np.array(v[-last:])
        out[k] = float(a.mean()) if len(a) else float('nan')
        out[k + '_std'] = float(a.std()) if len(a) else float('nan')
        out[k + '_hist'] = v
    out['n'] = len(hist['CL'])
    return out
