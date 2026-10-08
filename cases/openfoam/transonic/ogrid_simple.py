"""Write SIMPLE-transonic (laminar slip) dicts for the O-grid case."""
import os

D = 'D:/CodeProjects/cfd/cases/openfoam/ogrid_simple'
os.makedirs(D + '/system', exist_ok=True)
Umag, rho = 282.446, 1.1231
W = lambda s: open(D + '/system/' + s.split('/')[0], 'w', newline='\n')
open(D + '/system/controlDict', 'w', newline='\n').write(
    'FoamFile { version 2.0; format ascii; class dictionary; '
    'object controlDict; }\napplication rhoSimpleFoam;\n'
    'startFrom startTime; startTime 0; stopAt endTime; endTime 2000;\n'
    'deltaT 1; writeControl timeStep; writeInterval 2000;\n'
    'purgeWrite 0; writeFormat ascii; writePrecision 6; '
    'writeCompression off;\ntimeFormat general; timePrecision 6; '
    'runTimeModifiable true;\nfunctions {\n forces {\n'
    '  type forceCoeffs; libs ("libforces.so");\n'
    '  writeControl timeStep; writeInterval 50;\n'
    '  patches ("airfoil");\n'
    '  rho rhoInf; rhoInf %.4f; CofR (0.25 0 0);\n'
    '  liftDir (0 0 1); dragDir (1 0 0); pitchAxis (0 1 0);\n'
    '  magUInf %.3f; lRef 1.0; Aref 0.1;\n }\n'
    ' residuals { type residuals; libs ("libutilityFunctionObjects.so");\n'
    '  writeControl timeStep; writeInterval 50; '
    'fields (p U h); }\n'
    ' limT {\n  type limitTemperature; '
    'libs ("libfieldFunctionObjects.so");\n'
    '  writeControl timeStep; writeInterval 1;\n  min 200; max 500;\n }\n'
    ' surfCp {\n  type surfaces; libs ("libsampling.so");\n'
    '  writeControl writeTime; surfaceFormat raw; fields (p);\n'
    '  interpolationScheme cell;\n'
    '  surfaces ( foil { type patch; patches ("airfoil"); } );\n }\n}\n'
    % (rho, Umag))
open(D + '/system/fvSchemes', 'w', newline='\n').write(
    'FoamFile { version 2.0; format ascii; class dictionary; '
    'object fvSchemes; }\nddtSchemes { default steadyState; }\n'
    'gradSchemes { default Gauss linear; }\n'
    'divSchemes { default none; div(phi,U) Gauss upwind; '
    'div(phi,h) Gauss upwind; div(phi,e) Gauss upwind;\n'
    ' div(phid,p) Gauss upwind; div(phi,Ekp) Gauss upwind;\n'
    ' div(phi,K) Gauss upwind;\n'
    ' div((muEff*dev2(T(grad(U))))) Gauss linear;\n'
    ' div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear; }\n'
    'laplacianSchemes { default Gauss linear corrected; }\n'
    'interpolationSchemes { default linear; }\n'
    'snGradSchemes { default corrected; }\n')
open(D + '/system/fvSolution', 'w', newline='\n').write(
    'FoamFile { version 2.0; format ascii; class dictionary; '
    'object fvSolution; }\nsolvers {\n'
    ' p { solver GAMG; tolerance 1e-7; relTol 0.01; '
    'smoother GaussSeidel; }\n'
    ' "(U|h|e)" { solver smoothSolver; smoother symGaussSeidel; '
    'tolerance 1e-6; relTol 0.1; }\n}\n'
    'SIMPLE { nNonOrthogonalCorrectors 3; transonic yes;\n'
    ' residualControl { p 1e-4; U 1e-4; h 1e-4; } }\n'
    'relaxationFactors { fields { p 0.3; rho 0.05; } '
    'equations { U 0.5; h 0.5; } }\n')
print('wrote SIMPLE O-grid dicts')
