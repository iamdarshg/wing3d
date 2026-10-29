#!/bin/bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc >/dev/null 2>&1
cd /work/tune1 || exit 1
# robust compressible settings for a ~89 deg non-orthogonal BL mesh
cat > system/fvSchemes <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object fvSchemes; }
ddtSchemes { default steadyState; }
gradSchemes { default cellLimited Gauss linear 0.5; grad(p) cellLimited Gauss linear 0.5; }
divSchemes {
 default none;
 div(phi,U) bounded Gauss upwind;
 div(phi,K) bounded Gauss upwind;
 div(phi,Ekp) bounded Gauss upwind;
 div(phi,e) bounded Gauss upwind;
 div(phi,T) bounded Gauss upwind;
 div(phi,nuTilda) bounded Gauss upwind;
 div((nuEff*dev2(T(grad(U))))) Gauss linear;
 div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited 0.33; }
interpolationSchemes { default linear; }
snGradSchemes { default limited 0.33; }
fluxRequired { default no; p; }
wallDist { method meshWave; }
EOF
cat > system/fvSolution <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }
solvers
{
 p { solver GAMG; tolerance 1e-7; relTol 0.1; smoother DICGaussSeidel; nCellsInCoarsestLevel 500; }
 U { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-8; relTol 0.1; }
 e { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-8; relTol 0.1; }
 T { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-8; relTol 0.1; }
 nuTilda { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-8; relTol 0.1; }
}
rhoSimple { momentumPredictor no; nNonOrthogonalCorrectors 5; consistent no; rhoMin 0.05; rhoMax 8; }
relaxationFactors {
 fields { p 0.1; rho 0.2; }
 equations { U 0.3; e 0.5; nuTilda 0.5; T 0.5; }
}
EOF
cp -r /work/naca_su2_case3/0 . 2>/dev/null
cp /work/naca_su2_case3/constant/thermophysicalProperties \
   /work/naca_su2_case3/constant/turbulenceProperties constant/ 2>/dev/null
sed -i s/WING3D/wing/g constant/polyMesh/boundary
sed -i s/SYMMETRY/symmetry/g constant/polyMesh/boundary
sed -i s/FARFIELD/farfield/g constant/polyMesh/boundary
sed -i 's/^endTime .*/endTime         300;/' system/controlDict
sed -i 's/^writeInterval .*/writeInterval   100;/' system/controlDict
timeout 600 rhoSimpleFoam > log.run 2>&1
echo "rc=$?"
grep -cE '^Time = ' log.run
grep -E 'Cl:' log.run | tail -2
grep -A4 'FOAM FATAL' log.run | head -6
