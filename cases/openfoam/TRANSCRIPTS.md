# OpenFOAM validation transcripts (real runs, Docker)

Image: `opencfd/openfoam-default:2406`. All cases: simpleFoam (+SST unless
noted), forceCoeffs from solver logs. Full logs kept alongside each case.

## 1. Ahmed-like car, slant 25deg, Re ~3.1M (SST)

Case: `cases/openfoam/car`, log: `log.simpleFoam` (800 iters, converged flat).

```
Cd:  0.275103  0.253923  0.0211804  0
Cl:  0.0530313 0.0529288 0.000102576 0
```
Columns: Total, Pressure, Viscous. CD = 0.275 (pressure 0.254 + viscous
0.021). Published Ahmed-body 25deg slant: CD 0.26-0.28. VERDICT: VALIDATED.

## 2. Sphere, Re = 100, laminar

Case: `cases/openfoam/sphere`, log: `log.simpleFoam2` (800 iters).
Convergence: CD 3.45 -> 2.06 -> ... -> 1.106 (flat).

```
Cd:  1.10588  0.525472  0.580409  0
Cl: -0.00201957 (symmetric, ~0 as expected)
```
Schiller-Naumann at Re=100: CD ~1.09. Measured 1.106 (+1.5%).
VERDICT: VALIDATED.

## 3. NACA0012 AR6 wing, alpha = 4deg, Re 3e6 (SST)

Case: `cases/openfoam/wing`. Status: MESH-LIMITED.
- Coarse mesh (172k): CL froze at 0.10-0.11 (leading-edge suction
  unresolved; likely LE separation on coarse cells).
- LE-refined mesh (318k) exceeds the 395 MB available to Docker here
  (OOM). 152k sharp-TE mesh also stalls at CL ~0.11.
- alpha = 0deg (same mesh): CD 0.024, CL 0.043 (mesh asymmetry noise).

Expected (LL/VLM/panel): CL ~0.33-0.39. The RANS needs finer LE mesh
than this host's Docker RAM allows. Panel/VLM/LL mutual agreement
(CL 0.33-0.39) stands as the wing validation; RANS wing is an open
item for bigger hardware, NOT presented as truth.

## 4. NACA0012 2D tutorial mesh

The bundled airFoil2D tutorial was investigated and dropped: unknown
mesh provenance/orientation and ambiguous force normalization made it
unsuitable as a reference. 2D validation rests on Abbott/von Doenhoff
via the in-repo Hess-Smith panel (cl 0.46 vs 0.44 at 4deg).

## Environment notes (this host)

- Docker bind mounts broken (`Cannot allocate memory` on readdir);
  used `docker cp` in/out + named containers (ofwing/ofcar/ofsphere).
- 395 MB RAM visible to containers: keep RANS meshes < ~150k cells.
- CRLF in Windows-authored dicts breaks sed patterns; strip first.
- motorBike tutorial dicts (in-image) used as the version-compatible
  base for all handwritten dictionaries.
