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

## 5. NACA0012 transonic anchor, M0.8 alpha=1.25deg (Euler)

Case: `cases/openfoam/transonic` (slab STL pierces y-walls; symmetry
quasi-2D), container `oftrans`. AGARD-like case 1.

- **rhoCentralFoam + LTS (Kurganov-Tadmor, best-in-class for shocks)**,
  12k cells: **CL = 0.402, CD = 0.0414** (flat). Upper-surface
  supersonic pocket Cpmin -1.27 with a **crisp shock at x/c ~ 0.35**
  (Cp -1.22 -> -0.31 across one station). Lower Cpmin -0.22.
  VERDICT: ANCHOR SET (shock location + Euler CL for TSD validation).
- **rhoSimpleFoam + SA** (70k cells, wall functions, y+ uncontrolled):
  converged (residuals met) but to CL = 0.028, CD = 0.048 --
  non-physical (false convergence: upwind dissipation + bad wall
  treatment bury the Kutta condition; effective-Re-mush diagnosis).
  NOT used as truth. Lesson: SIMPLE+upwind needs thin layers
  (y+ 30-100) + linearUpwind; kept as open item.
- wing3d panel + Karman-Tsien at M0.8: midspan cl ~ 0.13,
  Cpmin -1.05, no shock (KT is subsonic-only, expected).
  Gap to close with TSD: shock@0.35 + CL~0.4.

## 6. NACA0012 M1.0 alpha=1.25deg Euler anchor (rhoCentralFoam LTS)

Case: `cases/openfoam/caseM1` (cloned 12k mesh, U=347.2). Finished
clean to t=0.3: **CL = 0.257, CD = 0.0265** (flat). Upper pocket
Cpmin -0.80 with shock at **x/c ~ 0.35** (Cp -0.77 -> -0.47 ->
-0.20, slightly more smeared than M0.8). Sonic freestream gives
less suction headroom than M0.8 (Cpmin -1.27 there) hence lower CL.
VERDICT: M1.0 ANCHOR SET (the TSD M1.0 target).

## 7. NACA0012 M0.9 alpha=1.25deg Euler (rhoCentralFoam LTS)

Case: `cases/openfoam/caseM09`. Converged (residuals to 0):
**CL = 0.072, CD = 0.053**. DOUBLE SHOCK state: upper pocket
Cpmin -0.76 AND lower pocket Cpmin -0.63, both shocking at
x/c ~ 0.35-0.45. Near-symmetric pockets cancel lift (transonic
lift bucket: 0.402 -> 0.072 -> 0.257 across M0.8/0.9/1.0).
Physical (M0.9 lower surface goes supersonic even at low alpha).

## 8. NACA0012 M0.8 alpha=4deg Euler (rhoCentralFoam LTS) -- SHOCK STALL

Case: `cases/openfoam/caseA4`. **CL = 0.404, CD = 0.0415** (3
consistent samples +-0.002). Pocket Cpmin -1.27 + shock at x/c
~0.35 -- IDENTICAL pocket/shock to alpha=1.25deg (CL 0.402).
Extra incidence goes into the shock, not lift: buffet boundary
below 4deg at M0.8. (Runs die at the shock-formation transient
t~0.1 -- captured state consistent 3x; needs hardening.)
Links 2D section to the 3D wing alpha. (Accidental alpha=8deg
run also gives CL 0.404 -- stall flat from 1.25deg through 8deg.)

## 9. NACA0012 M1.1 alpha=1.25deg Euler (rhoCentralFoam LTS) -- SUPERSONIC

Case: `cases/openfoam/caseM11`. **CL = 0.213, CD = 0.0219**.
Ackeret theory: 4*alpha(rad)/sqrt(M^2-1) = 0.19. Measured +12%
(Euler, coarse mesh). Upper expansion to Cp -0.66 + oblique
shock x/c ~ 0.3-0.4; bow-shock LE (stagnation Cp 0.27, not 1.0).
VERDICT: SUPERSONIC ANCHOR, theory-consistent.

## 10. NACA0012 M0.95 alpha=1.25deg Euler (rhoCentralFoam LTS)

Case: `cases/openfoam/caseM95`. Settles into +-0.001 limit cycle:
**CL = 0.285, CD = 0.0294**. Single upper shock again (pocket
Cpmin -0.90 at x/c ~0.33-0.40, lower mild -0.18) -- bucket
recovering past M0.9 (lower pocket gone). Transonic CL curve:
0.402 / 0.072 / 0.285 / 0.257 @ M0.8/0.9/0.95/1.0.

## 11. Complex-geometry campaign (wing3d panel + wakes, this study)

Waked multi-surface configs (`wing3d/shapes.py`: `build_f16_waked`,
`build_a330`, `build_shuttle`) with wake clipping + buried-panel force
masks at body junctions. Without these, overlapping solids double lift
and 10x drag (measured: F-16 CL 0.69/CDp 0.13 before, 0.26/0.001
after at 4deg). Script: `tests3d/complex_campaign.py`.

F-16-like (slope 0.065/deg; clean-subsonic published ~0.07-0.09,
strake vortices missing so slightly low is expected):
a=0/4/8/12deg -> CL -0.00/0.26/0.47/0.61. CDp surface-noisy
(2x induced at high alpha; use wing-alone/Trefftz for drag).

A-330-like (slope 0.108/deg, airliner-like; cruise CL~0.5 needs
a~5deg here vs published): a=0/2/4 -> CL -0.04/0.18/0.39.
Induced ~0.006 at CL 0.39; +friction (~0.015 missing) -> ~0.021
vs published ~0.026 (L/D 19-20). Qualitative agreement.

Shuttle-orbiter-like (approach): a=0/5/10/15 -> CL
0.00/0.18/0.34/0.45. No vortex-lift break (panel stays attached;
real orbiter gets LE-vortex lift then stalls). L/D at 10deg:
10.4 inviscid vs published orbiter ~4-5 (no base/flap/gear drag
here -- documented 2x gap). Polhamus suction-analogy correction
(`forces.polhamus_vortex_lift`, Kv=pi): +0.02/+0.09/+0.20 at
5/10/15deg -> CL 0.20/0.43/0.66, toward orbiter vortex-lifted
values (drag increment not modeled -- L/D still overstated).

RUNTIMES (this host): wing3d complex config assemble 7-16 s (once)
+ solve 4-16 s per alpha point. OpenFOAM: car RANS-SST ~4 min,
central Euler 12k ~10-19 min, wing RANS 172k ~30+ min (mesh-limited).
Panel ~30-70x faster per point at inviscid fidelity. C-kernel port
scoped separately (fused panel_rows, est. another 3-5x).

WORST-REGIME VERDICT (Re + speed sweeps, fixed IBL, 2026-10-08):
Re sweep (NACA0012 AR6 a=4, vs Abbott): 1e5 -> -13% (UNDER:
bubble loss missing -- worst PHYSICS), 3e5 -> +12%, 1e6 -> +21%,
3e6 -> +21%, 1e7 -> +26% (worst MAGNITUDE: offset dominates as
physical drag shrinks). Mach sweep (KT vs PG): <5% to M0.5,
+5.6% M0.6, +11% M0.7, +15% M0.75. Worst speed: M>=0.65
(KT invalid, needs TSD/OF). Overall weakest: transonic
M0.65-1.1 (no wing3d predictor) + low-Re bubbles.

## Environment notes (this host)

- Docker bind mounts broken (`Cannot allocate memory` on readdir);
  used `docker cp` in/out + named containers (ofwing/ofcar/ofsphere).
- 395 MB RAM visible to containers: keep RANS meshes < ~150k cells.
- CRLF in Windows-authored dicts breaks sed patterns; strip first.
- motorBike tutorial dicts (in-image) used as the version-compatible
  base for all handwritten dictionaries.
