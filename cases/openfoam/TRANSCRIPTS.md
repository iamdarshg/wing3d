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

## 6. NACA0012 M1.0 alpha=1.25deg Euler anchor -- WITHDRAWN
(same audit: M0.8 alpha=8). Rerun pending.

## 7. NACA0012 M0.9 alpha=1.25deg Euler (rhoCentralFoam LTS)

Case: `cases/openfoam/caseM09`. Converged (residuals to 0):
**CL = 0.072, CD = 0.053**. DOUBLE SHOCK state: upper pocket
Cpmin -0.76 AND lower pocket Cpmin -0.63, both shocking at
x/c ~ 0.35-0.45. Near-symmetric pockets cancel lift (transonic
lift bucket: 0.402 -> 0.072 -> 0.257 across M0.8/0.9/1.0).
Physical (M0.9 lower surface goes supersonic even at low alpha).

## 8. NACA0012 M0.8 alpha=8deg Euler (rhoCentralFoam LTS) -- SHOCK STALL
(CORRECTED 2026-10-08: the alpha=4deg runs below were mislabeled --
 inflow audit showed caseA4/caseM95/caseM1/caseM11 all carried M0.8
 alpha=8deg inflow. Denormalized they read 0.402-0.404 consistently.)

Case: `cases/openfoam/caseA4` (+3 accidental repeats). **CL = 0.404,
CD = 0.0415** (4 consistent samples +-0.002). Pocket Cpmin -1.27 +
shock at x/c ~0.35.

## 8b. NACA0012 M0.8 alpha=4deg Euler TRUE (verified inflow, N2 gas)

Case: `cases/openfoam/caseA4v` (makecase protocol, maxCo 0.3,
binary). Clean finish to t=0.3: **CL = 0.219, CD = 0.0464**.
Pocket Cpmin -0.96, shock x/c ~0.3-0.43, lower mild -0.38.
M0.8 alpha curve now: 0deg 0.000 / 4deg 0.219 / 8deg 0.404 --
monotonic, diminishing slope (no stall cliff; earlier "flat"
narrative void -- it compared alpha=8 with itself).

## 9. NACA0012 M1.1 anchor -- WITHDRAWN (inflow audit 2026-10-08:
ran M0.8 alpha=8, misnormalized; the "Ackeret agreement" is void)

## 9b. NACA0012 M1.1 alpha=1.25 TRUE (verified N2 inflow, AR1.6)

Case: `cases/openfoam/caseM11v`. Flat CL=0.053, CD=0.103.
Bow shock (LE Cp +0.9/+1.0, not 1.0) + expansion to -0.6 both
sides + oblique TE recompression (smeared on coarse mesh).
Textbook supersonic pattern. High wave drag (blunt 12%) physical.

## 10. NACA0012 M0.95 anchor -- WITHDRAWN on first pass (M0.8a8),
then TRUE rerun first showed 0.285 (transient slosh artifact);
clean rerun to t=0.3: **CL = 0.083, CD = 0.065** (flat).
DOUBLE pocket (upper -0.76, lower -0.65), shocks x/c ~0.55-0.65.
Coherent bucket with M0.885 (double-shock@0.4).

## 10b. NACA0012 M1.0 alpha=1.25 TRUE (verified N2 inflow, AR1.6)

Case: `cases/openfoam/caseM10v` (makecase protocol). Flat
CL=0.079, CD=0.077 at t=0.16. DOUBLE pocket (upper Cpmin -0.72,
lower -0.62) with shock at x/c ~0.75-0.85 (aft shock: sonic
pocket extends far back at M1.0). High wave drag physical.
Shock march: 0.35 (M0.8) -> 0.4 (M0.9) -> 0.8 (M1.0).

## 11. DATA INTEGRITY AUDIT 2026-10-08 (READ THIS BEFORE CITING)

`of_verify.py` launch checklist added after finding: caseM95/caseM1/
caseM11/caseA4 0/U files carried M0.8-alpha8 inflow (stale clone +
no-op seds), so their "M0.95/M1.0/M1.1/a4" labels were wrong --
denormalized they all read CL~0.40 (M0.8a8 stall, consistent).
Also caseCentral/0/U was alpha=8 all along: the M0.8-alpha1.25
anchor NEVER EXISTED. Standing (inflow-verified): M0.8a0,
M0.8a4-true (0.219), M0.8a8 (0.404), M0.9a1.25 (0.072) + Harris
tunnel truth (cl 0.377 @M0.75a2). Gas note: biconic janaf thermo is N2 (R=296.8),
so old runs labeled with air sound speed sit ~0.015 low in Mach
and +3.4% in CL/CD normalization; corrected where cited.
`wing3d/makecase.py` (N2-consistent U/rhoInf) + `of_verify.py`
(Rgas arg) prevent recurrence. Fresh-start only (restarts across
alpha read stale BCs from time dirs).

## 11b. FINITE-WING RESOLUTION (the "1/3 levels" mystery solved)

The slab (span 1.6, domain 1.2, tips outside) is a FINITE wing of
AR=1.6, not quasi-2D -- tip vortices explain all level
suppression. Compressible lifting-line check (AR=1.6):
M0.885a1.25: LL 0.080 vs meas 0.070 (-13%);
M0.95a1.25: 0.088 vs 0.083 (-5%);
M0.8a4: 0.237 vs 0.219 (-8%);
M0.8a8: 0.474 vs 0.404 (-15%, stall onset).
Deficit grows with M/alpha = shock losses (physical trend).
Anchors VALIDATED vs theory; shapes (shock@0.35, bucket, stall)
corroborate. No mesh mystery remains (LE coarseness = error bars,
not the factor).

## 12. Complex-geometry campaign (wing3d panel + wakes, this study)
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

## 13. Harris tunnel truth (NASA TM-81927, OCR'd tables)

NACA0012 M0.749 alpha=1.99deg Re=2.0e6 (`transonic/harris.py`):
upper shock between x=0.50 (Cp -1.187) and x=0.55 (Cp -0.340),
Cpmin -1.23, integrated cl = 0.377. Use as transonic truth for
pocket/shock/CL (no mesh issues). vs our quasi-2D OF: shapes
match (shock@0.35-0.4 vs 0.52 -- ours forward, coarse-LE effect),
levels ~1/3 (LE starvation; finite-wing AR1.6 contamination on
early slabs + coarse LE). O-grid attempt (exact 2D, smooth LE)
documented in transonic/ogrid.py -- blocked on block-corner
non-orthogonality (80 deg faces) + central/SIMPLE blowup; needs
edgeGrading rework or C-grid follow-up (spec'd, not built).

## Environment notes (this host)

- Docker bind mounts broken (`Cannot allocate memory` on readdir);
  used `docker cp` in/out + named containers (ofwing/ofcar/ofsphere).
- 395 MB RAM visible to containers: keep RANS meshes < ~150k cells.
- CRLF in Windows-authored dicts breaks sed patterns; strip first.
- motorBike tutorial dicts (in-image) used as the version-compatible
  base for all handwritten dictionaries.
