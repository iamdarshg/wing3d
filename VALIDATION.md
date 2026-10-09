# wing3d validation report

Date: 2026-10-06. All 3D, no strips. Suite: `tests3d/test_shapes.py` (9/9 green).

## 1. Public-data validation (no OpenFOAM needed)

| Case | wing3d | Reference | Verdict |
|---|---|---|---|
| Sphere Cp | rms 0.08 (12x24) | analytic 1-2.25sin^2 | GOOD, converges |
| Sphere drag | -0.0 | d'Alembert 0 | EXACT |
| Cylinder Cp | peak -2.62 | analytic -2.99 | FAIR (12% low) |
| Cube | stag ~1, CD~0, symmetric | d'Alembert/physics | GOOD |
| Ellipsoid | converges, CD~0 | physics | GOOD |
| NACA0012 2D cl(4deg) | 0.46 (own Hess-Smith) | Abbott/von Doenhoff 0.44 | GOOD |
| NACA0012 AR6 slope | 0.09/deg | tunnel 0.08-0.09, VLM 0.098, LL 0.082 | GOOD |
| NACA0012 AR6 CL(4deg) | 0.365 invisc / 0.364 coupled | LL 0.33, VLM 0.39 | GOOD (in range) |

Notes: finite TE (0.2% gap, physical) regularizes the sharp-TE corner
singularity; sharp TE overpredicts ~30%. Viscous coupling (streamline
IBL + transpiration) converges stably: -0.8% decambering, +drag, no
separation at Re 3e6/alpha 4deg. CD runs ~40% high (IBL friction +
residual TE drag; calibration item).

## 2. OpenFOAM RANS validation (SST, Docker, this study)

| Case | OpenFOAM RANS | Reference | Verdict |
|---|---|---|---|
| Ahmed-like car | CD 0.275 (conv. flat) | published Ahmed 0.26-0.28 | VALIDATED |
| Sphere Re=100 lam. | CD 1.106 (conv. flat) | Schiller-Naumann 1.09 | VALIDATED (+1.5%) |
| NACA0012 AR6 wing | CL 0.11, CD 0.022 | expected ~0.35/~0.018 | MESH-LIMITED (see below) |

Wing RANS details: 152-172k cells (level 3-4 + LE box level 6),
converged flat. Leading-edge cells (~0.01) too coarse for the nose
radius + faceted STL sute; flow misses suction peak (likely LE
separation). Needs >300k cells, exceeding the 395 MB available to
Docker here (OOM at 318k). 2D airfoil tutorial also unusable as-is
(unknown mesh orientation + Aref ambiguity).

## 3. Cross-method wing summary (NACA0012 AR6, alpha 4deg, Re 3e6)

- Lifting line: 0.33. VLM (in-repo): 0.39. wing3d inviscid: 0.365.
  wing3d coupled: 0.364. OpenFOAM RANS: 0.11 stalled (mesh-limited).
- Inviscid methods agree within ~10%; RANS needs finer mesh.

## 4. Transonic Euler anchors (rhoCentralFoam LTS, quasi-2D slabs, this study)

NACA0012 section, 12k cells each. Full transcripts in
`cases/openfoam/TRANSCRIPTS.md`, machine-readable data + plots in
the UI (`/validation`) and `cases/openfoam/transonic/summary.json`.

| M | alpha | CL | CD | Shock | Note |
|---|---|---|---|---|---|
| 0.8 | 0 deg | 0.0003 | 0.047 | symmetric | symmetry check; ~0.04 spurious floor |
| 0.8 | 1.25 deg | 0.069 | 0.049 | x/c 0.35 | TRUE rerun (LL 0.074, -7%) |
| 0.8 | 4 deg | 0.219 | 0.046 | x/c 0.3-0.43, Cpmin -0.96 | TRUE rerun, verified inflow |
| 0.8 | 8 deg | 0.404 | 0.0415 | x/c 0.35, Cpmin -1.27 | 5 consistent runs |
| 0.885 | 1.25 deg | 0.070 | 0.051 | double shock | gas-corrected, LL -13% |
| 0.95 | 1.25 deg | 0.083 | 0.065 | double shock@0.6 | TRUE rerun (0.285 was transient) |
| 1.0 | 1.25 deg | 0.079 | 0.077 | upper x/c 0.8, Cpmin -0.80 | M1.0 TARGET anchor |
| 1.1 | 1.25 deg | 0.053 | 0.103 | bow + oblique | supersonic, N2-verified |
| 0.75 | 2 deg (Harris tunnel) | 0.377 | — | x/c 0.5-0.55 | transonic truth (no mesh issues) |
| (All OpenFOAM rows: finite-wing AR1.6 slabs, N2-consistent, LL-validated) | | | | | |

Caveats: Euler (no boundary layers); coarse mesh (no grid-convergence
bars yet; spurious drag floor ~0.04 at M0.8); rhoSimpleFoam+upwind
false-converges on these cases (rejected as truth, documented);
runs die at shock-formation transients (chunked restarts used).

## 5. How to reproduce

- Panel: `python tests3d/validate_all.py`, `pytest tests3d/test_shapes.py`
- UI: `python -m wing3d.app` (localhost:8080, STL upload)
- OpenFOAM: `wing3d/openfoam.py` writes cases; images/containers as in
  this report (`ofwing` wing 172k, `ofcar` car 82k). Pull
  `opencfd/openfoam-default:2406`. Note: Docker bind mounts broken on
  this host; use `docker cp` in/out (documented in session).
- Plots: `cases/openfoam/polar_compare.png`, `sphere_compare.png`.

## 6. Known limitations / next steps

1. Wing RANS fidelity needs more RAM (finer LE mesh + layers).
2. Bluff bodies (car/F-16): panel is qualitative only (no separation);
   RANS is the reference (car validated).
3. CD ~40% high in coupled solver (IBL friction calibration; Michel
   transition guard + Thwaites lambda clamp + bubble model landed;
   bubble-loss (Horton) still missing at low Re).
4. High-alpha (>6deg) panel degrades (TE singularity grows).
5. TSD nonlinear unclosed: linear+TSD reproduces panel (1.4%), but
   shock capture needs a free (Kutta-updated) bound vortex --
   currently prescribed (pins linear state). OF anchors stand in.
6. Transonic anchors need grid-convergence bars (finer LE mesh).
7. Perf: C kernel LANDED (cffi fused panel_rows, validated 1e-11,
   2.5-3.4x end-to-end, solve 8.7x; WING3D_NO_C=1 fallback).

## 7. Bias directions (truth vs lies)

Mean absolute error ~8-10%. Direction matters more than magnitude:

| Check | Loss | Direction |
|---|---|---|
| Sphere Cp/Cd, car CD, 2D cl, wing CL | 1-5% | ~neutral |
| Cylinder suction peak | -12% | conservative (underpredicts) |
| Re drag, high Re | +12..+29% | CONSERVATIVE (overpredicts drag) |
| Re drag, low Re (bubbles) | -9..-13% | FLATTERING (lies: drag lower than reality) |
| KT transonic | +70..+289%, then collapse | broken (delusional then garbage) |
| Shock-fitted M0.65-0.9 | +-5..11% | mildly conservative |
| TSD/fullpot, OF anchors vs LL | -5..-15% | conservative (mesh starvation) |
| Complex CL slopes | -10..-20% | conservative (missing vortex/separation lift) |

Errors lean CONSERVATIVE (safe direction) except low-Re drag
(optimistic lie -- Horton bubble loss modeled, burst unvalidated)
and KT (rejected as truth above M0.65).

## 8. Complex geometries with known numbers

Geometries used: NACA0012 wings/slabs (AR6/AR16/AR12-proxy), sphere,
cube, cylinder, ellipsoid, Ahmed car, F-16-like, A330-like,
shuttle-like, F-5E-like (new). Overlapping solids need wake
clipping + buried-panel masks (else 2x lift/10x drag); sharp
strakes (LEX) break low-order panel globally -- modeled via
Polhamus instead (wedge outward-normal bug also fixed).

| Config | wing3d | Published | Verdict |
|---|---|---|---|
| F-5E-like slope | 0.06/deg (+Polhamus vortex at high α) | no hard public polar found (only blowing/spin studies) | qualitative; LEX analytic |
| Orbiter approach | L/D 10.4 inviscid | ~3 (19 deg glideslope, sourced) to ~4.5 max (cited) | 2-3x over (no base/flap/gear/trim drag) |
| A330-like | slope 0.108/deg, L/D ~18 est | L/D 19-20 claims | qualitative agreement |
| F-16-like | slope 0.065/deg | ~0.07-0.09 est | plausible, strakes missing |
| A350 / Bombardier Global | not built | estimates only, no hard data | skipped (low value vs A330) |
