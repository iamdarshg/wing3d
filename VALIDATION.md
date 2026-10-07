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

## 4. How to reproduce

- Panel: `python tests3d/validate_all.py`, `pytest tests3d/test_shapes.py`
- UI: `python -m wing3d.app` (localhost:8080, STL upload)
- OpenFOAM: `wing3d/openfoam.py` writes cases; images/containers as in
  this report (`ofwing` wing 172k, `ofcar` car 82k). Pull
  `opencfd/openfoam-default:2406`. Note: Docker bind mounts broken on
  this host; use `docker cp` in/out (documented in session).
- Plots: `cases/openfoam/polar_compare.png`, `sphere_compare.png`.

## 5. Known limitations / next steps

1. Wing RANS fidelity needs more RAM (finer LE mesh + layers).
2. Bluff bodies (car/F-16): panel is qualitative only (no separation);
   RANS is the reference (car validated).
3. CD ~40% high in coupled solver (IBL friction calibration).
4. High-alpha (>6deg) panel degrades (TE singularity grows).
