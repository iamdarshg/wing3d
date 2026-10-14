# M0.85-M1.1 High-Fidelity Campaign (2026-10-13)

## What was built
- **STL suite (100k+ verts each)**: naca0012 (flat plate of 12%), rae2822
  (supercritical, 3.1% camber), biconvex (8% symmetric arc), dfvlr_r4
  (~7.5% research). Each extruded slab, AR1.6, ~154k verts / 307k tris.
- **44 OpenFOAM Euler cases**: 4 geometries x 11 Mach points (0.85..1.1),
  rhoCentralFoam LTS, 86k-cell snappy meshes re-meshed per geometry.
- **4 TSD full-potential runs**: matched span/alpha, real per-geometry
  thickness profiles.

## Converged results (t=0.25, a=1.25deg, AR1.6)
| geo | M0.85 | M0.9 | M0.925 peak | M0.95 | M1.0 | M1.05 | M1.1 |
|---|---|---|---|---|---|---|---|
| naca0012 | 0.0728 | 0.0774 | 0.0799 | 0.0800 | 0.0628 | 0.0579 | 0.0560 |
| rae2822 | 0.2475 | 0.2730 | 0.2778 | 0.2592 | 0.1886 | 0.1073 | 0.0696 |
| biconvex | 0.0729 | 0.0803 | 0.0836 | 0.0802 | 0.0601 | 0.0551 | 0.0511 |
| dfvlr_r4 | 0.2674 | (rerun) | 0.2780 | 0.2688 | .. | .. | .. |

CD (naca0012): 0.036 -> 0.057 (peak wave drag at M0.925-0.95) -> 0.038 at M1.1.
Geometry spread at M0.85 is camber-driven: rae2822/dfvlr at 3.5x the
symmetric sections (measured camber 3.1% at x=0.4).

## TSD vs OF error percentages
| geo | M1.00 | M1.05 | M1.10 | note |
|---|---|---|---|---|
| naca0012 | -16.8% | -8.4% | -6.7% | symmetric profile, valid above sonic |
| biconvex | -15.9% | -10.3% | -4.6% | symmetric profile, valid above sonic |
| rae2822 | -74.1% | -55.3% | -31.8% | CAMBER NOT MODELED in TSD (feature gap) |
| dfvlr_r4 | -74.6% | -55.4% | -31.9% | CAMBER NOT MODELED in TSD (feature gap) |

TSD below M0.95 returns negative CL (degenerate branch) - DO NOT USE.
TSD above M1.05 for symmetric sections matches OF within 5-10%.

## Surrogate upgrade (wing3d/transonic.py)
The app predictor now uses the fine-mesh M-slice at a=1.25
(_ANCHOR_FINE, 11 Mach nodes) with legacy alpha-shape from the M0.8 row:
- M0.85 CL 0.0686->0.0728; M1.00 CD 0.0770->0.0419 (spurious coarse-mesh
  drag floor removed, -46%); M1.10 CL 0.0533->0.0560.
- Coverage extended 0.80->0.85 lower bound.

## Known limits (honest)
- Euler only (no BL/friction): CD is pressure+wave only.
- a=1.25 only in the fine table (alpha shape inherited from M0.8 legacy).
- shock_x/cpmin still legacy (not re-extracted from the fine campaign).
- TSD has no camber support: RAE2822/DFVLR-class sections need a
  section_model extension.
- dfvlr M0.900 thermo FPE (fixed with known-good thermo + maxCo 0.15),
  M0.975-M1.100 queued in the runner.
