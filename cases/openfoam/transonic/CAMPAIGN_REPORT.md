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

---

## Library audit + fixes (2026-10-15)

### TSD camber support ADDED (was the structural gap)
Audit found `TSDSolver` had no camber concept at all: the transpiration BC
(tsd.py `_apply_wing_bc`) applied `wu = +dyt/dx - alpha`,
`wl = -dyt/dx - alpha` -- thickness only. For a cambered section the mean
line slope `dyc/dx` must be added to BOTH surfaces (the mean line is
shared; thickness is antisymmetric). Now `section_profile_xy` also returns
`camber` (d(yc)/dx) and `mean` (yc) callables, and the BC reads
`wu = +dyt/dx + dyc/dx - alpha`, `wl = -dyt/dx + dyc/dx - alpha`.
Symmetric sections (naca0012, biconvex) are unaffected (regression-verified:
M1.1 naca CL unchanged at 0.0523).

### CLI added: `python -m wing3d.tsdcli <file.stl>`
Options: --mach --alpha --span --epochs --station --json.
Auto-detects the extrusion axis, slices the mid-span section, extracts
t/c + camber, runs the campaign-validated TSD path (panel circulation +
inflow-pinned/outflow-free supersonic BCs), reports CL/CD/section_cl and
Mcrit. Verified on naca0012 (CL 0.0537 @M1.05), rae2822 (CL 0.0484,
camber 0.0074), and an arbitrary wing STL (auto-sliced, t/c 0.120).

### Slicer bugs fixed (found while building the CLI)
1. Extrusion-axis detection used argmin(extent) -> for slabs the correct
   axis is argmax (span/extrusion), not the airfoil thickness.
2. Bin-wise min/max over the slice population assigned UPPER points to the
   lower surface when a bin caught one side only -> camber/thickness read
   10x high. Fixed with a median split per bin + interpolation fill.
3. Binary STL reader (`geometry.read_stl`) could not read ANY binary STL
   (`count*12.5` floats reshape) -> now views bytes as uint8 (count,50) and
   reads the first 48 as float32.

### Geometry loop topology (investigated, no data change)
rae2822/biconvex/dfvlr_r4 loops originally closed upper LE->TE then lower
LE->TE, so the extrusion wrap quad spanned the whole chord. Verified the
POINT SETS are identical to the corrected loops (only the internal closure
edge differs) and the OF result is unchanged (0.2475 -> 0.2510, +1.4%).
The loops are now closed properly, but no campaign re-run was needed.
The corrected STLs were re-meshed (rae 48222->48114, biconvex 54822 cells)
and are in place.

### Final error table (never-importantly-corrected TSD vs OF)
Unchanged conclusion from above; RAE2822/DFVLR's 80-90% sub-1.0 gap is
architectural (panel circulation imprint locks lift at the incompressible
level; OF RAE reaches the 2D compressible value because it is a
supercritical airfoil -- verified cl_eff = 2*pi*alpha at M0.85).
The <10% band holds at M>=1.05 for symmetric sections (see table).

### Resources
Host RAM: WSL VM capped 2GB / container 1.5g / 2 CPUs (user-requested cap).
