# Of-mesh pathway — body-fitted grids for the viscous anchor leg

Status: **permanent pathway** (set 2026-10-28). Replaces the ad-hoc
Cartesian + `snappyHexMesh` route, which could not resolve the LE suction
peak (measured: NACA0012 M0.80 α1.25 in the old template gave CL = 0.0666
against the published 0.352 — a 5.3x lift deficit).

## Primary: NASA TMR `hcf_wing` (structure O-grid, y+ target)

Source: NASA Langley `wing_release_072319` (TMR). Archived in the repo at
`cases/openfoam/transonic/3drans/src/wing_release/wing_release_072319/`
(gitignored as bulky assets; the compiled binary lives in the container at
`/work/tmr/wing_release/wing_release_072319/hcf_wing`).

Why this one:
- NASA-authored and NASA-validated; it is the generator behind the TMR
  reference grids (ONERA M6, NACA0012 wing).
- Supports our sections directly: `airfoil_type=1, nacaxx=12` for
  NACA00XX, or `airfoil_type=0, airfoil_data_file="<selig>.dat"` for
  rae2822 / dfvlr_r4 / biconvex.
- Boundary-layer clustering to a y+ target; structured O-grid
  (`igrid_type=5`) or prism/hex (`igrid_type=4`).
- Emits the reference quantities we need: reference area, MAC,
  aspect ratio (verified for M6: area 1.1531508, MAC 0.8016730, AR 3.7786).

### Verified conversion recipe (two non-obvious gotchas)

1. `ugrid_file_unformatted = F` in `input.nml`, so the generator writes a
   **formatted** `.p3d`. Its default `.ufmt` is binary and OpenFOAM cannot
   read it.
2. `plot3dToFoam -noBlank <grid>.p3d` — the `-noBlank` flag is required,
   otherwise the reader consumes an iblank array that is not present and
   dies at EOF.

Then `autoPatch 30 -overwrite`. The three resulting patches are
identifiable by geometry and must be renamed:

| patch | test | assign |
|---|---|---|
| wing | z-range thin (<= ~0.06) and x in [0, ~1.5], y in [0, semispan] | `wall` (noSlip for RANS) |
| symmetry | **all** points at y = 0 | `symmetry` (root plane) |
| farfield | everything else (radius ~ R_outer) | `patch` (freestream) |

Verified on the M6 grid: 41x65x49 nodes -> 122,880 cells, 124,865 points,
boundary auto0/auto1/auto2 = wing/symmetry/farfield. checkMesh: all
topology checks pass, max non-orthogonality 86.6 deg (average 34.0),
3 highly skew faces at the O-grid tip pole (known artefact, accepted).

## Fallback: pyHyp (MDO Lab)

`mdolab-pyhyp` — hyperbolic extrusion from a surface mesh. Used in the
MACH/ADflow stack; produces the same O-grid + viscous-layer topology and
exports CGNS. Adopt it if `hcf_wing` cannot represent a section (e.g. a
blunt/odd TE) or if the PLOT3D conversion ever regresses.

## Quick 2D validation meshes: gmsh / airfoilgmesher

For fast 2D checks against published airfoil data (`pip install
airfoilgmesher`, Gmsh-backed) — useful to validate a solver against
CL = 0.352 (NACA0012, M0.8, a1.25) before spending time on 3D.

## Runner

`tests3d/of_mesh.py` (host side) orchestrates:
`input.nml -> hcf_wing (docker) -> plot3dToFoam -noBlank -> autoPatch ->
patch rename by geometry -> case skeleton`.

Usage:
```
python tests3d/of_mesh.py --airfoil naca0012 --semispan 1.6 --nspan 32 \
    --nradial 64 --yplus 1 --re 6.0e6 --out /work/cases/naca_ar16
```

## Status 2026-10-28 (end of session)

Working: hcf_wing compiles and generates grids; input.nml -> .p3d
(formatted) -> plot3dToFoam -noBlank -> utoPatch -> patch
classification/rename (wing / symmetry / farfield, verified by bbox) ->
case skeleton (
hoSimpleFoam + SpalartAllmaras, N2 janaf thermo,
turbulenceProperties RAS) all run end-to-end via tests3d/of_mesh.py.

BLOCKER: the converted meshes are FOLDED -- checkMesh reports
  max aspect ratio 212769, max skewness 183, 265 highly skew faces
  (Failed 2 mesh checks)
and rhoSimpleFoam aborts at iteration ~6 on a janaf thermo inversion.
Identical numbers are produced with NASA's OWN tuned input
(input.nml_naca0012_round verbatim) as with our parameterisation, so this
is NOT a parameter problem -- it is the PLOT3D -> OpenFOAM conversion
(either plot3dToFoam mishandles this file's ordering, or an iblank/2D
convention). Next step: convert from the SU2 export instead
(wing_strct.1.su2 is written by the generator) with a direct SU2->polyMesh
writer, and re-run checkMesh. Do not trust any RANS result from this
pipeline until checkMesh passes.

## CORRECTION 2026-10-28 (retracts a wrong claim above)

The earlier statement "identical bad numbers with NASA's own input, so it is
the conversion" was WRONG -- it was a STALE FILE artefact. The generator
writes a fixed filename; `input.nml_naca0012_round` uses `igrid_type = 4`
(prism/hex, UNSTRUCTURED) and therefore writes no .p3d at all, so the .p3d
converted the second time was still the M6 grid from the earlier run.

Corrected facts:
* `igrid_type = 5` (structured) -> .p3d   (the M6 case works this way)
* `igrid_type = 4` (prism/hex)   -> .ugrid/.su2 only (NACA0012 built-in input)
* ALWAYS `rm -f wing_strct.1.p3d` before regenerating, and check the
  Reference area / Aspect ratio line -- it identifies which grid you got:
    M6      : area 1.1531508, AR 3.7786
    NACA0012: area 1.5000000, AR 3.0000
* M6 structured .p3d -> plot3dToFoam: topology checks PASS, but 265 highly
  skew faces (max skew 183) remain and rhoSimpleFoam aborts at iteration ~6.
  The folded-cell problem is real for the structured PLOT3D route.
* `tests3d/su2_to_of.py` (direct SU2 -> polyMesh) parses correctly
  (124865 points / 122880 hexes; markers WING3D 1920, SYMMETRY 3072,
  FARFIELD 1920) but yields 3072 illegal cells and 65 unused points,
  i.e. face-key collisions -- most likely degenerate quads on the symmetry
  plane (zero-thickness first layer). Needs a degeneracy guard: skip faces
  whose sorted key already exists but with a different vertex multiset, and
  drop zero-area quads before deduplication.

## Update 2026-10-28 (later): mixed elements + tip variants

* `tests3d/su2_to_of.py` now handles **tets/prisms/pyramids as well as
  hexes**. That was the cause of the "3072 illegal cells": those cells are
  prisms/hex-mixed (igrid_type=4), and a hex-only converter left them with
  no faces at all. With mixed support, face counts are correct
  (370560 faces, 6912 boundary = 1920+3072+1920 exactly matching markers).
* The SU2 route and the PLOT3D route give the SAME checkMesh numbers, so
  the poor quality is a property of the generated grid, not the reader.

### Grid quality matrix (all after successful conversion)

| config | cells | max aspect | max skew | faces >70 deg | checkMesh |
|---|---|---|---|---|---|
| M6, round tip, structured | 122880 | 212769 | 183.1 | 19738 | FAIL 2 |
| NACA0012, round tip, structured | 122880 | 212770 | 183.1 | 19738 | FAIL 3 |
| NACA0012, blunt tip, structured | 245760 | 17184 | 92.3 | 54272 | FAIL 2 |

Blunt tip fixes aspect/skew by ~12x but doubles the non-orthogonal face
count. rhoSimpleFoam aborts at iteration ~6 on all of them: pressure runs
to -950 kPa in a single cell (janaf inversion), i.e. one irreducibly bad
cell region rather than gradual divergence.

### What is NOT the problem (verified)
* not the reader (two independent readers agree)
* not negative volumes (cell volumes OK)
* not the patch naming/BC set (all patches covered, correct types)
* not the thermo file (identical to the working campaign case)

### Recommended next steps (in order)
1. Tune the generator for orthogonality: larger `nre`, smaller
   `stretching_tanh_towards_lete`, larger `R_outer` (M6 reference uses
   R_outer=100, nr_gs=8, nre=64). The production TMR grids are known good;
   our parameterisation is not.
2. If tuning stalls: switch to **pyHyp** (MDO Lab), which extrudes a
   hyperbolic O-grid from a surface mesh with explicit orthogonality
   control, and export CGNS -> OpenFOAM.
3. For 2D checks only: `airfoilgmesher` (Gmsh) is the fast route to a clean
   airfoil mesh to validate the solver against CL=0.352.
