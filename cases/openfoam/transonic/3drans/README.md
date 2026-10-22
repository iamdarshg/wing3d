# 3D RANS anchor assets and status

This directory contains the OpenFOAM 2406 3D compressible viscous anchor leg.
The TMR source assets are intentionally gitignored; retain their provenance
here rather than copying the large source bundle into Git.

## ONERA M6 provenance

Local sources are in `src/` (not tracked): `AileM6_with_sharp_TE.igs/.stp`,
`profile_root.dat`, `profile_onerad.dat`, `case_2308.dat`, and the TMR release
under `src/wing_release/wing_release_072319/`. The release includes
`hcf_wing_v5p0.f90`, `hcf_coarsening_v3p9.f90`, structured-grid namelists,
and the sharp-TE section `om6_wing_section_sharp.dat`. The locally generated
`stl/m6wing.stl` loft uses that section, root chord 1, semispan 1.4760179762,
30-degree LE sweep, and taper 0.5625159853. The small rounded CAD tip is
approximated as a capped loft tip; it is not an exact CAD tessellation.

## Current run status (2026-10-22)

`tools/build_of3d_cases.py` prepares an M6 coarse case in `cases/m6_coarse/`.
That generated case is local/ignored. It is queued in container `oftrans`
behind the active DFVLR `rhoCentralFoam` batch; no M6 mesh or force result is
available yet. Do not treat the case as a validated anchor until mesh
statistics, y+, force histories, and at least two refinements are recorded.
The requested `ctest_validate.py --suite flow --label of3d` currently reports
all four flow anchors PENDING. No probe values or mesh-uncertainty bands have
been inferred from the pending runs.
