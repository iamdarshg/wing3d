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

## Current run status (2026-10-23)

`tools/build_of3d_cases.py` prepares an M6 coarse case in `cases/m6_coarse/`.
The generated case is local/ignored and runs in container `oftrans` at
`/tmp/of3d_m6_queue/m6_coarse`. The owner-launched first attempt failed in
`snappyHexMesh` because OpenFOAM 2406 requires `minTwist` in
`meshQualityControls`; the active retry adds `minTwist 0.02` and reruns
`blockMesh -> snappyHexMesh -> checkMesh -> rhoSimpleFoam`. At restart,
`snappyHexMesh -overwrite` was running; no mesh or force result is yet
validated. Record the container logs, mesh statistics, y+, force histories,
and at least two refinement levels before changing anchor status. The flow
ctest has not been rerun in this restart; no probe values or mesh-uncertainty
bands are inferred from pending runs.
