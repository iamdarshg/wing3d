# OpenFOAM fidelity investigation (A7)

## Status: no root cause isolated yet

The supplied campaign result is reproducibly steady at the reported level: the last forceCoeffs row at `t=0.25` gives `Cl=0.06657627` (and `Cd=0.02333923`). This confirms the reported ~0.0666 result; it does not identify the cause. The comparison target remains the 2D inviscid `Cl=0.352` (Vassberg & Jameson; UT-IAS DG `0.35169`). A finite-span wing must not be compared directly to this 2D value.

## Audited inputs and hypotheses

- Case uses `rhoCentralFoam`, `localEuler`, Kurganov flux and vanLeer reconstruction. No controlled alternate-solver or refinement run was completed in this investigation, so H1/H2 remain **untested**, not refuted.
- Wall boundary conditions are `U: slip`, `p: zeroGradient`, `T: zeroGradient`; the force patch is `naca0012.stl`, with `liftDir=(0 0 1)`, `rhoInf=1.1231`, `magUInf=282.446`, and `Aref=1.6`. The configuration appears internally consistent, but no independent wall-pressure integration was performed; H3 is **not fully tested**.
- Thermophysical properties specify perfect-gas N2 with JANAF data and molecular weight `28.01348`; initial/farfield values are `p=100000 Pa`, `T=300 K`, `U=(282.379,0,6.162)`. The resulting speed is `282.446 m/s`, and `M=0.800` if `a=353.06 m/s` as reported. `p` is prescribed through `freestreamPressure` with `freestreamValue=100000`, i.e. static pressure, not total pressure. The isentropic total/static ratio at M=0.8, gamma=1.4 is `1.52434` (`p0≈152434 Pa`); this is not a BC mismatch because the configured value is static. H4 is **consistent with the intended setup**, but the actual solved-field Mach and JANAF sound speed were not independently sampled.

## Conclusion / action

No confirmed root cause or minimal fix can be responsibly stated from these checks. The observed coefficient is a serious fidelity failure; therefore the campaign anchors based on it should be treated as invalid and regenerated **after** a controlled solver/refinement test identifies and verifies a correction. Do not attribute the deficit to mesh, dissipation, or BCs on the evidence currently available. No campaign deliverables were modified.
