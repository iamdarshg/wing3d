# HANDOFF — the anchors are inviscid; here is exactly what must change

Date 2026-10-30. Applies to `cases/openfoam/transonic/freetip/` (the free-tip
AR1.6 anchor set). Read `freetip/ANCHORS_INVISCID.md` for the numbers.

## 1. What we have

Compressible **Euler** results (`rhoSimpleFoam`, N2 JANAF, perfect gas, energy
equation on) on a free-tip AR1.6 wing. `CL` is usable as an inviscid anchor.
`CD` contains **wave drag only** — no skin friction.

Three settings make it inviscid:

| file | entry | current | must become |
|---|---|---|---|
| `0/U` | wing patch (`.stl`) | `type slip;` | `type noSlip;` |
| `constant/turbulenceProperties` | `simulationType` | `laminar` | `RAS` + `SpalartAllmaras` (or `kOmegaSST`) |
| `system/snappyHexMeshDict` | layers | `nSurfaceLayers 0` | `nSurfaceLayers 8-15` (see §3) |

Plus new fields that do not exist yet: `0/nut`, `0/nuTilda` (or `k`,`omega`),
and `system/fvOptions` with `limitTemperature` for robustness.

## 2. Why changing ONLY the wall type is not enough (measured)

The current near-wall mesh is far too coarse for any viscous boundary condition:

```
min cell volume = 2.15e-6 m^3   ->  wall cell height h ~ 0.0129 m  (1.3% chord)
Re = 6e6, U = 282 m/s, nu = 1.65e-5 m^2/s, u_tau ~ 10.9 m/s
y+ ~ h*u_tau/nu ~ 0.0129 * 10.9 / 1.65e-5  ~  8,500
```

* target for a **resolved** BL: `y+ ~ 1`  -> `h ~ 1.5e-6 m` (≈ **8,500x finer** than now)
* target for **wall functions**: `y+ ~ 30-100` -> `h ~ 4.5e-5 m` (≈ **300x finer** than now)

On the current mesh a `noSlip` wall evaluates `du/dy` across one 1.3%-chord cell,
i.e. across the entire boundary layer. The resulting wall shear is ~160x too
small (tau_computed ~ 0.8 Pa vs tau_true ~ 133 Pa at Cf = 0.003), so friction
drag would be badly under-predicted and the lift shift would be unphysical.
**A wall-type change alone produces a wrong viscous answer, not a viscous
answer.**

## 3. The ordered fix (XYZ)

**X. Wall functions first (fast, robust, defensible for CL/CD).**
1. `0/U`: wing -> `noSlip`.
2. `0/T`: wing -> `zeroGradient` (adiabatic).
3. `system/snappyHexMeshDict`: `nSurfaceLayers 8`, `expansionRatio 1.2`,
   first layer sized for `y+ ~ 50` (`h ~ 4.5e-5 m`), `minThickness 1e-6`,
   `nSmoothThickness 10`, `maxFaceThicknessRatio 0.5`.
   The earlier failure was layers collapsing at the **exposed free tips**
   (`deltaCoeffs -> inf`, killing every implicit solver). Mitigate by reducing
   the layer count near the tips (separate `layers {}` entry for the tip
   regions), raising `featureAngle`, or capping `maxThicknessToMedialRatio`.
4. `constant/turbulenceProperties`: `RAS { model SpalartAllmaras; }`.
5. New `0/nut`: wall `nutUSpaldingWallFunction` (or `nutkWallFunction`),
   farfield `calculated`.
6. New `0/nuTilda`: wall `fixedValue 0`, farfield `inletOutlet` ~ 3-5 nu.
7. `system/fvOptions`: `limitTemperature { min 200; max 4000; }`.
8. Verify with the `yPlus` function object that the wall y+ actually lands in
   30-300; if it does not, resize the first layer and repeat.

**Y. Resolved BL (y+ ~ 1) once X is stable.** Same as X but first layer
`h ~ 1.5e-6 m`, `nSurfaceLayers 15-20`; needs a finer surface mesh and will be
~2-3x more expensive. Use this for the final published anchors.

**Z. Re-run the whole set and report the viscous-vs-inviscid delta.**
| geometry | M | CL_inviscid | CL_viscous | dCL | CD_inviscid | CD_viscous | dCD |
(rows: dfvlr 0.85/0.925/1.00, rae2822 0.85/0.925/1.00, naca0012 0.80)

**That dCL/dCD table is the deliverable** — it is the size of the physics the
TSD potential solver omits, and it decides whether TSD-vs-OF comparisons are
meaningful at all.

## 4. Known traps (already paid for once)

* The repo `fvSchemes` is the **rhoSimpleFoam** variant (`steadyState` ddt +
  `SIMPLE`/`rhoSimple` blocks + relaxationFactors). Do **not** run
  `rhoCentralFoam` on these dicts: `"rho.*" { solver diagonal; }` plus
  `steadyState` gives the rho equation a **zero diagonal** -> FPE.
  (rhoCentralFoam also cold-starts badly: janaf T inversion. If you want it,
  converge rhoSimpleFoam first then `mapFields`.)
* `controlDict` in these cases had `deltaT 1`, which makes the steady loop do
  zero steps — verify `deltaT`/`endTime`/`writeInterval` before running.
* `rm -rf [0-9]*/` deletes `0/` (the initial fields). Restore with
  `git checkout -- cases/openfoam/transonic/freetip/`.
* Background jobs die on SSH exit — use `setsid nohup ... < /dev/null &`.
* VM: `gcloud compute ssh cfd-worker-1 --zone=asia-south1-c`. **Stop the VM
  when idle** (`gcloud compute instances stop ...`) — it bills while RUNNING.

## 5. Current owners

* Agent A18 is already working this (§3 X and Y) on the VM.
* The TSD side (A15 sweep, A17 Cp-estimator decontamination) is separate and
  unaffected: it needs the `dCL/dCD` numbers from §3 Z, not the raw inviscid rows.
