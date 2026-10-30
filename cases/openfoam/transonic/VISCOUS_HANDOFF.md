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
4. **Turbulence model — chosen on evidence, not preference: k-ω SST primary,
   Spalart–Allmaras as cross-check.** Rationale (Kaya 2025, *CMES* 145(1),
   doi 10.32604/cmes.2025.072000 — RAE2822 at M=0.74, Re=2.7e6 with
   **rhoCentralFoam/TSLAeroFoam in OpenFOAM**, i.e. our exact toolchain):

   | model | measured finding |
   |---|---|
   | **k-ω SST** | "most robust across all cases"; shock position within **<0.02c** of experiment; held accuracy even at higher y⁺ |
   | Spalart–Allmaras | "broadly comparable to k-ω SST"; more robust numerically (one equation) |
   | k-ω Realizable | strongly y⁺-sensitive; shock displaced downstream |

   The gap is small — SST is chosen for the best *shock position*, not because
   SA is poor (NASA TM-110446 validated both on RAE2822; NASA TMR uses SA-neg
   for the ONERA M6 reference). **Report both if budget allows** — that is the
   turbulence-model half of the anchor band, handled exactly like the
   rhoSimple/rhoCentral spread was.
   `constant/turbulenceProperties`: `simulationType RAS; RAS { model kOmegaSST; }`.
   SST wall treatment: prefer **y⁺ ≲ 1**; SST with wall functions is less
   accurate than SA with wall functions, so use SA at the X-stage wall-function
   mesh and switch to SST on the Y-stage resolved mesh.
5. New `0/nut`: wall `nutUSpaldingWallFunction` (or `nutkWallFunction`),
   farfield `calculated`.
6. **SST fields** `0/k` (wall `kqRWallFunction` or `fixedValue` ~1e-10;
   farfield ~1e-6 U²) and `0/omega` (wall `omegaWallFunction`; farfield
   ω ≈ 5–10 U/L). **SA instead**: `0/nuTilda` wall `fixedValue 0`,
   farfield `inletOutlet` ~ 3–5ν.
7. `system/fvOptions`: `limitTemperature { min 200; max 4000; }`.
8. Verify with the `yPlus` function object that the wall y+ actually lands in
   30-300; if it does not, resize the first layer and repeat.

**Y. Resolved BL (y+ ~ 1) once X is stable.** Same as X but first layer
`h ~ 1.5e-6 m`, `nSurfaceLayers 15-20`; needs a finer surface mesh and will be
~2-3x more expensive. **This is where k-ω SST belongs** (SST with wall
functions is weak; SST on a resolved y⁺≈1 mesh is what the Kaya 2025 study
found to give shock positions within 0.02c). Use this for the final published
anchors, with SA run once as the model-spread cross-check.

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

---

## 6. A18 post-mortem — prism-layer route BLOCKED (2026-10-30)

**Agent A18 was killed 2026-10-30 after ~90 minutes in mesh trials, having
produced no viscous case.** Recorded here so nobody repeats the loop.

### What it did
Eight snappy trials under
`~/cfd/cases/openfoam/transonic/freetip_viscous/trial_{t1,t2,t3,t4,t8,t9,ta,tools}`,
varying layer count / first-layer thickness / mediation parameters.

### What it found (the finding that matters)
**Every trial produced a mesh IDENTICAL to the inviscid one.**

| trial | cells | mesh | nSurfaceLayers req. | min volume |
|---|---|---|---|---|
| t1, t2, t3, t4, ta | **136,158** | OK | 4 | **2.1495911e-06** |
| inviscid reference | 136,158 | OK | 0 | 2.1495911e-06 |

Same cell count and same min volume to 8 digits => **zero layers were actually
added**. snappy prints `Layers added in = 3.5-7.5 s` regardless, which is
misleading: the free-tip layer collapse is manifesting as a silent no-op, not a
crash. `checkMesh` therefore reports `Mesh OK` on what is really the inviscid
mesh, and a solve on it would silently look "viscous" while being Euler.

**Do not trust `Mesh OK` as evidence that layers were added.** The check is:
`nSurfaceLayers > 0` AND the cell count / min volume must CHANGE from the
inviscid mesh. If they are identical, the layer request was a no-op.

### Unblock options (ordered)

**A. Refine instead of layer (recommended; ~30-45 min to first viscous solve).**
Skip prism layers. Raise snappy surface refinement by 3-4 levels so the wall
cell goes 1.3e-2 -> ~5e-5 (y+ ~ 30-100), then `noSlip` + SA with
`nutUSpaldingWallFunction`. A hex-refined near-wall region is less elegant than
prisms but is valid for wall functions and dodges the tip-collapse failure
entirely. Cost: mesh grows 4-8x (~0.5-1M cells), ~1 h/case on 8 cores.
This is the configuration the Kaya 2025 study calls acceptable for steady
SA/SST.

**B. `hcf_wing` prism grid + SU2 conversion (~2-3 h, better quality).**
`igrid_type=4` produces prism boundary layers on a hex core with a y+ target,
which is exactly the desired topology. `tests3d/su2_to_of.py` now handles
mixed prisms/tets/hexes correctly (that was the "3072 illegal cells" fix), so
the conversion route exists. Open risk: the folded-cell/orthogonality problem
that made the *structured* route unusable may or may not recur for prisms.

**C. If both fail**: accept that a proper BL mesh needs a body-fitted generator
(`pyHyp`, MDO Lab: hyperbolic extrusion with explicit orthogonality control,
CGNS -> `cgnsToFoam`). That is the standing recommendation in
`tools/of_mesh/README.md`.

### Where the state is
Trials left on the VM at
`~/cfd/cases/openfoam/transonic/freetip_viscous/trial_*/` for inspection.
The VM is idle (`load 0.03`), no jobs running. **Stop it when done**
(`gcloud compute instances stop cfd-worker-1 --zone=asia-south1-c`) — it bills
while RUNNING.
