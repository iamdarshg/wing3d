# MASTER HANDOFF — wing3d TSD vs OpenFOAM
**Written 2026-10-30.** Everything needed to resume this program from a cold
start, including the viscous+turbulent anchor work that is the current
priority. Read order: §1 status → §2 the viscous spec → §3 the traps.

---

## 0. North star (durable, set 2026-10-22)

1. **Speed**: wing3d TSD must be **3–5× faster than OpenFOAM** per
   (geometry, Mach) case at matched accuracy.
2. **Accuracy**: our CL uncertainty bound — derived from **our** hyperparameters
   (grid {80,120,140,160} × epochs {8,16,24} × mask state × k-imprint) — must
   lie **within OpenFOAM's own bound**, ideally centred on truth.
3. **Claim a bound, not a point.** Every deliverable reports CL ± band.
4. **Honesty rule**: physically unsupportable fixes that improve the score are
   logged as BAD BUT HELPS with their measured cost; anything whose removal
   lowers error is removed permanently.

---

## 1. Verified state (what is true today, with numbers)

### 1.1 The solver (`wing3d/`)

Transonic full-potential, irrotational by construction, with three rotational
terms added:

| term | commit | measured effect |
|---|---|---|
| shock entropy (non-isentropic density) | `06a4bdd` `WING3D_ENTROPY=1` | rae2822 M0.85: **−40.3% → −17.2%**; dfvlr M0.85 −4.0 → −2.1 |
| wall BL (integral BL as transpiration) | `1ecfa9e` `WING3D_BL=1` | **live but ~neutral** (+0.03…+0.5% CL) — see 1.4 |
| coupled bound-vortex Kutta | `4d53213`,`9fcf99b`,`ae137d8`… `WING3D_KUTTA=6/7` | structurally **correct** (hybrid row degrades to the kernel at 2.6e−16; φ-liveness 25.5 vs 0.0), but the **Mach response is 5× too weak** → CL stays flat |

Other fixes that matter:
* `c46ce46` — replaced the **broken 3D Morino panel imprint** with a Prandtl
  lifting-line imprint. That single fix moved the mean error 72.8% → 25.1%.
  (The panel solve returned 2D-like circulation at every span: CL +0.63 /
  **−7.50** / −1.28 at span 1.6/3.2/6.0.)
* `9699435` — Prandtl–Glauert capped at M_crit (it was inflating lift as M→1).
* `66142d9` — removed a `dy` double-count on the branch cut (`muw` is already
  Γ₂D; the repo's own `CL = 2∫muw dy/span` proves it).
* `WING3D_CPDECONTAM=1` (A17) — decontaminated the surface-Cp estimator:
  **dz spread 22.5% → 2.7% (8.3×)** on grids whose sheet is within half a cell
  of z=0; rigid 0.25dz offset −9.7% → −0.1%.
  **This is the single largest numerical win in the project** — it removed a
  10–27% discretisation floor that had been masking every physical result.

### 1.2 Current TSD accuracy (against the free-tip anchors, ±5.7% band)

| geometry | M | TSD CL | anchor CL | err |
|---|---|---|---|---|
| dfvlr_r4 | 0.850 | 0.3143 | 0.271861 | **+15.6%** |
| dfvlr_r4 | 0.925 | 0.3185 | 0.262529 | +21.3% |
| dfvlr_r4 | 1.000 | 0.3287 | 0.226938 | **+44.8%** |

(Against the *old* fine-mesh anchors the same runs gave −3.8/−12.9/+56.3%; the
anchor changed, not the solver. Numbers in
`cases/openfoam/transonic/freetip/ANCHORS_INVISCID.md`.)

**The failure mode is unambiguous: CL is flat across Mach while the truth rises
then collapses (−16.5% from M0.85→M1.00). The solver does not predict
transonic lift loss. That is the #1 open physics problem.**

### 1.3 Anchors (OpenFOAM) — and the two defects found

* **Defect 1 (fixed)**: the historic campaign slabs had symmetry on the **z**
  faces, so the wing tips sat in `freestream`/`freestreamPressure` planes —
  tips clamped by a *pressure* boundary, neither an isolated wing nor a clean
  periodic cascade. Rebuilt as **free-tip AR1.6** in
  `cases/openfoam/transonic/freetip/` (wing y=±0.8, domain y=±4, all outer
  faces farfield, no symmetry). Tip unloading proven (dCp drops 52% at η=0.81
  vs the clamped control which *gains* 96% at the boundary).
* **Anchor band ±5.7%** (span ±5.2% ⊕ domain ±2.1% ⊕ time ±1.0%), measured.
  The 22% rhoSimple-vs-rhoCentral spread was reclassified as a **definition
  choice**, not uncertainty.
* **Defect 2 (OPEN, the priority)**: these anchors are **compressible EULER**,
  not viscous. See §2.

### 1.4 Why the BL term is inert (architectural, understood)

With `free_k=False` the bound circulation is *prescribed*; a boundary layer can
redistribute pressure but cannot reduce the bound circulation — the main route
by which wall vorticity lowers CL. Measured: BL moves CL by +0.03…+0.5%.
It can only bite once the circulation is genuinely free (§3 P2).

### 1.5 Infrastructure

* **GCP `cfd-worker-1`** — `asia-south1-c`, `e2-standard-8` (8 vCPU/32 GB),
  80 GB disk, ~$0.27/hr **on-demand** (preemptible quota is 0 in this project).
  OpenFOAM 2406 (ESI, `/usr/lib/openfoam/openfoam2406`), gfortran, python3+numpy/scipy.
  **Stop when idle**: `gcloud compute instances stop cfd-worker-1 --zone=asia-south1-c`.
  Runbook + 7 gotchas: `tools/of_mesh/vm/README.md`.
  The local Docker/WSL container was stopped 2026-10-29 because it was
  saturating the host (WSL 1.8 GB → 254 MB; free RAM 0.79 → 2.17 GB).
* **Mesher**: NASA TMR `hcf_wing` compiles and generates grids (M6 verified:
  ref area 1.1531508, MAC 0.8016730, AR 3.7786), but both conversion routes
  (PLOT3D and a hand-written SU2→polyMesh) produce **folded cells**
  (max skewness 183, 265 bad faces) — the TMR *coarse* grids are genuinely that
  poor. Blocked; see §3 P4.

---

## 2. PRIORITY 1 — viscous + turbulent anchors (the user's explicit ask)

**The anchors are inviscid.** Current wall treatment: `0/U` wing = `slip`,
`simulationType laminar`, `nSurfaceLayers 0`. Full detail and the ordered fix
in **`cases/openfoam/transonic/VISCOUS_HANDOFF.md`**. Summary:

### 2.1 Why a wall-type change alone is NOT enough (measured)

```
wall cell h ≈ 1.3% chord (min volume 2.15e-6 m^3)
Re 6e6, U 282 m/s, ν 1.65e-5, u_τ ≈ 10.9 m/s
y+ ≈ h·u_τ/ν ≈ 8,500
resolved BL needs y+ ≈ 1     → h ≈ 1.5e-6 m  (8,500× finer)
wall functions need y+ ≈ 30–100 → h ≈ 4.5e-5 m (300× finer)
```
On today's mesh a `noSlip` wall differentiates across one cell that contains the
**entire** boundary layer: τ_computed ≈ 0.8 Pa vs τ_true ≈ 133 Pa → friction
**~160× too low**. A wall-type-only change yields a *wrong* viscous answer.

### 2.2 The ordered plan (X → Y → Z)

**X — wall functions (fast, robust, defensible for CL/CD).**
1. `0/U` wing → `noSlip`; `0/T` wing → `zeroGradient` (adiabatic).
2. `system/snappyHexMeshDict`: `nSurfaceLayers 8`, `expansionRatio 1.2`,
   first layer sized for `y+ ≈ 50` (h ≈ 4.5e-5 m), `minThickness 1e-6`,
   `nSmoothThickness 10`, `maxFaceThicknessRatio 0.5`.
   **The known failure**: layers collapse at the **exposed free tips** →
   `deltaCoeffs → ∞` → every implicit solver dies. Mitigate with a separate
   thinner `layers {}` entry for the tip regions, a larger `featureAngle`, or a
   capped `maxThicknessToMedialRatio`.
3. `constant/turbulenceProperties`: `simulationType RAS; RAS { model SpalartAllmaras; }`.
   **This is the turbulent part the user wants.** (`kOmegaSST` is the alternative.)
4. New `0/nut`: wall `nutUSpaldingWallFunction` (or `nutkWallFunction`),
   farfield `calculated`.
5. New `0/nuTilda`: wall `fixedValue 0`, farfield `inletOutlet` ≈ 3–5ν.
6. `system/fvOptions`: `limitTemperature { type limitTemperature; min 200; max 4000; }`.
7. Verify with the `yPlus` function object that wall y+ lands in **30–300**;
   resize the first layer and repeat if not.

**Y — resolved BL (y+ ≈ 1)** for the published anchors: first layer
h ≈ 1.5e-6 m, `nSurfaceLayers 15–20`, finer surface mesh, ~2–3× cost.

**Z — the deliverable**: re-run the set and publish the
**viscous-vs-inviscid ΔCL/ΔCD** table (dfvlr 0.85/0.925/1.00, rae2822
0.85/0.925/1.00, naca0012 0.80). That Δ is the size of the physics the
potential solver omits, and it decides whether any TSD-vs-OF comparison is
meaningful.

**Who**: agent A18 was KILLED 2026-10-30 after ~90 min of mesh trials that
produced **zero layers** (every trial mesh was bit-identical to the inviscid
one: 136,158 cells, min volume 2.1495911e-06) while snappy still printed
`Layers added`. **Do not repeat the prism-layer loop**; use option A
(refine instead of layer) or B (`hcf_wing` igrid_type=4 + `su2_to_of.py`).
Full post-mortem: `cases/openfoam/transonic/VISCOUS_HANDOFF.md` §6.
**Turbulence note**: at Re ≈ 6×10⁶ the real boundary layer is turbulent, so
`laminar` is unphysical *even with a no-slip wall* — SA/SST is required, not
optional.

---

## 3. Remaining priorities

### P2 — TSD transonic lift loss (the #1 physics gap)
CL is flat (0.314/0.319/0.329) while truth rises then collapses
(0.272/0.263/0.227). Established:
* Circulation must become **genuinely free**: A13's hybrid tangency row is
  structurally correct and φ-coupled; what remains is that the Mach signal is
  ~5× too weak and non-monotone, swamped by non-circulatory content.
* The field's sheet-response error is a **chordwise resolution** effect
  (field/kernel ratio 1.59 → 1.48 → 1.39 at nx 80/120/160, **no** dz_core
  dependence) — so it converges with `nx`, i.e. it is cells, not a rewrite.
* Next concrete step: run the coupled Kutta at **nx ≥ 160** and with the
  A17 decontaminated estimator, and check whether CL now falls monotonically
  ~15–25% from M0.85→M1.00. If it does, the solver is genuinely transonic.

### P3 — planarity / non-planar geometry (A16 findings)
| effect | ΔCL/CL |
|---|---|
| dihedral 5°/10°/20° | −0.31 / −1.26 / −5.0% |
| spanwise arch h=0.1c | −0.43% |
| **current model's plane-thickness/snap error** | **±10–27%** (now largely
  fixed by `WING3D_CPDECONTAM`) |

Implemented and bit-identical at zsurf=0 (`set_surface`, per-(i,j) sheet/source,
`WING3D_NONPLANAR`), but **dihedral currently fails**: warping the branch cut
splits the dipole across two faces and contaminates the ±1-plane estimators
(mode 7 walks to negative circulation). Per-span sections (taper/sweep/crank/
centrebody) **work** — a rectangular wing built that way reproduces stock at
ΔCL = 0.0 — but a strongly swept planform has an **inconsistent frame** between
the absolute-x chordwise stations and the local-chord kernel cells
(CL_circ collapsed to 0.0034). Fix identified: body-fitted sheet.

### P4 — mesher (blocked)
`hcf_wing` generates NASA-quality grid *definitions*, but our conversions fold
cells (skew 183). Options, in order: (a) tune `nre`/`R_outer`/`stretching`
toward the production TMR regime; (b) **pyHyp** (MDO Lab, hyperbolic extrusion,
explicit orthogonality control, CGNS → `cgnsToFoam`); (c) for 2D only,
`airfoilgmesher` (Gmsh). Do not trust any RANS result from this path until
`checkMesh` passes.

### P5 — speed target (3–5×)
**Never measured.** Needed: wall-clock per (geometry, Mach) for TSD vs OF at
*matched accuracy*, plus the UQ band (grid × epochs × mask × k-imprint = 48
combos). The fast-epoch solver (`tests3d/tsd_fast.py`: coloured-FD Jacobian +
GMRES) reaches **12 s/epoch at g80 vs ~60 s** for campaign iterations, but
stalls at epochs 1+ because the mask-refreshed operator is pathologically
non-normal — that is the blocker to the 30–45× path identified in the perf
audit.

---

## 4. How to run things

```powershell
# VM
gcloud compute ssh cfd-worker-1 --zone=asia-south1-c
# (copy files: gcloud compute scp <local> cfd-worker-1:/home/iamdarshg/<name> --zone=asia-south1-c)
#   NOTE: '~' is NOT expanded by pscp. Use absolute paths.
bash /home/iamdarshg/vm_status2.sh      # progress table

# TSD (local, ~300 MB)
$env:WING3D_IMPRINT='ll'; $env:WING3D_ENTROPY='1'
python -X utf8 -c "... run_sup(M, span=1.6, epochs=8, free_k=False, ...)"

# validation harness
python -X utf8 tests3d/validation/ctest_validate.py --suite geometry
python -X utf8 tests3d/su2_to_of.py <grid.su2> <case_dir>   # blocked: folded cells
```

## 5. Traps already paid for (do not repeat)

1. `scp` does not expand `~` → absolute `/home/iamdarshg/...`.
2. **Meshes and STLs are not committed** — each case needs
   `constant/triSurface/<name>.stl` copied from `cases/openfoam/transonic/stl/`
   before `snappyHexMesh`; the STL name is read from the case's own
   `snappyHexMeshDict`.
3. Background jobs die on SSH exit even with `nohup` → `setsid nohup … < /dev/null &`.
4. The repo `fvSchemes` is the **rhoSimpleFoam** variant
   (`ddtSchemes { default steadyState; }` + `SIMPLE`/`rhoSimple` blocks). Run
   `rhoSimpleFoam` on them. `rhoCentralFoam` gets a **zero diagonal** for rho
   (`"rho.*" { solver diagonal; }` + steadyState) → FPE on step 1, and
   cold-starts badly (janaf T inversion; converge rhoSimpleFoam then `mapFields`).
5. `controlDict` in the freetip cases had `deltaT 1` → the steady loop takes
   **zero steps** (`Starting time loop` then `End`).
6. `rm -rf [0-9]*/` **deletes `0/`** (the initial fields). Restore with
   `git checkout -- cases/openfoam/transonic/freetip/`.
7. Windows PowerShell mangles `$(...)`, `$env:`, `|` and heredocs passed
   through `gcloud compute ssh --command`. **Write scripts to a file, scp them,
   run them** — that is the only reliable pattern.
8. Do not run heavy container/VM work while the host is low on RAM; the user's
   box has ~15.2 GB total and browsers alone hold 2–3 GB.

## 6. File map

| path | what |
|---|---|
| `wing3d/tsd.py` | solver: residual, transpiration BC, entropy, Kutta modes, Cp decontamination |
| `wing3d/ibl.py` | integral boundary layer (Thwaites/Head) |
| `tests3d/tsd_m_sup.py` | `run_sup`: epoch loop, state vector, Kutta coupling |
| `tests3d/tsd_coupled.py` | `build_case`: panel/imprint setup |
| `tests3d/validation/` | ctest harness, `published_refs.json`, `known_anchors.json` |
| `cases/openfoam/transonic/freetip/` | free-tip AR1.6 anchors (**inviscid**) |
| `cases/openfoam/transonic/freetip/ANCHORS_INVISCID.md` | the numbers, marked inviscid |
| `cases/openfoam/transonic/VISCOUS_HANDOFF.md` | **the viscous+turbulent spec** |
| `cases/openfoam/transonic/BAD_BUT_HELPS.md` | unsupportable-but-helpful ledger |
| `cases/openfoam/transonic/SOLVER_FINDINGS_2026-10-27.md` | root-cause doc |
| `tests3d/a9..a13_kutta_report.md` | the Kutta development trail |
| `tools/of_mesh/` | mesher pathway + VM runbook + scripts |
| `tests3d/of_mesh.py` | host-side mesher orchestrator |
