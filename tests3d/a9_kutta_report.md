# A9 — free bound-vortex circulation closed by Kutta: what was built and what it measured

Status: **negative result on the specified formulation, with a structural
diagnosis.** The circulation *can* be made an unknown (the plumbing works, the
Newton system converges), but a Kutta row built on the TE loading cannot *select*
it when the sheet's loading shape is prescribed. Numbers below; nothing is tuned.

Gates: `WING3D_KUTTA=1` (formulation as specified), `=2` (classical free wake),
`=3` (bound sheet + tangency). Default 0 = baseline, bit-identical.

## 1. Formulation implemented

`wing3d/tsd.py`
* `flat_plate_jump(x) = (2/pi)[sqrt(x(1-x)) + arcsin sqrt(x)]` — the cumulative
  flat-plate loading, `F(0)=0`, `F(1)=1`. Its derivative is the classic
  LE-singular `A(x) = (2/pi)sqrt((1-x)/x)` normalised to `int A dx = 1`, so the
  sheet amplitude **is** the section circulation. The branch cut carries the
  *potential jump*, so what is applied on the wing is the integral `F`, not `A`
  (this also makes the TE jump equal the trailed wake jump, i.e. Kutta/Helmholtz
  by construction).
* `gamma_b[j]` (per-span circulation) is a new unknown — named `gamma_b`, not
  `Gamma`, because `self.gamma` is the ratio of specific heats and a case
  collision would be silent. It is applied as the bound sheet `F(x)*gamma_b[j]`
  on the wing and continued unchanged down the wake (Helmholtz, constant in x;
  the `_circff` farfield blending is disabled so nothing re-pins it).
* Blood-out gate: `_apply_wing_bc` drops `slope_c - alpha` (thickness only) in
  mode 1, keeping `+slope_t` / `-slope_t`.
* `tests3d/tsd_m_sup.run_sup`: state vector `y = [phi[1:].ravel(), gamma_b]`,
  both unpacked per residual call; the stale-state clearing at entry resets
  `gamma_b`; rows `kutta_residual(upstream=2) == 0` appended. Initial guess from
  Prandtl lifting line (`tests3d/tsd_coupled.lifting_line_gamma`, the same
  function the `WING3D_IMPRINT=ll` path now calls, so the two cannot drift).
  The build_case (LL-imprint) interior is kept as the initial phi.

## 2. Row weight: measured, not 10

`tests3d/tsd_kutta_weight.py` on the real g80x32x32 grid, M0.3:

| quantity | value |
|---|---|
| `A_scale` = dR_phi/dphi (median node diagonal, 5 nodes, k0±1) | 6.33e3 |
| `C_scale` = dR_kutt/dphi (extrapolation planes k0±1, k0±2) | 2.10e1 |
| `B_scale` = dR_phi/dgamma (sheet source) | 2.90e3 |
| direct dR_kutt/dgamma at frozen phi | **0** (the rows read phi only) |
| upper bound `C*B/A` | 9.6 |
| **converged coupled gain** `|dCp(TE)/dgamma|` (two NK solves) | **1.18e2** (mode 1), **4.72** (mode 2) |

The Kutta rows have no direct `gamma` term (they are built from phi), so the
honest measure is the *converged* gain, not a frozen-phi derivative. Weight
`w = A_scale/|gain|` = 54 -> **50** (modes 1/3), and 1.34e3 -> **1500** (mode 2,
much weaker channel). Env override `WING3D_KUTTA_W`.

## 3. Ladder A — low Mach, lifting-line check (M0.3, AR1.6, g80x32x32, ep4)

Reference: LL imprint CL 0.0545 (AR1.6) / 0.2178 (AR6); LL CL3 from the same
formula = 0.0637 (symmetric) / 0.2698 (dfvlr_r4 section, camber included).

| mode | section | Gamma mid (LL) | Gamma final | CL | verdict |
|---|---|---|---|---|---|
| 1 | symmetric | 0.0407 | **-0.0092** | **-0.0080** | collapse |
| 1 | dfvlr (camber) | 0.1722 | -0.0090 | -0.008 | collapse |
| 3 | symmetric | 0.0407 | -0.0092 | -0.0080 | collapse |
| 3 | dfvlr (camber) | 0.1722 | +0.0974 (rising) | 0.1489 | wrong magnitude |
| 2 | symmetric | 0.0407 | -0.268 → -0.449 → -0.679 | -0.26 → -0.41 → -0.69 | runaway |

The Kutta rows converge (Kmax 7.5e-7 … 1.6e-6), so these are **converged roots,
not convergence failures**.

## 4. Why — the Γ-scan (mode 3, symmetric, φ converged at each fixed Γ)

`tests3d/tsd_kutta_scan.py` maps the closure function the coupled solve looks
for (`tests3d/a9_gscan_M0.300_m3.json`):

| Γ/Γ_LL | Γ mid | CL | CL from 2∫γdy/S | dCp(TE) mid |
|---|---|---|---|---|
| -1.0 | -0.0407 | -0.0461 | -0.0637 | -0.0179 |
| -0.5 | -0.0203 | -0.0197 | -0.0319 | -0.0059 |
| 0.0 | 0 | +0.0064 | 0 | +0.0052 |
| +0.5 | +0.0203 | +0.0325 | +0.0319 | +0.0166 |
| +1.0 | +0.0407 | **+0.0587** | +0.0637 | +0.0282 |

Two things follow.

* The sheet's *physics* is right: at Γ = Γ_LL it reproduces the lifting-line CL
  (0.0587 vs 0.0545 reference / 0.0637 LL).
* The *closure* is wrong. `dCp(TE)` is monotone in Γ (slope +0.57 per unit Γ)
  and crosses zero at Γ ≈ **-0.009**, which is exactly where every coupled solve
  landed. Cause: `kutta_residual(upstream=2)` samples three chord cells upstream
  of the TE, x ≈ 0.87, where a Kutta-satisfying loading is still ~40% of its peak
  (`A(0.87)=0.42`). Requiring dCp = 0 there is not the Kutta condition — it is
  "zero the aft loading", and the only Γ that does that is ≈ 0 (here slightly
  negative, to cancel the small transpiration TE loading of +0.0052). Moving the
  probe to `upstream=0` (x≈0.98, `A=0.16`) only rescales the bias; it cannot
  remove it, because a prescribed shape that vanishes at x=1 has a positive
  loading at *every* sample inside the chord. This is also why the older free-k
  closure found its root at k≈0.5 instead of 1 (`tests3d/diag_kroot_out.txt`):
  the same row removed half the sheet.

## 5. Why — the other three failure modes

1. **Nothing breaks the symmetry except camber.** With thickness-only
   transpiration (mode 1) the discrete problem is symmetric, so `dCp(TE)=0` holds
   at Γ=0 for any Kutta-satisfying shape: the row has a root at zero with no
   information to select the physical amplitude. Camber breaks it, but then the
   row picks its biased amplitude (Γ=0.097 vs LL 0.172).
2. **The row is span-diagonal.** The selected Γ is flat to ±3% root-to-tip
   (0.974/1.000/0.978 of the mid value vs an elliptic 1.000/0.300). The jump
   source is local to the cut, so station j's TE loading is set by gamma_j: the
   closure picks a *local sheet amplitude*, not a circulation distribution.
   Elliptic loading cannot emerge.
3. **A free wake alone cannot lift.** The wake cut is the trailed half of a
   Helmholtz pair; this solver's own decomposition measures body-only CL +0.267,
   wake-only -0.196. Mode 2 therefore unloads instead of lifting and the TE row
   finds a mirror root in reversed flow (Γ -0.27 → -0.45 → -0.68, CL -0.26 →
   -0.41 → -0.69), which then diverges.
4. **The farfield is not the culprit.** `WING3D_FARFIELD` off vs wake at mode 0
   is bit-identical (CL 0.06515 both), because `run_sup` zeroes every non-inflow
   boundary row: the farfield array only seeds the initial guess.


## 5. Ladder B — transonic CL vs Mach (formulation as specified)

`WING3D_KUTTA=1 WING3D_ENTROPY=1`, dfvlr_r4, span 1.6, g80x32x32, ep8
(`tests3d/a9_dfvlr_r4_M0.*_k1_entropy_ep8.json`, `tmp_a9_B1.log`):

| M | anchor (fine 125k) | CL | err % | Γ_mid (LL sew = 0.2366) | converged? |
|---|---|---|---|---|---|
| 0.850 | 0.3210697 | 0.02852 | **-91.1** | 0.0121 | yes (\|F\|max 3.4e-4) |
| 0.925 | 0.3629209 | 0.04460 | **-87.7** | 0.0214 | no (\|F\|max 2.9) |
| 1.000 | 0.3037588 | 0.05422 | **-82.2** | 0.0291 | no (\|F\|max 4-6, nsup 67k) |

CL does *not* lose lift at transonic Mach here — because there is no circulation
left to lose. What is left is the residual thickness/camber transpiration lift,
and it drifts *up* with M (0.029 → 0.045 → 0.054). The signature the task asks
for (a fall past M_crit) is absent, and the honest reading is the intended one:
**a monotone rise with a 82-91% level error means the circulation is still not
free** — it is now not merely pinned but destroyed by the closure.

## 6. Ladder C — BL interplay (M0.85, ep8, entropy on)

| BL | CL | Δ vs BL0 | Γ_mid final | δ*/c max | sep |
|---|---|---|---|---|---|
| 0 | 0.02852 | — | 0.0121 | — | — |
| 1 | 0.02964 | **+3.9% (+0.0011 CL)** | 0.0132 | 1.1e-3 | 0 |

Not neutral any more (A8 measured BL0 -2.1% / BL1 -2.2%, i.e. 0.1 pp), but the
+0.0011 CL is a change in the *residual* lift, not viscosity acting on a
physical circulation — the circulation it could act on has already collapsed to
5% of the lifting-line value. Viscosity is live; there is simply no lift left
for it to reduce. δ*/c max = 1.1e-3 is small because the BL is being marched on
a surface flow with almost no suction peak (Cfmax 0.39, nsup 182).

## 7. Ladder D — legacy free_k, and why it agrees

`WING3D_KUTTA=0` (imprint sheet at k=1) + `free_k=True`, M0.85, ep4, entropy:

| ep | k_mid | CL |
|---|---|---|
| 0 | 0.095 | 0.0381 |
| 1 | 0.128 | 0.0384 |
| 3 | 0.128 | 0.0388 |

free_k still collapses, now to **k ≈ 0.13** (A8 recorded ~0.02) with CL 0.0388
against 0.3143 at k=1. It is the same failure mode as the new closure and the
same cause: the row `kutta_residual(upstream=2)` samples the loading at x≈0.87
where a Kutta-satisfying loading is still 40% of peak, so its root is "no aft
loading" (k→0), not "smooth TE flow" (k=1). free_k's k≈0.13 and mode 1's
Γ/Γ_LL≈0.05 are two measurements of one structural defect.

## 8. Verdict

* The plumbing works and is committable: the circulation *is* an unknown, the
  states/rows/re-seeding/farfield/weights are all in place, and every solve
  converges (Kutta rows to 1e-6…1e-7) instead of blowing up.
* The **closure** does not work, and it cannot work with a prescribed loading
  shape: the Kutta condition constrains the loading at the TE, a Kutta-shaped
  sheet satisfies it for *every* amplitude, and the only usable equation left
  (`dCp` at 87% chord, or the free wake alone) has its root at Γ≈0 or on a
  reversed-lift mirror branch. Measured span-flat Γ confirms the row is
  span-diagonal: it sets a local amplitude, not a circulation distribution.
* **Is the circulation genuinely free? No.** It became an unknown but not a
  free parameter: the closure removes it rather than selecting it.
* **Does CL now lose lift at transonic Mach? No** — CL rises 0.029 → 0.054 while
  sitting 82-91% below the anchors, which is the signature of *no* circulation,
  not of a free one.

What would actually close it (the item this work isolates): a phi-**dependent**
tangency row, i.e. the wing BC must be an equation in phi (the sheet's induced
normal velocity balancing `alpha + dy_c/dx`) rather than a fixed source, plus
the Morino wake condition `Γ = Δφ(TE)`. That is the VLM+TSD coupling already
described in `wing3d/tsd.py` lines 17-21; nothing in the present architecture
can substitute for it. A second, cheaper improvement measured here: move the
Kutta probe to the TE cell (`upstream=0`) — it shrinks the bias but cannot
remove it, so it is not a fix on its own.

Files: `wing3d/tsd.py` (`flat_plate_jump`, `gamma_b`, `kutta_mode`,
`kutta_thick_only`), `tests3d/tsd_m_sup.py` (state vector, rows, LL seeding, φ
retention), `tests3d/tsd_coupled.py` (`lifting_line_gamma`), tools
`tests3d/tsd_kutta_weight.py`, `tests3d/tsd_kutta_scan.py`,
`tests3d/tsd_a9_ladder.py`, `tests3d/tsd_a9_ff.py`. Results
`tests3d/a9_*.json`, logs `tmp_a9_*.log`. Nothing outside my four files + these
new test tools was touched; `cases/openfoam/**` and
`tests3d/validation/*.json` untouched.
