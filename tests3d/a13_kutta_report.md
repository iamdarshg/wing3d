# A13 — hybrid row: `w_field(phi) − (L·E − D)·mu`, i.e. the field's own
# downwash with the estimator's bias removed, so the mu block depends on phi

Status: **the mu block is now genuinely phi-coupled** — `d(row)/d(phi)` is
25.5 vs 3.2 for `d(row)/d(mu)`, where mode 6's was *identically zero* — and
the construction is exact in the sense that matters: fed the analytic sheet
that `E` models, the mode-7 rows reproduce the mode-6 root to **2e-16**.
`gamma_mid` is no longer bit-identical across Mach and CL does move. But the
two acceptance levels are **not** met and the reason is measured, not
speculative: **test A lands −26.1%** (0.0340 vs A12's 0.0460, uniform
downwash offset ≈0.31α) and **test B's required monotone 15–25% transonic
fall does not appear** (CL −2.6%/+1.7% at M0.925/1.00, non-monotone;
`gamma_mid` 0.02975/0.02831/0.02960). Section 6 decomposes both into a
non-circulatory mean-plane component of the field that `B`, built only from
free-space quantities, cannot cancel. No tuning was done.

Gates: `WING3D_KUTTA=7` (new; modes 0–6 untouched), everything else as A12:
`WING3D_KUTTA_TE`, `WING3D_KUTTA_W` (8.9e3 inherited, see §1.4),
`WING3D_KUTTA_KERNEL_*`, `WING3D_ENTROPY`, `WING3D_BL`.

## 1. The construction

### 1.1 What is `w_field_equivalent(mu)`

The task's fix is `row = w_field(phi) − w_field_eq(mu) + D.mu − g`. The
executable content is in defining *which* operator `w_field_eq` is.
`w_plane()` (A10) is not a generic "downwash of a sheet": it is a
**specific estimator geometry**,

    wu = (phi[k0+1] − phi[k0]) / (z[k0+1] − z[k0])
    wl = (phi[k0−1] − phi[k0−2]) / (z[k0−1] − z[k0−2])
    wav = (wu + wl)/2

A one-sided difference of a smooth potential *is* the normal velocity at the
cell centre, so the estimator reads `∂φ/∂z` at **two z**, and on g80 those
are **not symmetric**: `+0.00736c` above the sheet and `−0.02456c` below
(`tsd._stencil_heights`). `w_field_eq` must be built at exactly those
heights, with the same two-sided average, or the height bias (which is most
of the A11 bias) does not cancel. So `E = dlm_field_equivalent()`
(`tsd.py`) evaluates the *same horseshoe lattice* with
`ckernels.loop_vel_block` at `(x_c[iu], y_c[jy], +0.00736c)` and
`(…, −0.02456c)`, averages, and maps segment→cell→mu with the same `S`,`C`
matrices mode 6 uses. Rows are the field's own wing nodes, because that is
where `w_plane` lives.

The lattice collocates at the 3/4 chord between stations (A12, validated
there — it must not move), so `L = dlm_obs_matrix()` is the (nc·ng,
n_node·ng) linear map node-lattice → control points (bilinear, linear
extrapolation at the LE cell). `L` is applied **identically** to `w_field`
and to `E_node`, so it cancels in `B` by construction — only the
*difference* has to be smooth, not either side alone.

**`B = L·E_node − D`** is a fixed matrix: grid- and kernel-determined, no
solution enters (verified: `B` is identical to 4 significant digits on the
AR1.6 and span-40 lattices, i.e. it is local sampling bias, not lattice
content). The row is

    row = w_field(phi) − B·mu − (dyc/dx − α)      (row nc−1: Kutta, as modes 5/6)

linear in `mu`, dependent on `phi` only through `w_field`. Documented
choice: this is the task's **fixed-B** form, and the reason is in §1.2 — in
the ideal limit the row collapses *exactly* onto mode 6, so nothing about
A12's validated chordwise shape can be lost; every deviation is a property
of the real field.

### 1.2 The bias-cancellation check (no TSD solve; `tsd_a13_probe.py`)

Exact 2D flat-plate doublet distribution on the lattice,
`mu = π·α·F(x)`, read back through each operator (span 40, mid-span, in
units of α):

| x_cp | .032 | .096 | .168 | .243 | .321 | .400 | .481 | .563 | .646 | .729 | .813 | .896 | .978 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **D·mu/α** (A12 kernel) | −1.13 | −0.94 | −1.02 | −1.03 | −1.03 | −1.03 | −1.03 | −1.03 | −1.03 | −1.03 | −1.03 | −1.03 | −1.02 |
| **E·mu/α** (estimator) | +1.43 | +0.82 | +0.30 | +0.01 | −0.17 | −0.30 | −0.40 | −0.49 | −0.57 | −0.65 | −0.73 | −0.81 | −0.89 |
| **B·mu/α** (what mode 7 subtracts) | +2.56 | +1.77 | +1.33 | +1.04 | +0.86 | +0.73 | +0.63 | +0.54 | +0.45 | +0.38 | +0.30 | +0.22 | +0.13 |

The estimator's chordwise bias against the kernel is **2.4α of tilt over the
chord** (the A11 3.7x effect, now with a mechanism: the two unequal stencil
heights), and `B` is exactly that number.

**Exactness:** with the field replaced by its analytic equivalent
(`w_plane := L·E·mu`), the mode-7 Newton solve returns
`max|mu_hybrid − mu_kernel| = 2.6e−16` (span 40) and `2.1e−17` (AR1.6),
residual 8e−17. The assembly (E, L, B, the 3/4-chord map, the signs) is
therefore *provably* the mode-6 equation plus nothing.

**Liveness:** `d(row)/d(phi)` = **25.46** (central difference on the four
stencil nodes of the mid control point), `d(row)/d(mu)` = **3.17**,
ratio **8.0**; the same measurement on mode 6's rows is **0.000e+00** —
exactly zero, by construction, which is A12's block-triangular defect.

**Regression:** mode-6 rows untouched — span-40 `CL_circ = 0.11685 =
0.852·2πα` (A12: 0.852), AR1.6 `CL/α = 2.22/rad` (A12 recorded 2.24,
Helmbold 2.20; the 1% is the standalone-vs-ladder grid, `ny=64` here).

### 1.3 Live coupling in the coupled Jacobian

`B` is not small: `max|B| = 4.19`, `max row-sum|B| = 8.64` vs `max|D| =
13.36` (AR1.6, g80). At convergence the tangency rows reach 1.5e−7 and the
Kutta row 1.7e−7 with `nsup = 0`, so the coupled solve is consistent, not
merely nominally coupled.

### 1.4 Row weight

`kutta_w = 8.9e3` is inherited from mode 6/4 (same velocity units, same
`A_scale` construction). Unlike mode 6 this weight *can* now move the root
(the mu block depends on phi); it is declared here as conditioning
inheritance, **not** tuned, and the answer's sensitivity to it is not
measured in this round.

## 2. Test A — M0.3 symmetric, AR1.6 (dfvlr thickness, profile=None)

| quantity | mode 6 (A12) | **mode 7** |
|---|---|---|
| **CL** | **0.04601** | **0.03400 (−26.1%)** |
| 2Σ(Γ dy)/S | 0.0489 | 0.0339 |
| Γ_mid | 0.0298 | 0.0208 |
| tip/mid | 0.64 | **0.61** (LL 0.56) |
| tangency residual | 6.6e−9 | 1.5e−7 |
| Kutta row | 1.2e−8 | 1.7e−7 |
| \|F\|max / nsup | 1.2e−4 / 0 | 2.5e−3 / 0 |
| d(row)/d(phi) | **0** (structural) | **25.5** (dμ 3.2, ratio 8.0) |

Epochs: 0.0272 / 0.0302 / 0.0305 / 0.0324 / 0.0325 / 0.0339 / 0.0340 /
0.0340 — converged (8 epochs, 718 s, `a13_…_ep8_noprof.json`).

Chordwise loading at mid-span (x from .043 to .999):

    gamma .053 .036 .028 .022 .017 .013 .011 .008 .007 .005 .004 0.000 0.000

LE peak, monotone fall, `gamma_TE = 0.000` (Kutta), i.e. **A12's shape is
preserved exactly** and the elliptic span survives (tip/mid 0.61 vs 0.64).

**Verdict A: shape/span/Kutta PASS, coupling liveness PASS, level FAIL
(−26.1%, criterion ±3%).**

## 3. Test B — transonic dfvlr, ENTROPY=1, ep2

| M | CL | err vs anchor | Γ_mid | tip/mid | nsup | \|F\|max |
|---|---|---|---|---|---|---|
| 0.850 | 0.06378 | −76.5% | 0.02975 | 0.81 | 87 | 4.2e−2 |
| 0.925 | 0.06210 | −76.4% | 0.02831 | 0.92 | 279 | 8.0e−1 |
| 1.000 | 0.06318 | −72.2% | 0.02960 | 0.82 | 66 871 | 4.9 |

Anchors 0.27186/0.26253/0.22694 (older 0.2719/0.2625/0.2269 unchanged).

* **`Γ_mid` now VARIES with Mach** (−4.8% then +4.6% from M0.85) where mode 6
  was `0.0618001 / 0.0617992 / 0.0617953` — bit-identical to 5 digits. The
  defect A12 exposed is structurally gone.
* **CL still does not fall:** −2.6% then +1.7%. Mode 6 was +1.2%/+5.9%. The
  required monotone 15–25% collapse is **not reproduced**; the response that
  exists is 5x too small and non-monotone.
* M1.00 is not a converged field (66 871 of 80 896 nodes supersonic, |F| 4.9,
  as in A12): no Mach statement at M1.0 is possible until the supersonic
  branch converges, independent of the closure.

**Verdict B: FAIL on the required signature; PASS on "the circulation now
sees the Mach number".**

## 4. Test C — BL on top of the hybrid row (M0.85, ep2, `WING3D_BL=1`)

CL **0.06380** vs 0.06378 (+0.03%), Γ_mid 0.02972 vs 0.02975, δ*/c max
1.16e−3, w correction 1.4e−4, no separation. **Still neutral**, and the
mechanism is now identifiable in the code: the displacement enters the
transpiration **antisymmetrically** (`wu += ddu`, `wl −= ddl`,
`_apply_wing_bc`), the same way thickness does, so it cancels in
`wav = (wu+wl)/2` — the quantity the tangency row reads. Its size (1.2e−3 vs
α = 2.2e−2, 5%) is also 6x below the ~0.3α offsets of §6. The BL is in the
field but not in the row that sets the circulation.

## 5. Test D — cambered dfvlr, M0.3

Cambered finite-span estimate (A12): `cl2 = 0.519` → `CL3 = 0.191`.

| run | CL | vs 0.191 | cambered/symmetric |
|---|---|---|---|
| mode 6 (A12) | 0.1065 | 56% | 2.31x |
| **mode 7** | **0.0632** | **33%** | **1.86x** |

Aft-loading shape (mid-span loading, x .043→.999):

    +.086 +.044 +.016 −.007 −.021 −.017 −.018 −.011 +.036 +.123 +.192 0.000 0.000

The aft peak at x=0.83 **persists and sharpens** (A12: +0.181 at x=0.83 with
a positive mid-chord); the mid-chord loading is now *negative*. The
camber-driven deficit is the aft mean-line demand interacting with the
mean-plane estimator (§6), not a Mach effect — the same run at M0.85 gives
0.0638 (+0.9%), so the level is Mach-flat in the cambered family.

## 6. Why the levels moved: the row decomposition (measured)

The stored record contains the row's three pieces at the mid span; the
decisive quantity is `resid_phys = w_field − w_eq`, the part of the field's
mean-plane velocity that the analytic sheet of the same mu does **not**
explain (units of α):

| case | x .043 → .833 |
|---|---|
| M0.3 symmetric | −0.36 −0.39 −0.36 −0.34 −0.32 −0.31 −0.31 −0.32 −0.32 −0.32 −0.31 |
| M0.3 cambered | −0.88 −0.53 −0.31 −0.12 −0.01 −0.05 −0.12 −0.29 −0.81 −1.65 −2.22 |
| M0.85 | −1.56 −0.62 −0.22 +0.02 +0.11 −0.02 −0.17 −0.49 −1.26 −2.31 −2.81 |

1. **The estimator's chordwise tilt IS cancelled.** Individually,
   `w_field/α` and `w_eq/α` each swing ≈0.8α along the chord; their
   difference at M0.3 symmetric is flat to **0.14α** (spread/chord vs 2.4α
   for the raw bias of §1.2). That is the mechanism the task asked for, and
   it works.
2. **What survives is a uniform scale offset.** At M0.3 symmetric the
   residual is a near-constant ≈−0.31α, which is algebraically an extra
   downwash: mu is scaled by 1/(1+0.31) = 0.764, predicting
   `0.0460·0.764 = 0.0351` — measured 0.0340. **The entire −26% of test A is
   this one number**: the discrete field's mean-plane downwash response to
   the same doublet sheet is ≈1.36x the free-space kernel's. `B` cannot see
   it: `B = L·E − D` is built exclusively from free-space quantities (the
   horseshoe kernel and its own sampling), so a *normalisation* error of the
   discrete operator's Green's function — box farfield, 1st-order
   x-upwinding, `dz_core` anisotropy — passes straight through. (Same family
   as the 2D lattice's own 0.852-of-2πα deficit, which is discretisation,
   not physics.)
3. **Camber/transonic add a chordwise ramp on top.** `resid_phys` goes from
   flat (−0.31±0.07) to a −0.9…−2.2α ramp with camber, at **M0.3**, so the
   ramp is not Mach-driven; M0.85 makes it ~40% larger at the LE. This is
   the non-circulatory (thickness + camber + asymmetric supersonic
   response) mean-plane velocity that the `(wu+wl)/2` cancellation removes
   only in the linear limit. It pulls the mid-chord loading negative and the
   integral CL down, which is why B's non-monotone ±5% is all the Mach
   signal that can reach CL: the physical signal is buried under a
   non-circulatory term 5–10x larger that also changes with Mach.

## 7. Verdict

* **The hybrid row is implemented correctly.** Bias operator `B` fixed and
  grid-determined, ideal-field limit reduces to mode 6 to 2e-16, `d(row)/d(φ)`
  25.5 vs 3.2 for `dμ` (mode 6: exactly 0), converged solves with tangency
  1e−7 and the Kutta row intact.
* **Circulation is now physically shaped AND Mach-coupled — in structure.**
  A12's validated shape (LE peak, monotone fall, `gamma_TE=0`, tip/mid 0.61)
  is retained, and `Γ_mid` genuinely responds to Mach.
* **But the coupling is not physical in magnitude.** Test A is −26.1%
  (−0.31α uniform → the discrete Green's function's 1.36x normalisation),
  and test B's CL is flat and non-monotone (±2–5%) because `w_field` is
  dominated by a non-circulatory mean-plane component (uniform + camber
  ramp, 0.9–2.8α at M0.85) that `B` — free-space by construction — cannot
  remove. Test C remains neutral by the same route (antisymmetric entry,
  cancels in `wav`). Test D persists at 33% of the finite-span estimate with
  the aft-loading collapse sharpened.
* **Nothing was tuned.** The level is where the physics the row now reads
  puts it.

**Next step that would settle it** (not done here): make the reference
*field* part of the bias operator — subtract a fixed **thickness-only**
field's mean-plane `w` (one extra solve with the circulation sheet
disabled) as `B' = B + w_thick`, which removes the non-circulatory response
*including* its compressible part, and then test whether the uniform 1.36x
offset separately shrinks under a `dz_core`/`nx` refinement (discretisation)
or not (physics). Until then the residual freedom is a genuine
model/reference difference of the discrete mean-plane estimator, measured —
not a bug in the row.

Files: `wing3d/tsd.py` (`_stencil_heights`, `dlm_field_equivalent`,
`dlm_obs_matrix`, `dlm_hybrid_bias`, `hybrid_tangency_residual`,
`hybrid_coupling_report`, mode-7 gate), `tests3d/tsd_m_sup.py` (mode-7 rows,
`tang_rows`, row decomposition in the record), `tests3d/tsd_a13_probe.py`
(the bias exhibit, the exactness and liveness checks, the mode-6
regression), `tests3d/tsd_a13_ladder.py` (tests A–D, recorder `one7`).
Results `tests3d/a13_*.json`; logs `tmp_a13_*.log`.
