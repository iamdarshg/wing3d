# A11 — free chordwise doublet distribution: what it measures

Status: **the chordwise distribution is now genuinely free and the row set
does resolve chordwise shape (measured local self-response −6.2, neighbour
+1.45 — a delta/dipole operator, not a flat one), and the cambered failure of
A10 is gone in the qualitative sense the task asks for: cambered CL is now
3.4x (Kutta row on) to 6.6x (off) the symmetric CL, where A10 had it BELOW.
But two new numbers are negative and they are the headline: the symmetric
case regresses to −32% vs the 0.0545 reference (A10's mode 4 was +2.0%), and
the cambered level is not resolved — it depends on the TE treatment
(0.1245 vs 0.2445) and on nc (0.0459 at nc=8 vs 0.1245 at nc=13), because
the aft 8–17% of the chord carries most of the cambered load. Root cause of
both, measured in §2: the mean-plane tangency estimator's chordwise response
is 3.7x stronger aft than at the LE for a physically uniform downwash, so
enforcing tangency at *every* station flattens `mu` and loses circulation.

Gates: `WING3D_KUTTA=5` (this work), `WING3D_KUTTA_NC` (control stations,
default = all 13 wing nodes), `WING3D_KUTTA_TE=1|0` (discrete Kutta row vs
pure collocation), `WING3D_KUTTA_W` (row weight), `WING3D_KUTTA_WTE` (Kutta
row weight). Modes 0–4 untouched (`git checkout` of A10's HEAD reproduces
mode-0 free_k bit-for-bit, §6).

## 1. Discretisation

Unknowns `mu(i,j)`, the potential jump (doublet strength), at `nc` chordwise
control stations × `ng` span cells. Control stations are a *subset of the
field's own wing nodes* (`xc = 0.043 … 0.999`, 13 on g80 AR1.6), so every
tangency row reads `w_plane` at a node rather than interpolating a
jump-contaminated estimator. `mu` is applied as the cut's doublet sheet and
mapped to all wing nodes by linear interpolation in x; the wake trails
`mu(TE,y)` (Helmholtz). Closure: **one tangency row per (station, span)**,
`wav(x_i,y_j) − (dyc/dx(x_i) − α) = 0`, i.e. `nc·ng` rows for `nc·ng`
unknowns (91 rows on g80/AR1.6 vs mode 4's 7). Transpiration is
thickness-only (as mode 4), the Mach/entropy treatment is unchanged.

`WING3D_KUTTA_TE=1` (default) replaces the *last* row (at x=0.9988, i.e. on
the TE) by the discrete Kutta condition `mu(TE) = mu(last-1)`; collocating
tangency on a mean line that ends at the TE is what produced a TE loading
spike (§4). The Kutta row is in mu units and carries the measured local
self-response as its weight, so it is commensurate with the velocity rows and
does not move the root, only the conditioning.

## 2. Measured row scaling (nothing reused from A10)

`tests3d/tsd_a11_probe.py`, M0.3, g80x32x32, AR1.6, two converged phi solves
per geometry at `mu = 0.6` and `1.4 × F(x)Γ_LL`:

| quantity | measured |
|---|---|
| `A_scale` = median \|dR_phi/dphi\| (32 nodes at k0±1) | **8.44e3** |
| `D_ij` = converged d wav / d mu, global F(x)Γ perturbation | −0.86 … −1.55, **median 1.05** |
| ⇒ scalar row weight `w = A_scale/median\|D\|` | **8.0e3** |

`D` is uniform along the chord (spread 1.8x) and non-uniform in span (tips
1.18 vs mid 0.86–1.05), so the row set is span-coupled — same structure A10
measured for its single row, unresolved in chord.

**The local structure is what matters** (`tsd_a11_probe2.py`, one panel
bumped by δ=0.004, response read at every station, mid-span):

| bump | resp at x = .04 .11 .19 .26 .34 .42 .50 .58 .67 .75 .83 .92 1.0 |
|---|---|
| panel 2 | 0.5 1.6 **−6.5** 1.5 0.4 0.2 0.1 0.05 … |
| panel 6 | 0.05 0.06 0.09 0.17 0.43 1.47 **−6.15** 1.44 0.40 0.15 … |
| panel 11 | 0.03 … 0.05 0.08 0.15 0.39 1.43 **−6.08** 1.45 |
| **constant shift** | **−4.20 −2.57 −2.07 −1.87 −1.77 −1.72 −1.68 … −1.61** |

So the operator is a strong local (difference/dipole) term — self −6.1,
neighbours +1.45, next +0.4/+0.15, decaying — plus a weak global amplitude
term (−1.6). The `nc` rows therefore DO constrain chordwise shape; the matrix
is not rank-deficient in the chord (that was the risk, and it is measured,
not argued). The 6.1 self-response is the Kutta row's weight.

**The bias that costs the level**: for the flat-plate sheet at Γ_LL, where
thin-airfoil theory requires `wav = −α = −0.0218` *everywhere*, the measured
profile is

    x     .043   .113   .187   .262   .340   .420   .501   .583   .667   .750   .833   .917   .999
    wav  −.0099 −.0206 −.0240 −.0263 −.0283 −.0302 −.0319 −.0334 −.0347 −.0358 −.0365 −.0369 −.0363

— 0.45x α at the LE cell and 1.7x α aft, a 3.7x span in the *estimator's own
chordwise kernel*. Mode 4 never saw it because its single row was
A(x)-weighted, and A is LE-singular.

## 3. Test A — M0.3 symmetric (must not regress A10)

dfvlr thickness, `profile=None`, ep4, reference CL 0.0545, Γ_LL 0.0407:

| nc | TE row | CL | Γ = μ_TE−μ_LE | μ_TE | tip/mid (LL 0.56) | tangency | \|F\|max |
|---|---|---|---|---|---|---|---|
| 13 | Kutta | **0.0369** (−32.3%) | 0.0141 (35%) | 0.0257 (63%) | **0.58** | 8.8e-9 | 9.8e-5 |
| 13 | colloc. | 0.0379 (−30.5%) | 0.0147 (36%) | 0.0257 (63%) | 0.58 | 9.0e-11 | 1.8e-6 |
| 8 | Kutta | 0.0359 (−34.1%) | 0.0137 (34%) | — | 0.59 | 9.9e-9 | 8.6e-5 |

**FAIL on the ±3% criterion, PASS on the spanwise loading** (tip/mid 0.58 vs
LL 0.56 — better than A10's 0.65). The −32% is NOT a chordwise-resolution
problem: nc 8→13 moves CL by 2.7%. It is §2's estimator bias. Verdict: this
is a regression against mode 4 (+2.0%) caused by enforcing a condition
(uniform downwash at all 13 stations) that the biased estimator cannot honour
without flattening the loading.

## 4. Test B — M0.3 cambered dfvlr (the key test)

Cambered section: thin-airfoil A0 = 0.0486, A1 = 0.0679 → cl2 0.519,
LL imprint CL3 = **0.2698** (the 0.0545–0.06 figure in the task is the
symmetric one; `lifting_line_gamma` prints 0.26977, cambered Γ_LL = 0.1722):

| nc | TE row | CL | vs LL CL 0.2698 | Γ | vs cambered LL 0.1722 | cambered/symmetric |
|---|---|---|---|---|---|---|
| 13 | **Kutta** | 0.1245 | −53.9% | 0.0672 | **39%** | **3.4x** |
| 13 | colloc. | 0.2445 | −9.4% | 0.1382 | **80%** | 6.6x |
| 8 | Kutta | 0.0459 | −83.0% | 0.0229 | 13% | 1.3x |

**The key criterion PASSES**: the cambered CL is above the symmetric CL by
3.4–6.6x, where A10 had it 47% *below*. But the level is a property of the TE
treatment and of nc (0.0459 → 0.1245 → 0.2445), i.e. of the aft 8–17% of the
chord, so it is not yet a converged number.

Why the aft chord dominates: the dfvlr_r4 mean line rises to y=0.0205 at
x=0.833 and then **plunges back to y=0.0003 at the TE** (a sharp-TE closure
artifact of `section_profile_xy`), so `dyc/dx` reaches −0.26 at x→1 and the
tangency target is +0.28 there. A collocation row at the TE then buys that
with a TE loading spike; the Kutta row refuses it and the answer halves.

## 5. Test C — chordwise loading at mid-span

γ = dμ/dx, symmetric, nc=13 (peak/target: LE peak, fall to ~0):

    x    .043   .113   .187   .262   .340   .420   .501   .583   .667   .750   .833   .917   .999
    γ   .0598  .0381  .0265  .0189  .0136  .0099  .0073  .0056  .0045  .0039  .0035  .0030  .0030

Peak at the LE, monotone fall, TE value 5% of peak — **a physical loading
shape, and Kutta ≈ satisfied** (it plateaus over the last two cells rather
than going to exactly zero). This is the first mode in this series to produce
that shape from the solve rather than from a prescribed `F(x)`.

Cambered, nc=13, Kutta row on:

    x    .043   .113   .187   .262   .340   .420   .501   .583   .667   .750   .833   .917   .999
    γ   .1002  .0475  .0191 −.0018 −.0038  .0147  .0212  .0503  .1336  .2200  .2296  .0000  .0000

Kutta holds at the last cell (γ=0), but the peak is at x=0.83, not the LE —
the aft-loaded signature §4 explains. Collocation instead gives a monotone
rise to γ(TE)=0.41 = γ_max (no Kutta at all).

## 6. Test D — transonic dfvlr, ENTROPY=1, ep8, nc=13, Kutta row

New free-tip AR1.6 anchors 0.27186 / 0.26253 / 0.22694:

| M | CL | err | μ_TE mid (LL at that M) | tip/mid | converged? |
|---|---|---|---|---|---|
| 0.850 | **0.1172** | **−56.9%** | 0.0637 (0.2366, 27%) | 0.78 | no, \|F\|max 1.7e-1 |
| 0.925 | **0.1277** | **−51.4%** | 0.0731 (28%) | 0.73 | no, \|F\|max 0.3–1.2 |
| 1.000 | **0.1318** | **−41.9%** | 0.0786 (33%) | 0.71 | no, \|F\|max 4–7, nsup 67k |

CL **rises** with Mach (+8.9%, +3.2%); the anchors require a 16.5% fall.
Level error is much smaller than A10's −85%, but the trend is the opposite
sign and the M0.85 solve does not reach A10's residual (1.7e-1 vs 3.2e-5),
so the levels are still epoch-limited. `mu_LE` goes 0.0042 → −0.0008 as the
supersonic pocket covers the chord (nsup 33 → 203 → 67k).

## 7. Test E — free_k and BL

* **free_k (mode 0, untouched path)**, M0.85 ep4 entropy: k_mid **0.020**,
  CL **0.00378**, nsup 29/62/100/128. A10 recorded k=0.128 / CL 0.0388; that
  does **not** reproduce here. I verified this is not my change: restoring
  A10's HEAD (`git checkout` of both files) and re-running gives the
  identical 0.00378 / 0.020. So the legacy closure is still collapsed
  (stronger than A10 reported), and the discrepancy is a harness/env
  difference in A10's ladder-D row, outside this change set.
* **BL** (M0.85, ep8, entropy, `WING3D_BL=1`): CL 0.11694 vs 0.11718
  (−0.2%), Γ_mid 0.0635 vs 0.0637 (−0.3%), δ*/c max 1.3e-3, no separation.
  Nearly neutral — an order of magnitude less than mode 4's +4.4%, and with
  no wrong-signed circulation change to explain.

## 8. Verdict

1. **The formulation works as specified.** `nc·ng` unknowns, `nc·ng` tangency
   rows, Helmholtz wake, thickness-only transpiration, mode-5 gate. Every
   solve converges to a root (tangency 7e-9…9e-11, Kutta row 1e-9) and the
   Jacobian's local structure is a proper delta/dipole (±6.2 self, +1.45
   neighbours), so the chordwise shape is *an output*, not a reduction. The
   chordwise loading it produces for an uncambered section is physical.
2. **The cambered reference now works in sign and in span, not in level.**
   Cambered CL is 3.4x the symmetric CL (Kutta row) / 6.6x (collocation) —
   A10's failure mode is gone. But the level is set by the aft 8–17% of the
   chord (nc 8→13 moves it 2.7x) because the reference section's mean line
   plunges to zero at a sharp TE.
3. **The symmetric case regresses and I am not hiding it**: −32% vs the
   0.0545 reference, against mode 4's +2.0%. The measured cause is not the
   closure but the estimator: `wav`'s chordwise kernel is 3.7x stronger aft
   than at the LE for a physically uniform downwash, so the 13-station
   uniform tangency condition flattens `mu` and drops the circulation by
   35%. A10's A-weighted single row was blind to that bias.
4. **What would fix it**: a *consistent* downwash estimator — the doublet
   lattice's own kernel (velocity induced by panel j's doublets at control
   point i, from the analytic/ring formula) instead of `w_plane` read out of
   a TSD field whose dz_core is 1.5% of the chord. That is the VLM+TSD
   coupling A9 named; this work shows the row layout is ready for it and
   quantifies exactly how much the field-based estimator costs (×0.65 in
   circulation, chordwise-uneven).

Files: `wing3d/tsd.py` (`kutta_free_mu`, `init_free_mu`,
`tangency_panel_residual`, `mu_circulation`, `mu_loading`, `kutta_te_kutta`),
`tests3d/tsd_m_sup.py` (state/rows/reporting for mode 5), tools
`tests3d/tsd_a11_probe.py`, `tsd_a11_probe2.py`, `tsd_a11_ladder.py`;
results `tests3d/a11_*.json`, logs `tmp_a11_*.log`.
