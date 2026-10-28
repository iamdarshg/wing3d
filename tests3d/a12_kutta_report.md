# A12 — analytic doublet influence kernel for the coupled Kutta closure

Status: **the tangency rows are now the classical doublet-lattice equation
`D.mu = dyc/dx - alpha`, with `D` built from the same Biot-Savart kernel the
codebase already uses for the Morino panel solve, and it is VALIDATED: fed
the exact 2D flat-plate doublet distribution it returns `w = -0.94..-1.03
alpha` along the chord (A11's field estimator returned `-0.45..-1.7`), and
it reproduces the classical finite-wing lift slope (2.24/rad for AR1.6 vs
Helmbold/DATCOM 2.20/rad). The symmetric M0.3 level improves from A11's
−32% to **−15.6%** and the chordwise loading is physical, but it does **not**
reach the 0.0545 anchor: −15.6% is reported as a *reference difference*
(section 4 explains why, with the numbers), not tuned away. The cambered
section is physical in shape and sign (2.3x the symmetric CL) but its level
is set by the TE treatment because this reference mean line plunges to zero
at the trailing edge. And mode 6 has a structural cost that the task's test
D exposes: because the mu block is linear and phi-independent, **the
circulation has no Mach dependence at all** — CL is Mach-flat 0.105/0.100/
0.106 at M0.85/0.925/1.00.

Gates: `WING3D_KUTTA=6` (new; mode 5 and 0–4 untouched),
`WING3D_KUTTA_TE=1|0`, `WING3D_KUTTA_KERNEL_WKP/WLEN` (wake panels/length),
`WING3D_KUTTA_KERNEL_PG=1` (Prandtl-Glauert-stretched kernel, **ablation
only**), `WING3D_KUTTA_NC` as before.

## 1. The kernel, and where it came from

A constant doublet panel of strength `mu` is equivalent to a vortex ring of
circulation `mu` around its edges (Katz & Plotkin 10.4; the identity the
Albano–Rodden doublet lattice is built on). For a planar sheet, the velocity
a vortex ring induces at a point **in its own plane** is purely normal to
that plane, so the classical doublet-lattice influence coefficient *is* the
z-component of the ring kernel — no new singular integral is needed.

Reuse: `wing3d/ckernels.py:loop_vel_block` (C, with the pure-numpy fallback
`wing3d/solver.py:_loop_vel_rows`), i.e. exactly the Biot-Savart source
`solver.assemble` uses for the Morino doublet system. `tsd.dlm_kernel()`
assembles horseshoes from it.

Formulation (all of `tsd._dlm_panels`, `tsd.dlm_kernel`,
`tsd.dlm_tangency_residual`):

* Cells: chordwise between **stations** (`[x_{i-1}, x_i]`, `x_{-1}` = LE),
  spanwise centred on a station (edges at the midpoints, tips at ±span/2).
  The unknown loading of cell `i` is `Gamma_i = mu_i - mu_{i-1}` with
  `mu(LE) = 0` **by construction**, so `mu(TE) = sum Gamma` is the section
  circulation.
* Each cell is one horseshoe: bound filament at its 1/4 chord running +y,
  legs trailing downstream (4 panels, growth 1.6, total 1.6 span). The
  trailed wake is therefore **inside** the kernel, and the closed-loop
  orientation makes adjacent span cells leave `Gamma_{l-1}-Gamma_l` on the
  shared leg.
* Control points: 3/4 chord of each cell, mid-span of each span cell
  (classical 1/4–3/4). Row `nc-1` is the discrete Kutta `mu(TE) =
  mu(last-1)` (A11's better variant). `D` is `(nc*ng)^2`, i-major to match
  `mu_b.ravel()`; column map `Gamma = (I - shift) mu`.
* Coupling: mu **stays a Newton unknown** (same state layout as mode 5).
  The phi rows depend on mu (the cut jump and the trailed wake); the mu
  rows do not depend on phi, so the coupled Jacobian is *block-triangular*
  and the Krylov solve converges the mu block in one step (measured: the
  tangency residual is 8e-9..1e-12 from epoch 0 on, and CL moves 0.0461 →
  0.0460 over the epochs). That is the honest description of the loop the
  task allowed: [phi with mu] → [linear mu solve] is exact, not iterated.

Two formulations were measured and rejected, and the numbers are the
reason the code looks the way it does:

| form | flat-plate check `w/alpha` | 2D CL / 2πα |
|---|---|---|
| `mu` primary, constant doublet panel reaching x=0 | −0.95..−1.04 mid-chord, **+0.41 at the LE row** | 0.82 |
| same + sub-panel stack (SUBX=4) | **+22 mid-chord** (1/r² is not integrable: a midpoint stack does not converge to the kernel) | — |
| **horseshoe, loading primary (adopted)** | **−0.94..−1.03** (LE −1.13) | **0.852** |

(`2D CL / 2πα` for the adopted form is 0.852; against the AR40 finite-span
value `2πα*40/42` it is 0.895. The span-40 test still has tips.)

Matrix structure (AR1.6, 13x7): chord neighbours +1.74/+1.67, span
neighbours +0.81/+0.81, `cond(D incl. Kutta) = 83`. No tuning parameter is
in the physics: `kutta_w` scales the whole mu block (row scaling of a linear
block cannot move its root) and nothing else is fitted.

## 2. Test A — M0.3 symmetric, AR1.6 (dfvlr thickness, profile=None, ep4)

Reference 0.0545, `Gamma_LL` 0.0407, LL CL 0.0639.

| quantity | mode 6 |
|---|---|
| **CL** | **0.04601 (−15.6%)** |
| 2 sum(Gamma dy)/S (circulation CL) | 0.0489 |
| Gamma_mid | 0.0298 (73% of Gamma_LL) |
| tip/mid | **0.64** (LL 0.56; A11 0.58) |
| tangency residual | 6.6e-9 |
| Kutta row | 1.2e-8 |
| \|F\|max | 1.2e-4, nsup 0 |

Epoch stability: 0.0461 / 0.0460 / 0.0460 / 0.0460. Ablation
`WING3D_KUTTA_TE=0` (tangency in the last cell, no Kutta row): CL 0.0473,
i.e. the discrete Kutta cell costs **2.8%** — it is *not* what sets the
level.

**Verdict on A: improves on A11 (−32% → −15.6%) but FAILS the ±3%
criterion.** Section 4 attributes the remainder with numbers.

## 3. Test C — chordwise loading at mid-span (symmetric, mode 6)

    x     .043   .113   .187   .262   .340   .420   .501   .583   .667   .750   .833   .917   .999
    gamma .080   .052   .038   .029   .023   .018   .015   .012   .009   .007   .004  0.000  0.000

LE peak, monotone fall, `gamma_TE/gamma_max = 0.000` — **Kutta satisfied**.
Cambered (same nx, profile on):

    x     .043   .113   .187   .262   .340   .420   .501   .583   .667   .750   .833   .917   .999
    gamma .123   .061   .028   .008   .001   .013   .022   .043   .098   .167   .181  0.000  0.000

Aft-loaded: the peak is at x=0.83, not the LE. That is the reference
section's mean line, which rises to y=0.0205 at x=0.833 and then **plunges
back to 0.0003 at the TE** (`section_profile_xy`), so `dyc/dx = -0.26` near
x=1 and the aft cells are asked for positive loading while the Kutta row
demands zero loading on the last one.

## 4. Where the −15.6% comes from (measured, no tuning)

1. **The kernel is at the classical finite-wing answer.** `CL_circ/alpha`
   = **2.24/rad**; Helmbold/DATCOM for a rectangular AR1.6 wing is
   `2*pi*AR/(2+sqrt(AR^2+4)) = 2.20/rad` (+2%). The harness's own lifting
   line (`AR/(AR+2)` with PG) gives **2.93/rad** — that elliptic formula is
   the outlier for a rectangular planform, and the imprint derived from it
   (0.0639) is 31% above the classical value.
2. Span refinement of the kernel alone (7 → 16 → 32 span cells, standalone)
   changes CL_circ by **−1.1% / −4.6%**, so the 7-cell span grid is not the
   limit.
3. In the 2D limit the same lattice gives 0.852 of `2*pi*alpha` (0.895 of
   the AR40 finite-span value); the discrete Kutta cell accounts for 0.028
   of that (measured above) and the LE cell — 4.3% of chord carrying a
   singular loading — for most of the rest.
4. The field adds nothing to the level: the field's pressure CL is **94%**
   of `2 sum(Gamma dy)/S` (0.0460 vs 0.0489), and the coupled circulation is
   within 1% of the standalone kernel's, so phi and the kernel are
   consistent — the kernel, not the coupling, sets the amplitude.
5. What the free-space kernel does **not** contain: the finite box's
   farfield, thickness-induced plane downwash, and the reference's real-wing
   effects. A10's +2.0% came from a closure that read those field effects
   through a chordwise-biased estimator (the same row gave 11% of the
   cambered lift, so +2% there was not evidence of correctness).

## 5. Test B — M0.3 cambered dfvlr

Cambered LL: CL 0.2698, `Gamma` 0.1722 (A0=0.0486, A1=0.0679, cl2 0.519).

| run | CL | vs cambered LL | Gamma_mid | vs LL Gamma | cambered/symmetric |
|---|---|---|---|---|---|
| **mode 6, Kutta row** | **0.1065** | −60.5% | 0.0618 | **36%** | **2.31x** |
| mode 6, `KUTTA_TE=0` (colloc) | **0.2394** | **−11.3%** | 0.1464 | **85%** | 5.20x |
| A11 mode 5, Kutta | 0.1245 | −53.9% | 0.0672 | 39% | 3.4x |

A **consistent finite-span estimate** for this section, from the same
classical slope the kernel reproduces: `cl2 = 0.519`, `CL3 =
0.519 * (2.20/2pi) / beta(0.3) = 0.191`. Mode 6 with the Kutta row lands at
0.1065 = 56% of that; with collocation it lands at 0.2394 (+25% over the
estimate) — and in both cases the *circulation* is what the TE treatment
decides (36% vs 85% of the cambered LL). This is A11's diagnosis, now with
the kernel bias removed: the residual freedom is entirely in the aft 8% of
the chord, because this section's mean line demands its largest loading
there and then closes to zero.

**Verdict on B: sign, ordering and shape PASS (cambered is 2.3x symmetric,
LE-and-aft loaded, Kutta satisfied); the level does not converge** — it is
a property of the reference mean line at the TE, not of the kernel.

## 6. Test D — transonic dfvlr, ENTROPY=1, ep2, anchors 0.27186/0.26253/0.22694

| M | CL | err | Gamma_mid | tip/mid | nsup | \|F\|max | converged? |
|---|---|---|---|---|---|---|---|
| 0.850 | 0.10497 | **−61.4%** | 0.0618 | 0.72 | 73 | 4.3e-3 | tangency yes, field no |
| 0.925 | 0.10037 | **−61.8%** | 0.0618 | 0.72 | 1212 | 4.8e-2 | no |
| 1.000 | 0.10628 | **−53.2%** | 0.0618 | 0.72 | 67717 | 4.5 | no (supersonic pocket on the chord) |

**CL does NOT fall with Mach: it is flat (0.105 / 0.100 / 0.106), and
`Gamma_mid` is bit-identical at all three Mach numbers** — the geometric
kernel does not know the Mach number, so the circulation it sets is
Mach-independent, and the only Mach response left is the field's pressure
integration (a few %).

Ablation `WING3D_KUTTA_KERNEL_PG=1` (x → x/beta, beta capped at 0.3), same
grid and epochs:

| M | CL | err |
|---|---|---|
| 0.850 | 0.14151 | −48.0% |
| 0.925 | 0.16278 | −38.0% |
| 1.000 | 0.18105 | −20.2% |

PG makes the *level* better at M1.0 by luck and the *trend* worse (CL rises
+28% across the range where the anchors fall 16.5%). It is reported and
**not used**: picking it would be tuning, and a linear similarity cannot
produce a transonic lift collapse. (The `a12_*.json` files on disk are the
**default geometric kernel** — reproduced twice, bit-identical; the PG
ablation ran over the same filenames and its numbers live in
`tmp_a12_Dpg.log` and the table above.) **The structural cause is section
1's block-triangular Jacobian**: because the mu rows do not depend on phi,
the field's nonlinearity can never feed back into the circulation. Any
mode-6 answer at M≥0.85 is the incompressible circulation riding on a
transonic pressure field.

## 7. Test E — free_k and BL

* **free_k (mode 0, untouched)**, M0.85 ep4 entropy: k_mid **0.020**, CL
  **0.00378** — reproduces A11 exactly (0.020 / 0.00378). The legacy
  closure is still collapsed.
* **BL** (M0.85, ep2, entropy, `WING3D_BL=1`, mode 6): CL 0.10499 vs
  0.10497 (**+0.02%**), Gamma_mid 0.0618 unchanged, delta*/c max 1.24e-3,
  no separation. Neutral, as in A11 (−0.2%); with a frozen, kernel-set
  circulation there is no mechanism left for the displacement thickness to
  act on.

## 8. Verdict

1. **The kernel is right.** It came from the constant-doublet = vortex-ring
   identity, built from the codebase's own Biot-Savart kernel
   (`ckernels.loop_vel_block`). Validation: exact flat-plate distribution →
   `w = -0.94..-1.03 alpha` (vs the field estimator's −0.45..−1.7); 2D
   limit 0.852 of `2 pi alpha` (0.895 of the AR40 finite-span value,
   including the discrete Kutta cell); AR1.6
   `CL/alpha = 2.24/rad` vs the classical 2.20/rad. The A11-measured 3.7x
   chordwise bias is gone: the tangency equation is now chordwise
   consistent.
2. **The coupled Kutta does produce physical lift for both sections**: LE
   peak, monotone fall, `gamma_TE = 0.000`, elliptic span (tip/mid 0.64 vs
   LL 0.56), cambered 2.31x symmetric with an aft-loaded shape. Sign,
   ordering, spanwise shape and Kutta are all right, and the whole mu block
   is a well-conditioned (cond 83) linear solve with residuals 1e-8..1e-12.
3. **But the level is the classical free-space thin-plate value, 15.6%
   below the OpenFOAM anchor** (2.24/rad vs the anchor's 2.50/rad), and the
   harness's lifting-line imprint (2.93/rad) is 31% above the classical
   slope — the gap is a reference/model difference, measured in section 4,
   not a defect I can remove without reading the field's box and thickness
   effects back in. I did not tune it.
4. **Mode 6's structural cost is test D**: a phi-independent mu block means
   no transonic circulation response. CL is Mach-flat and cannot fall with
   Mach. That is a design limit of the *analytic* kernel, not of the row
   layout. Fixing it requires the tangency equation to be evaluated *in the
   field* — which is mode 5's row — with the A11 bias removed by a
   consistent kernel, i.e. a hybrid `w_field - w_field_kernel + D.mu` row.
   That hybrid is the next step this work defines, and I have not spent it
   here.

Files: `wing3d/tsd.py` (`_dlm_panels`, `dlm_kernel`, `dlm_tangency_residual`,
`mu_circulation` mode-6 branch, mode-6 gate), `tests3d/tsd_m_sup.py`
(mode-6 rows/reporting), `tests3d/tsd_a12_kernel.py` (standalone kernel
validation), `tests3d/tsd_a12_ladder.py` (ladders; imports A11's harness so
grid/alpha/anchors are identical). Results: `tests3d/a12_*.json`,
`tests3d/a11_*_k6_*.json`; logs `tmp_a12_*.log`.
