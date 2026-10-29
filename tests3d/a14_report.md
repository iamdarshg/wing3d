# A14 — BL symmetry, `B' = B + w_ref`, and the dz/nx scaling of the "1.36x"

Box saturated (0.4-0.7 GB free), so every number here is a **low-resolution
diagnostic**: absolute CL is NOT comparable to the g80x32x32 anchors; only
deltas and trends are. Mode-6 and the B-column (ref off) reproduce A13
bit-for-bit (0.04602 at g80 `ep2`; 0.06378/0.06210/0.06318 on the ladder), so
the harness is intact.

Gates: `WING3D_KUTTA=7`, `WING3D_KUTTA_THICK` (new: B'), `WING3D_THICK_EP`,
`WING3D_BL`/`WING3D_BL_SCALE` (unchanged).

## 1. BL: the antisymmetric entry was right; the non-cancelling part is small

`_apply_wing_bc` solves `wu = s/2 + v`, `wl = -s/2 + v`. A displacement
surface `y_u = y_c+y_t+d*_u`, `y_l = y_c-y_t-d*_l` therefore enters as

    source/2 = (ddu+ddl)/2  -> +above, -below -> CANCELS in wav  (correct:
                               a symmetric thickening has no mean-plane w)
    camber   = (ddu-ddl)/2  -> +both         -> the ONLY BL term wav sees

A13 read the antisymmetric entry as a bug; it is the correct sign (the lower
surface moves -z), and the non-cancelling term already existed. The split is
now explicit (`bl_wav_terms()`); the rewritten BC is **bit-identical** to the
old (`max|dR| = 0.0`). Measured at M0.3 symmetric, 80x24x24, ep2:
`d*/c max 1.16e-3`, `wav/alpha = +0.0258` (2.6% of alpha), source/2 3.9e-3.

| `BL_SCALE` | CL | vs scale 0 |
|---|---|---|
| 0 (== BL off) | 0.05682 | — |
| 1 (physical) | 0.05709 | **+0.47%** |
| 50 | 0.06782 | +19.4% |

So `dCL/d(delta*) > 0`: the physical viscous term **raises** CL by 0.5%.
"Viscosity reduces lift" is not a mechanism this formulation contains — the
lift-losing viscous effects (wake displacement, Kutta-point shift, effective
incidence loss) have no representative in a mean-plane doublet row; the only
route is upper/lower asymmetry, which at positive alpha adds camber.

## 2. B' = B + w_ref: removes thickness, does not remove the offset

`thickness_reference()` (tsd_m_sup) solves the SAME operator/BCs/farfield with
the bound+wake sheets disabled, and `w_ref = L.wav_ref` shifts every hybrid
row except the Kutta row. In the continuum a source sheet has `wav = 0`, so
this is a pure discretisation reference, not a tuned correction.

M0.3 symmetric, 80x24x24 ep4:

| | CL | resid_phys/alpha (mid) | wf/we |
|---|---|---|---|
| B (A13) | 0.06300 | -0.425 (-1.26..-0.12) | 1.65 |
| B' | 0.06211 | -0.350 (-1.19..-0.20) | 1.77 |

Transonic (80x32x32, ep2, ENTROPY=1):

| M | CL B | CL B' | resid B | resid B' | |F|max B' | nsup |
|---|---|---|---|---|---|---|
| 0.850 | 0.06378 | 0.06757 | -2.81..+0.11 | -2.38..+0.01 | 1.4e-1 | 80 |
| 0.925 | 0.06210 | 0.06869 | -2.99..+0.84 | -3.90..+1.51 | 9.2e-1 | 292 |
| 1.000 | 0.06318 | 0.08769 | -2.19..-0.14 | -1.36..+0.30 | 4.9 | 67472 |

B' removes the thickness *shape* (`w_ref/alpha` = -0.88 at x=0.04c, ~-0.1
mid-chord) but leaves the flat band: mean resid moves only -0.425 -> -0.350,
the field/kernel ratio gets worse (1.65 -> 1.77), and CL moves -1.4%, inside
the epoch wiggle (0.0595 -> 0.0621 across epochs, |F|max ~1e-2). **The flat
residual is the SHEET response, not thickness content** — so A13's own next
step cannot recover the level, and it does not. At M1.0 the reference solve
itself does not converge (|F|max 4.4), so that entry is void.

## 3. The 1.36x is an x-direction discretisation that shrinks with nx

Direct probe (`tsd_a14_norm.py`): prescribed exact sheet `mu = pi a F(x)`,
phi-only solve (converged |F| <= 1e-4), no mu coupling, no Kutta row:

    ratio = median(w_field / (E.mu))  over |E.mu| > 0.25 alpha, mid span

| grid | nc | z_up/c | z_lo/c | dz_below | ratio |
|---|---|---|---|---|---|
| 40x16x16 | 6 | +0.0203 | -0.0856 | 0.0406 | 1.648 |
| 56x32x32 | 9 | +0.0074 | -0.0246 | 0.0147 | 1.574 |
| 80x32x24 | 13 | 0.0 | -0.0585 | 0.0263 | 1.568 |
| 80x32x32 | 13 | +0.0074 | -0.0246 | 0.0147 | 1.590 |
| 80x32x48 | 13 | 0.0 | -0.0209 | 0.0101 | 1.622 |
| 120x32x32 | 19 | +0.0074 | -0.0246 | 0.0147 | **1.480** |
| 160x32x32 | 26 | +0.0074 | -0.0246 | 0.0147 | **1.394** |

* **nx refinement shrinks it monotonically**: 1.590 -> 1.480 -> 1.394 for
  nx = 80/120/160 (-12.6% for a 2x refinement), i.e. it is discretisation,
  converging slowly (extrapolates to ~1.2, not to 1).
* **dz_core/nz does not drive it**: 1.568/1.590/1.622 at nz = 24/32/48.
* It is **not a scalar**: symmetric cases give 1.57-1.65, cambered 0.40-0.64.
  A factor that flips sign with loading shape is an operator/shape mismatch,
  which is exactly why a fixed matrix `B'` cannot remove it.
* Coupled check: CL(mode7)/CL(mode6) at ep2 is 0.709 (nx=56) and 0.675
  (nx=120) — it does **not** close as the ratio shrinks, and these ep2 states
  are not tightly converged. No 1.36 factor was hard-coded anywhere.

## 4. Mach ladder with B': the fall still does not appear

0.06757 / 0.06869 / 0.08769 at M0.85/0.925/1.00 — CL **rises** with Mach
(+1.7% then +27.7%); vs anchors 0.27186/0.26253/0.22694 that is -75.2/-73.8/
-61.4% (unchanged to within a percent when the reference is off). The
non-circulatory content (`resid_phys/alpha` spanning -2.4..0.0 and
-3.9..+1.5) is still 5-10x the physical ~5% Mach signal, and B' does not
shrink it at 0.925 (it widens it). M1.00 is not a converged field (67472/80896
supersonic, |F| 4.9).

**Verdict:** the BL split is correct and its CL effect is +0.47% (sign
proven by the 50x liveness test); B' removes thickness but not the flat
offset, because that offset is the discrete sheet response; that response is
an nx discretisation (1.59 -> 1.39 as nx 80 -> 160) with no dz_core
dependence and no single value; and the transonic lift-loss slope still does
not appear — the circulation is set by a discrete downwash field whose
non-circulatory content is an order of magnitude larger than the Mach effect.

Files: `wing3d/tsd.py` (`bl_wav_terms`, `set_thickness_reference`,
`dlm_thickness_reference`, hybrid-row subtraction), `tests3d/tsd_m_sup.py`
(`thickness_reference`, records), `tests3d/tsd_a14_ladder.py`,
`tsd_a14_run.py`, `tsd_a14_norm.py`, `tsd_a14_analyze.py`, `tsd_a14_bl.py`;
results `tests3d/a14_*.json`, `a14_norm_M0.300.json`.
