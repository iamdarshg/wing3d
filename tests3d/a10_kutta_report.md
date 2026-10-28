# A10 — tangency row + Morino wake condition: what it measures

Status: **the closure row is now LIVE in both phi and Gamma (unlike A9's
dCp(TE), it has an O(1) converged gain and a non-uniform spanwise gain),
and for a symmetric section it converges to an elliptic loading. But its
amplitude is set by a chordwise REDUCTION that the formulation does not
contain: one equation per span cell cannot represent the chordwise
tangency problem, and the two defensible reductions of it differ by a
factor of −2.3 in sign.** Numbers below; nothing tuned (alpha 1.25°,
Re, transition, grid, coefficients all fixed).

Gates: `WING3D_KUTTA=4` (this work), `WING3D_MORINO=1` (Morino wake
condition, default on for mode 4), `WING3D_TANG_MODE=pt|gal` (chordwise
reduction). Modes 0–3 unchanged, default 0 = baseline bit-identical.

## 1. Row formulation

Unknowns stay `[phi, gamma_b]` — one circulation per span cell, with the
prescribed chordwise shape `F(x)` (cumulative flat plate, `F(1)=1`) and
thickness-only transpiration. Replacing A9's `dCp(TE)` row:

    row_j = wav(x_cp, y_j) - (dyc/dx - alpha)             ('pt', x_cp=0.75c)
    row_j = int A(x) [wav - (dyc/dx - alpha)] dx / int A  ('gal', A=dF/dx)

`wav = 0.5*(wu+wl)`, `wu = (phi[k0+1]-phi[k0])/dz`, `wl = (phi[k0-1]-phi[k0-2])/dzm`
(`tsd.TSDSolver.w_plane`). The wing is a potential discontinuity on the
cut, so `phi` is only a *stored* single-valued array; every estimator of
the plane normal velocity is contaminated by the jump. Measured, M0.3,
g80x32x32, AR1.6, mid-span (k0=15, dz_core 0.0147):

| estimator | value at gamma=Gamma_LL | note |
|---|---|---|
| wide extrapolation (planes k0±1, k0±2) | **+0.982 = 24.4 gamma** | = jump/dz_core: plumbing, 100x the physical downwash scale |
| first-order, upper side only | +0.0755 | one-sided, still carries ±J/2 |
| **symmetric first-order pair (`wav`)** | **−0.0358 = −0.82 gamma** | ±J/2 cancels to first order: the physical downwash scale |
| mass through the cut, `(dphi−J)/dz` | **−0.0285** | ≈ the transpiration; nearly independent of gamma |

The symmetric pair is the only estimator whose gamma-response is O(1).
`dphi(TE) - gamma_b = -5e-4` and `d(wav)/dGamma = 0` at frozen phi (both
rows read phi), so — as in A9 — only the *converged* gain is meaningful.

## 2. Measured row weight (not 50/1500)

`A_scale = dR_phi/dphi = 7.24e3`, `B_scale = dR_phi/dgamma = 4.62e3`
(median/max frozen-phi Jacobian entries, g80x32x32, M0.3). Converged
gains from three phi solves at gamma = 0, 1.0, 1.5 Gamma_LL
(`tsd_a10_probe3.py`):

| row | converged dR/dGamma per span station | weight |
|---|---|---|
| pt (0.75c) | −1.155 −0.959 −0.856 −0.822 −0.834 −0.907 −1.052 | 7.24e3/0.82 = **8.0e3** |
| gal | −0.815 −0.699 −0.635 −0.614 −0.620 −0.664 −0.747 | 7.24e3/0.62 = **1.09e4** |

The gain is O(1) (A9's was 1.18e2 *only* because dCp(TE) was dominated by
the spurious jump gradient) and it is **non-uniform in span (tips 1.16 vs
mid 0.82)**, i.e. the row is not span-diagonal.

**Morino wake condition** (`kutta_morino`): the wake doublet is set to
`phi_up(TE)-phi_lo(TE)` measured from the field instead of a copy of
`gamma_b`. Measured residual `Morino d = -1.9e-4` at `gamma=0.0248` (0.8%),
and `d(te_jump)/dgamma = 0.985`: the cut-cell *defines* the jump as the
face flux, so the condition is very nearly an identity and carries almost
no amplitude information (it is kept because it is the correct Helmholtz
statement, and reported because it is not the thing that closes the
system).

## 3. Ladder A — low Mach M0.3, AR1.6, g80x32x32, ep4

`WING3D_KUTTA=4`, symmetric section (`a10_dfvlr_r4_M0.300_k4_ep4_noprof.json`),
'pt' row. Gamma trajectory (mid-span; LL 0.0407, LL CL 0.0545):

| ep | Gamma_mid | CL | CL_circ | tangency | |F|max |
|---|---|---|---|---|---|
| 0 | 0.0242 | 0.0390 | 0.0399 | 1.8e-8 | 5.1e-4 |
| 1 | 0.0248 | 0.0397 | 0.0408 | 2.4e-9 | 1.4e-4 |
| 2 | 0.0248 | 0.0397 | 0.0408 | 5.0e-11 | 1.0e-6 |
| 3 | 0.0248 | 0.0397 | 0.0408 | 5.9e-12 | 4.5e-7 |

A converged root (the row residual falls to 1e-11, not a convergence
failure), at **61% of Gamma_LL and 73% of the LL CL** — a real positive
lift, elliptic loading (test B), but 27% low. The row's control-point
balance is exact: `wav = -0.0028 - 0.82*0.0248 = -0.0231` against the
required `-alpha = -0.0218`.

Cambered dfvlr section, same run (`..._k4_ep4.json`): Gamma walks to
**-0.0075** and CL to **-0.0052** — collapse. Cause measured, not
guessed: `dyc/dx(0.75c) = +0.0251` for this section, so the single 0.75c
station's geometric slope `(0.0251 - 0.0218) = +0.0033` is ~zero
incidence, and the row correctly (for that station) returns ~zero
circulation.

## 4. Test B — spanwise loading

`Gamma/Gamma_mid` across the 7 span cells (y = 0.664, 0.455, 0.231, 0, ...),
symmetric case:

| y/semi-span | 0.83 | 0.57 | 0.29 | 0 | 0.29 | 0.57 | 0.83 |
|---|---|---|---|---|---|---|---|
| this row (k4) | 0.63 | 0.85 | 0.96 | 1.00 | 0.98 | 0.88 | 0.67 |
| A9 (dCp row) | 0.974 | 1.000 | 1.000 | 1.000 | 1.000 | 0.978 | 0.978 |
| lifting line | 0.558 | 0.823 | 0.957 | 1.00 | 0.957 | 0.823 | 0.558 |

The row is no longer span-flat (A9's ±3%) — the loading is now
**elliptic in sign and shape** (tip/mid 0.65 vs 0.56 LL), because the
row's gain varies across the span. It is still under-elliptic (0.65 vs
0.56) because the response is dominated by the *local* 2D chordwise
kernel (gain 0.82 at mid) and only ~15% of it is the trailed 3D
contribution. Note the middle of the span is reproduced to 3%
(0.98/0.96 vs 0.96/0.96); the error is concentrated in the last span
cell at each tip, where the whole 3D tip closure lives in one cell.

## 5. Where it breaks — the chordwise reduction

One equation per span cell cannot represent a chordwise tangency
condition. The two defensible reductions give, for the same dfvlr state
(probe3, M0.3, gamma-scan roots `G* = -R(gamma=0)/gain`):

| reduction | root Gamma_mid | vs cambered LL 0.1722 |
|---|---|---|
| pt at 0.75c | **-0.0069** | negative: 0.75c camber slope (+0.025) ≈ alpha |
| gal (A-weighted) | **+0.0175** | 10%: the flat-plate A(x) weight is LE-singular (A(0.043)=4.25) and the LE camber slope is +0.145, so the section's real camber lift (thin-airfoil A1 = 0.148) is not what this weight samples |

So the amplitude is a property of the *reduction*, not of the physics:
reduction choice moves the root from −0.007 to +0.018. That is the
structural limit of the formulation, and it is the honest answer to
"is circulation now free": **it is an unknown with a well-conditioned,
phi-dependent, span-coupled equation — but the equation is
under-determined in the chordwise direction, so the selected amplitude
is set by an arbitrary reduction.**

What would remove it: a free chordwise jump (`gamma(i,j)` per wing cell,
nc*ng unknowns with nc*ng tangency rows — a doublet-lattice/Morino
system). Then the chordwise loading shape is the solution, Kutta enters
through the wake jump, and no reduction is needed. That is a
redesign of the sheet representation, not another row weighting.

## 6. Ladder B/C/D — transonic, BL, free_k

PLACEHOLDER_B
PLACEHOLDER_C
PLACEHOLDER_D

## 7. Verdict

* The tangency row is genuine: it depends on phi, it has an O(1)
  converged gain in Gamma, it is span-coupled, and the coupled Newton
  system converges to 1e-11 residual instead of running away.
* For a **symmetric** section it produces a physical, elliptic, stable
  solution: Gamma 0.0248 (61% of LL), CL 0.0397 (73% of 0.0545),
  tip/mid 0.65 (LL 0.56).
* For a **cambered** section the answer is set by the chordwise
  reduction (pt: −0.007, gal: +0.018, truth 0.172 for dfvlr), because one
  row per span cell cannot carry the chordwise tangency problem.
* The **Morino wake condition** is correct Helmholtz physics and is now
  enforced from the field, but measured it is nearly an identity
  (d/dgamma = 0.985) and is *not* what closes the system.
* The wide-extrapolation normal velocity that a naive tangency row would
  use reads `24.4*gamma` (jump/dz_core) — 100x the physical downwash
  scale. Any pressure/gradient-sampled closure built on the stored
  single-valued phi plane is dominated by that contamination; only the
  symmetric first-order pair across the cut cancels it.

Files: `wing3d/tsd.py` (`w_plane`, `te_jump`, `tangency_residual`,
`kutta_tangency`, `kutta_morino`, `tang_mode`), `tests3d/tsd_m_sup.py`
(mode-4 rows/history), `tests3d/tsd_a10_probe{,2,3}.py`,
`tests3d/tsd_a10_ladder.py`, results `tests3d/a10_*.json`, logs
`tmp_a10_*.log`.
