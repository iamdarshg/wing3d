# A10 — tangency row + Morino wake condition: what it measures

Status: **the closure row is now LIVE in both phi and Gamma (unlike A9's
dCp(TE), it has an O(1) converged gain and a non-uniform spanwise gain).
On an uncambered section it converges to a physical, stable, elliptic
loading: Gamma 0.0366 (90% of Gamma_LL), CL 0.0556 vs the 0.0545
reference. On a cambered section the amplitude is set by a chordwise
REDUCTION that the formulation does not contain: one equation per span
cell cannot represent the chordwise tangency problem, and the two
defensible reductions of it differ by a factor of −2.3 in sign.**
Numbers below; nothing tuned (alpha 1.25°, Re, transition, grid,
coefficients all fixed).

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

(The ladder runs below used w = 8.9e3, the pt-row value, which was the
number available when they started; the gal row's own measured weight is
1.09e4. A 22% row weight changes the Newton conditioning only, not the
root — the residual histories confirm it, tangency 4e-12 at w=8.9e3.)

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

`WING3D_KUTTA=4` (default 'gal' row), symmetric section
(`a10_dfvlr_r4_M0.300_k4_ep4_noprof.json`). Gamma trajectory (mid-span;
LL 0.0407, LL CL 0.0545):

| ep | Gamma_mid | CL | CL_circ | tangency | |F|max |
|---|---|---|---|---|---|
| 0 | 0.0354 | 0.0541 | 0.0585 | 8.8e-9 | 5.3e-4 |
| 1 | 0.0366 | 0.0555 | 0.0603 | 3.2e-9 | 5.3e-5 |
| 2 | 0.0366 | 0.0555 | 0.0603 | 5.1e-10 | 6.5e-6 |
| 3 | **0.0366** | **0.0556** | 0.0604 | 4.1e-12 | 1.4e-7 |

**Pass on the symmetric section**: `Gamma_mid` 0.0366 = **90% of
Gamma_LL**, CL = **0.0556 vs 0.0545 (+2.0%)**, from a converged root
(tangency residual 4e-12, |F| 1.4e-7). The pressure CL sits 8% below the
circulation CL (0.0604) — the sheet's lift is slightly over-predicted by
`2*int(gamma)dy/S`, the same 8-10% the A9 gamma-scan showed at
Gamma_LL (0.0587 pressure vs 0.0637 circulation).

The same run with the single-point ('pt') row instead gives
`Gamma_mid` 0.0248 (61%), CL 0.0397 (73%) — the chordwise reduction moves
the answer by 33%, see section 5.

Cambered dfvlr section, same settings (`..._k4_ep4.json`): Gamma walks
0.0178 -> 0.0188 and CL to **0.0296**, i.e. **11% of the cambered LL
0.1722** — and, wrongly, *less* than the symmetric case: the A-weighted
reduction samples the section's +0.145 leading-edge camber slope (table
in section 5), so its real camber lift (thin-airfoil A1 = 0.148) is not
what the row sees.

## 4. Test B — spanwise loading

`Gamma/Gamma_mid` across the 7 span cells (y = 0.664, 0.455, 0.231, 0, ...),
symmetric case:

| y/semi-span | 0.83 | 0.57 | 0.29 | 0 | 0.29 | 0.57 | 0.83 |
|---|---|---|---|---|---|---|---|
| this row (k4) | 0.638 | 0.856 | 0.962 | 1.000 | 0.978 | 0.885 | 0.673 |
| A9 (dCp row) | 0.974 | 1.000 | 1.000 | 1.000 | 1.000 | 0.978 | 0.978 |
| lifting line | 0.558 | 0.823 | 0.957 | 1.000 | 0.957 | 0.823 | 0.558 |

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

## 6. Ladder B — transonic dfvlr 0.85/0.925/1.00, span 1.6, g80, ep8, ENTROPY=1

New free-tip AR1.6 anchors (B1 rebuild): 0.27186 / 0.26253 / 0.22694.

| M | CL (this closure) | err vs anchor | Gamma_mid | tip/mid | converged? |
|---|---|---|---|---|---|
| 0.300 | 0.0296 | (no anchor) | 0.0188 | 0.65 | yes, \|F\| 9e-8 |
| 0.850 | **0.03865** | **−85.8%** | 0.0242 | 0.63 | yes, \|F\| 3.2e-5 |
| 0.925 | **0.03975** | **−84.9%** | 0.0234 | 0.64 | no, \|F\|max 3.6 |
| 1.000 | **0.02564** | **−88.7%** | 0.0132 | 0.65 | no, \|F\|max 6.8 (nsup 62k) |

CL sequence: 0.0296 (M0.3) → 0.0387 (M0.85) → 0.0398 (M0.925) → 0.0256
(M1.0). The required signature is a monotone fall 0.272 → 0.263 → 0.227
(−16.5%). What is measured is: **+31% compressibility lift growth up to
M0.85, +2.8% to M0.925, then a −35.5% collapse from M0.925 to M1.00** —
a transonic lift loss of twice the anchor's slope but starting one Mach
late, and with Gamma falling 0.0234 → 0.0132 (−44%) as the supersonic
pocket covers the chord (62k/82k nodes). A frozen circulation cannot do
that; this closure can, and does. But at **−85% level error** the
absolute agreement is not there, and M0.925/M1.0 are not converged
(\|F\|max 3.6/6.8, Kmax 5e-6/4e-5) — reported as such, not as physics.

Epoch stability: at M0.85 CL wanders 0.0307 → 0.0453 → 0.0387 over the 8
epochs while each epoch's inner Newton solve is converged (\|F\| 1e-4,
tangency 1e-9). The non-stationarity is the epoch re-freeze of the
supersonic mask/entropy operator (nsup 21 → 179), not the closure row.

## 7. Ladder C — dfvlr M0.85 + WING3D_BL=1

| BL | CL | Gamma_mid | dstar/c max | sep | Δ |
|---|---|---|---|---|---|
| 0 | 0.03865 | 0.0242 | — | — | — |
| 1 | 0.04034 | 0.0259 | 1.11e-3 | 0 | **+4.4% (+0.0017), Γ +7.0%** |

The BL now acts on a *physical* circulation (it moves Gamma itself by
7%, where A9's 5%-of-LL circulation barely responded), but the sign is
the wrong one for a viscous correction: displacement thickening
*increases* the circulation here by 7% instead of reducing it. With
dstar/c = 1.1e-3 (a Re 1.8e7 turbulent BL at 14% chord-worth of
thickening is plausible only for the short chordwise extent that is
resolved) and no separation, this is a live but small and
wrong-signed viscous effect — a finding, not a calibration.

## 8. Ladder D — legacy free_k (mode 0, the A9 row, unchanged)

| ep | k_mid | CL |
|---|---|---|
| 0 | 0.095 | 0.0381 |
| 1 | **0.128** | 0.0384 |
| 2 | 0.127 | 0.0386 |
| 3 | 0.128 | 0.0388 |

k = **0.128**, not ~1: the legacy `dCp(TE)` closure is untouched by this
work and still collapses (A9 measured the same 0.128). The healthy k≈1
is not reachable with that row; it is the same structural defect A9
diagnosed, and mode 4 does not repair it.

## 9. Verdict

* The tangency row is genuine: it depends on phi, it has an O(1)
  converged gain in Gamma (−0.61…−1.16 per span station), it is
  span-coupled, and the coupled Newton system converges to a 1e-11
  residual instead of running away.
* For a **symmetric** section it is a *pass*: Gamma 0.0366 (90% of
  Gamma_LL), **CL 0.0556 vs 0.0545 (+2.0%)**, loading elliptic
  (tip/mid 0.66 vs 0.56 LL) — test A and test B in one run.
* For a **cambered** section the selected amplitude is set by the
  chordwise reduction, not by physics (pt: −0.007, gal: +0.019, truth
  0.172), because one row per span cell cannot carry a chordwise
  tangency problem; the camber's real contribution (thin-airfoil
  A1 = 0.148) is not what an A(x)-weighted row samples (A-weighted
  camber slope +0.127, dominated by the +0.145 LE cell).
* **Morino wake condition**: correct Helmholtz physics, now enforced
  from the field rather than prescribed, but measured nearly an identity
  (d(te_jump)/dgamma = 0.985, residual 2e-4) — it is *not* what closes
  the system.
* Anything sampled from the stored single-valued phi plane is dominated
  by the jump: the wide extrapolation reads 24.4·gamma = jump/dz_core,
  100x the physical downwash scale. Only the symmetric first-order pair
  across the cut cancels the ±J/2 offset; that measurement alone
  rescues a tangency row from being meaningless.
* **Is circulation now free? Partly.** It is a real unknown with a
  well-conditioned, phi-dependent, span-coupled equation and a physical
  solution for an uncambered section — but the equation is
  under-determined in the chordwise direction, so any cambered section's
  amplitude is an artefact of the reduction. The fix is a free
  chordwise jump (`gamma(i,j)` per wing cell with nc·ng tangency rows —
  a doublet-lattice/Morino system), not another row weight.

Files: `wing3d/tsd.py` (`w_plane`, `te_jump`, `tangency_residual`,
`kutta_tangency`, `kutta_morino`, `tang_mode`), `tests3d/tsd_m_sup.py`
(mode-4 rows/history), `tests3d/tsd_a10_probe{,2,3}.py`,
`tests3d/tsd_a10_ladder.py`, results `tests3d/a10_*.json`, logs
`tmp_a10_*.log`.
