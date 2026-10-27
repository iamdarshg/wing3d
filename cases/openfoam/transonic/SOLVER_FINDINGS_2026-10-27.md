# Solver findings — 2026-10-27 (root-cause hunt)

## 1. The panel imprint was broken for any span (ROOT CAUSE of the "+45%")

`build_case` solved a 3D Morino doublet panel system to get the circulation
imprint (`wake_gamma['gamma'] = res['muw']`, `wing_jump`). That solve is
ill-conditioned and returns **2D-like circulation at every span**:

| span | AR | panel CL (M0.85, a1.25) | muw range |
|---|---|---|---|
| 1.6 | 1.6 | +0.626 | +0.13 … +0.39 |
| 3.2 | 3.2 | **−7.499** | −5.02 … −1.14 |
| 6.0 | 6.0 | **−1.281** | −1.02 … −0.08 |

Truncated SVD at rcond 1e-4 removes the blow-up but leaves CL≈0, i.e. the
physical lift is not in the well-conditioned subspace → formulation defect,
not conditioning. Consequence: **every TSD CL in the campaign was built on a
garbage imprint**, and the TSD's CL flipped sign/exploded with span
(+0.52 / −3.75 / −1.21 at span 1.6/3.2/6.0, M0.3 where there is no shock).

Fix (`c46ce46`): `WING3D_IMPRINT=ll` — Prandtl lifting line
(thin-airfoil section cl, finite-AR factor AR/(AR+2), elliptic γ).
Result: CL +0.055 / +0.117 / +0.218 at AR1.6/3.2/6.0 — monotone, physical.

## 2. Branch-cut jump was missing the span-cell width

The wake/wing jump was applied as `Gamma/dz` without the span-cell width
`dy`. `muw` is circulation PER UNIT SPAN (see `tests3d/diag_camber.py`:
`CL = 2*int(muw dy)/span`). Fixed: jump = γ·dy (`c46ce46`).

## 3. Prandtl-Glauert applied uncapped

`cl2/beta` inflated lift as M→1 (CL 0.39/0.56/0.72 **rising** with Mach at
span 6). PG is only valid below critical Mach; capped at M_crit=0.72
(`9699435`), which is the standard bounded-lift-growth rule.

### Error after 1+2+3 (span 6.0 vs the 2D OF anchors, like-for-like)
| M | before | after |
|---|---|---|
| 0.850 | +45.0% | −10.6% |
| 0.925 | +35.5% | −18.7% |
| 1.000 | +137.8% | +46.1% |
| mean |err| | 72.8% | **25.1%** |

## 4. The lift is NOT coming from circulation — it comes from transpiration

Split test at M0.85, span 6.0, LL imprint:
- with transpiration: CL = 0.2953
- without transpiration (alpha=0, profile=None): CL = **0.0080**

So the imprint contributes ~0.008 and the **wing transpiration BC does all
the lift**. Consequences:
- `free_k` can never work as wired: it scales an imprint that carries no
  lift, and the wing tangency BC does not depend on k, so the Kutta row
  (`TE dCp = 0`) has its only root at k→0 (measured CL 0.039–0.052).
- The TSD CL is therefore Mach-flat (0.287/0.295/0.311) while the truth
  rises then collapses (0.321/0.363/0.213 at M1.0).
- Keeping the wake-consistent farfield (`WING3D_FARFIELD=wake`) changes
  little (−10.8/−17.8/+40.6%), confirming circulation is not flow-generated.

**Implication:** the transonic lift loss cannot be reproduced by any knob on
the current architecture. It needs the structural bound-vortex unknown the
audit described (γ_b solved WITH phi and closed by Kutta, with the tangency
BC and the bound vortex coupled — not two independent lift sources), or a
non-isentropic (entropy) density correction.

## 5. The M1.00 anchor is itself suspect

The dfvlr fine band was terminated at M0.975, so M1.00 is only the 86k
corpus value (0.2126). The 86k branch runs **~28% low** vs fine at every
measured point (fine/86k = 1.20–1.33). Corrected M1.00 fine estimate ≈ 0.272
→ our +46% becomes ≈ +14%. The real M1.00 error is probably ~14%, not 46%.
