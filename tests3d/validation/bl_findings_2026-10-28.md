# Viscous-inviscid coupling (WING3D_BL) — anchor measurement, 2026-10-28

Agent A8. Protocol: TSD full potential, span 1.6, g80x32x32, ep8, lifting-line
imprint (`WING3D_IMPRINT=ll`), entropy correction (`WING3D_ENTROPY=1`), identical
configuration for the BL0 and BL1 rows. Only `WING3D_BL` differs.

## What is implemented

* Integral BL (`wing3d/ibl.py`) marched per span strip on BOTH surfaces from the
  inviscid edge velocity: Thwaites laminar (exact `theta^2 = 0.45 nu/Ue^6 int Ue^5 ds`,
  White `l(lambda)`/`H(lambda)` fits), Michel/short-bubble transition, Head
  entrainment turbulent closure (RK45) with Ludwieg–Tillmann `Cf`.
* Edge velocity from the inviscid surface gradient at the wing plane, densified
  13 -> 161 chordwise points (the TSD wing plane has only 13 chordwise stations;
  differentiating `Ue` on those directly makes Thwaites' lambda oscillate).
* Transition **fixed at x/c = 0.05** (`WING3D_BL_XTR`), Re_c = **1.8e7**
  (`WING3D_RE`) — derived from the anchor case freestream (T=300 K, p=100 kPa,
  U=300 m/s, chord 1 m -> rho=1.123, mu_Sutherland=1.848e-5).
* Coupling: displacement surface `y_s = y_c ± y_t ± delta*(x)` -> added transpiration
  `u_e d(delta*)/dx` in `_apply_wing_bc` (added to the upper, subtracted from the lower),
  one BL update per nonlinear epoch (8 total), under-relaxed 0.5.
* Separation (`H > 2.4` or `Cf <= 0`) recorded per strip in `s.bl_sep`, never clamped.
* Gate: `WING3D_BL=1`. Default off; the BL0 column below reproduces the previously
  reported errors, so the gate is behaviour-preserving.

IBL sanity check (flat plate, Ue=1, Re=1.8e7): theta(1) = 1.235e-3 vs the
0.036 Re^-0.2 correlation 1.274e-3 (-3%); Cf = 2.115e-3 vs 0.0592 Re^-0.2 = 2.095e-3.
So the BL itself is not the weak link.

## Results (CL, error vs known_anchors)

| case | M | anchor | BL0 (final ep) | BL1 (final ep) | BL0 mean last 4 ep | BL1 mean last 4 ep | dCL | epoch sd (1 sigma) |
|---|---|---|---|---|---|---|---|---|
| dfvlr_r4 | 0.850 | 0.3210697 | 0.31430 (-2.11%) | 0.31480 (-1.95%) | 0.31502 | 0.31428 | **-0.24%** | 0.25% |
| dfvlr_r4 | 0.925 | 0.3629209 | 0.31850 (-12.23%) | 0.31880 (-12.15%) | 0.31853 | 0.31817 | **-0.11%** | 0.23% |
| dfvlr_r4 | 1.000 | 0.3037588 | 0.32870 (+8.22%) | 0.32830 (+8.09%) | 0.32900 | 0.32875 | **-0.08%** | 0.10% |
| rae2822 | 0.850 | 0.2509912 | 0.20770 (-17.24%) | 0.20980 (-16.42%) | 0.20727 | 0.20775 | **+0.23%** | 0.58% |

Mean |error|: BL0 **9.96%** -> BL1 **9.98%**. Two anchors improve (dfvlr M1.000,
rae2822), two worsen (dfvlr M0.850/0.925); every delta is inside one sigma of the
epoch-to-epoch scatter, so **the effect is not resolved above the solver's noise**.

The "final epoch" column is quoted because it is the campaign's convention, but it is
misleading here: the rae2822 BL1 run's last epoch jumps (CL 0.2070 -> 0.2098 while
delta*/c jumps 1.03e-3 -> 1.68e-3), i.e. that +0.99% "improvement" is an unconverged
coupling step, not a measurement. The last-4-epoch mean is the defensible statistic.

## Physical direction

delta*/c max is only **1.1–1.7e-3** at Re_c = 1.8e7, so the added transpiration slope
is ~0.002 against alpha = 0.0218 rad: the term perturbs the geometric normal velocity
by ~10%, and the measured CL response is a few tenths of a percent at most.

Expected direction was "viscosity reduces lift and weakens the peak". Measured: 3 of 4
cases lose a *little* lift (consistent signs, -0.08 to -0.24%), rae2822 gains
(+0.23%, and its final-epoch value +0.99% is unconverged); `Cp_min` (upper, centerline)
becomes slightly MORE negative in all four (+0.010, +0.115, +0.013, +0.360 in |Cp|).
So on this evidence the displacement acts mainly as added camber/aft thickening with a
frozen circulation, not as the classical lift-reducing displacement.

## Why the magnitude has to be small — structural, not a bug

The campaign runs with `free_k=False`, i.e. the bound jump `k` and the wake `Gamma`
are **prescribed** (lifting-line imprint, k=1). The BL can only redistribute pressure;
it cannot reduce the bound circulation, which is the dominant channel by which wall
vorticity lowers CL on a wing. A BL that does not see a Kutta/`Gamma` feedback cannot
deliver a percent-level CL correction.

Sensitivity check (`WING3D_BL_SCALE`, diagnostic only, never a default): scaling the
displacement slope by 20 (delta*/c -> 0.024, already ~15x the physical value) changes
CL by -0.54% and Cp_min by -0.09. Extrapolating that sensitivity, closing a 12% CL gap
would need the displacement slope scaled ~400x, i.e. delta*/c ~ 0.7 — half a chord of
displacement thickness. No physically admissible wall vorticity is available to this
model.

## Honest verdict

1. The term is implemented correctly (verified against flat-plate integral-BL
   correlations) and it is live (the diagnostic gain moves both CL and Cp_min).
2. At the anchor configuration it moves CL by **<= 0.24%**, i.e. ~40x smaller than the
   2–17% anchor errors, and inside the epoch noise. It does **not** close the gap.
3. Wall vorticity is therefore not the missing physics at these conditions — consistent
   with Re_c = 1.8e7 producing delta*/c ~ 1e-3.
4. **Caveat on the grade itself**: the four anchor cases in `cases/openfoam/transonic/
   campaign/` are `simulationType laminar` with `type slip` walls — they are inviscid
   Euler-equivalent solutions, graded as if viscous. A correct viscous term is being
   scored against an inviscid reference. The genuinely viscous anchor leg
   (`transonic/3drans`, M6 SA-neg) is still `of_anchor_pending`.
5. Nothing here was tuned to the anchors: Re, transition, relaxation and alpha were
   fixed before the runs and left at their physical values.

## Separation

No separation was flagged in any of the 8 runs (`s.bl_sep` empty at every epoch,
both surfaces, all spans), so no lift contribution was limited or zeroed anywhere.
At Re 1.8e7 with transition at 5% the computed shape factor stays below the H = 2.4
threshold even downstream of the shock on the coarse 13-station chordwise grid.
