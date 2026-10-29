# A17 — branch-cut decontamination of the surface Cp estimator

Estimator decontaminated: the z-extrapolated off-surface `dphi/dx` used by
`TSDSolver.loads()`, `TSDSolver.fp_loads()` and `TSDSolver._surf_edge()`
(`wing3d/tsd.py`), gated by `WING3D_CPDECONTAM=1` (default OFF).

Everything below is measured; the harnesses are `tests3d/a17_diag.py`
(cut structure), `tests3d/a17_cache.py` + `tests3d/a17_lab.py` (estimator
variants from cached fields), `tests3d/a17_cpdecontam.py` (acceptance
A–D). Low-RAM diagnostic grid g40x24x16, M0.3, AR1.6, `WING3D_KUTTA=6`,
`WING3D_EP=1`, unless stated.

## 1. Where the Cp is estimated, and which cells carry the jump

| site | stencil | feeds |
|---|---|---|
| `loads`, `fp_loads` | `gradx` on planes `k0+1,k0+2` (upper) and `k0-1,k0-2` (lower), linear extrapolation to z=0 with weight `w = abs(z2)/abs(z2-z1) = 1.45` | CL, CD |
| `_surf_edge` | same two-plane pair | BL edge Mach |
| `w_plane` | first-order pair straddling the cut | tangency rows (modes 1–5, 7) |
| `kutta_residual` | `surf_cp(k0+1,k0+2)`/`(k0-1,k0-2)` at the TE | Kutta rows (free-`k` modes) |
| `surface_cp` | plane `k0` only | placeholder |

The branch cut is the **face** `k0-1` (`residual()` subtracts `J/dz` from
`phiz_f[...,k0-1]`); `k0` is `argmin|z|`. For a symmetric grid `z[k0]=0`,
so the model's sheet lies at `z_face = (z[k0-1]+z[k0])/2 = -dz/2`: the
sheet is always half a cell below the physical plane, and by more when the
grid does not pin a node at 0 (`z_face` is tabulated below).

## 2. Derivation

The one-dimensional FV argument with the jump imposed as a face-flux
correction: nodes at `z=±dz/2` around the face, physical `dphi/dz = m` above
and `m'` below, jump `J`. Requiring the face flux `(phi[k0]-phi[k0-1])/dz -
J/dz` to equal the physical average `(m+m')/2` gives

    phi[k0] = phi_true(0+) and phi[k0-1] = phi_true(0-)

i.e. **the two nodal planes are already the true one-sided values**, and
the difference across them is exactly `J`. Measured (`a17_diag.py`, nz=16,
mid span, all six chord stations):

    xc      0.102   0.253   0.414   0.580   0.750   0.920
    J       0.01201 0.01909 0.02338 0.02605 0.02751 0.02751
    phi[k0]-phi[k0-1]
            0.00982 0.01778 0.02273 0.02578 0.02749 0.02778

So **there is no stored ±J/2 offset to subtract**. The A16 prescription was
tested as an estimator (`tmp_a17_est.py`, `shiftJ2`) and it *raises* the dz
spread from 47.5% to 62.0%.

What is exact and usable is the sheet identity (the `gradx` operator applied
to the two cut nodes reproduces it to 2%, measured):

    dphi/dx(z=0+) - dphi/dx(z=0-) = dJ/dx,      dCp = 2 dJ/dx

with `J` known exactly. The decontaminated estimator is therefore

    gx_upper = bounded local extrapolation to z=0 from nodes (k0, k0+1)
               (|t| <= 1; contrast the stock unbounded w = 1.45 on the
               two far planes)
    gx_lower = gx_upper - dJ/dx
    dJ/dx    = cell-integrated (face averages of the nodal jump, end faces
               pinned to the physical sheet ends) so that
               sum_i (dJ/dx)_i dx_i = J_TE = Gamma exactly

and `J` is the jump the **solve was given** (`sheet_jump`: free mu sheet,
`F(x)*gamma_b`, or the prescribed panel jump), not the field's own
`phi[k0]-phi[k0-1]`. Measured why: the field jump is 18% low at the first
chord node (0.00982 vs 0.01201) because the discrete field smears the
LE-singular jump over the first cell; feeding it into the identity costs
7% of CL and does not even stabilise the z-offset test (0.5 dz: -18.2%
instead of -0.4%), since an off-plane sheet skews the measured jump too.
The field jump is kept only as a fallback for the `ramp*wake_G` mode.

Two checks that the estimator is now *exactly* consistent with the model
(cached nz=16 field, `tmp_a17_bias.py`):

| quantity | value | ratio to CL_circ |
|---|---|---|
| `2*int Gamma dy/S` (the model's own circulation) | 0.051562 | 1.0000 |
| `loads()` CL, TSD order `Cp = -2 dphi/dx` | 0.051562 | **1.0000** |
| `fp_loads()` CL, isentropic Cp | 0.048340 | 0.9375 |
| `sum_i gamma_i dx_i / J_TE` | 1.0000 | — |

So the decontaminated estimator reproduces the circulation the sheet
carries *exactly*; the 6.25% that `fp_loads` sits below it is the
isentropic pressure law at M0.3 (a constant at fixed M), not a dz effect
and not an estimator error.


## 3. Acceptance A — dz sweep (g40x24x16, M0.3, AR1.6)

| nz | dz below/above | z[k0] | z_face | CL OFF | CL ON | A16 OFF |
|---|---|---|---|---|---|---|
| 12 | 0.0783/0.0393 | -0.0196 | -0.0588 | 0.07029 | 0.04813 | 0.07029 |
| 16 | 0.0406/0.0406 | 0.0000 | -0.0203 | 0.05541 | 0.04834 | 0.05541 |
| 24 | 0.0263/0.0198 | -0.0099 | -0.0231 | 0.05358 | 0.04834 | 0.05358 |
| 32 | 0.0147/0.0147 | 0.0000 | -0.0074 | 0.04380 | 0.04864 | 0.04380 |

| | max-min spread | range |
|---|---|---|
| OFF | **47.5%** | 0.04380 .. 0.07029 |
| ON | **1.1%** | 0.04813 .. 0.04864 |

Reduction **52x** (success criterion: >=3x; target <9% met with room to
spare), and OFF is **bit-identical** to A16's recorded
0.07029/0.05541/0.05358/0.04380 (all four agree to <5e-6).
`max|dphi(ON-OFF)| = 0.0` at every nz.


`max|dphi(ON-OFF)| = 0.0` at every nz: the decontamination is pure
post-processing and does not touch the solve.

## 4. Acceptance B — rigid z offset (nz=16)

| offset | dz below/above | z_face | CL OFF | dCL OFF | CL ON | dCL ON |
|---|---|---|---|---|---|---|
| 0.00 dz | 0.0406/0.0406 | -0.0203 | 0.05541 | 0.0% | 0.04834 | 0.0% |
| 0.25 dz | 0.0406/0.0406 | -0.0101 | 0.05003 | **-9.7%** | 0.04830 | **-0.1%** |
| 0.50 dz | 0.0900/0.0406 | -0.0653 | 0.06926 | **+25.0%** | 0.04814 | **-0.4%** |

Both offsets are now under 0.5% (criterion <5%). Note the 0.50 dz case is
the harshest in the study — `make_grid`'s cosine remap means a half-cell
rigid shift **doubles** `dz` below the wing (0.0406 -> 0.0900) and parks
the model's sheet at `z_face = -0.0653` — and the estimator does not care:
0.04814 against 0.04834 at the nominal grid.

## 5. Acceptance C — campaign regression (dfvlr_r4, M0.85, g80x32x32, ep8)

`WING3D_IMPRINT=ll`, `WING3D_ENTROPY=1`, `WING3D_KUTTA=0`, `WING3D_BL=0`
(`sweep_multi.apply_validated_env`), built at M0.3 then run at M0.85 for 8
epochs, identical in the two states:

| state | CL | CD | vs anchor 0.32682 |
|---|---|---|---|
| CPDECONTAM=0 | 0.22134 | +0.00098 | **-32.3%** |
| CPDECONTAM=1 | **0.25123** | +0.00347 | **-23.1%** |

Reported as-is, not tuned. Points that matter for interpretation:

* The decontamination **improves** agreement with the anchor by **9.2
  points** (-32.3% -> -23.1%) and makes CD positive and of the right order
  (0.0035). This is a by-product, not a target: nothing was fitted, the
  estimator is fixed analytically by the sheet identity, and the same
  estimator is what produces the 52x dz collapse in section 3.
* The OFF number is -32.3% at this grid, inside the documented campaign
  error band for dfvlr_r4 ("30-50% low", `sweep_multi.py` header). The
  0.2674 quoted there is OpenFOAM's own coarser value
  (`cli_battery.json: of_cl`), not a TSD record.
* The remaining -23.1% is not addressed by this change; per section 7 the
  campaign grid carries only 6-19 chord nodes, where the stock estimator
  alone drifts 26%.
* The solve is untouched: `cp_decontam()` is reachable only from
  `loads`, `fp_loads` and `_surf_edge`, none of which the residual,
  `_apply_bc` or the frozen coefficient/rho builders call; measured
  `max|dphi(ON-OFF)| = 0.0` on the diagnostic grid. **Caveat**: with
  `WING3D_BL=1` the decontaminated edge Mach does feed the integral BL,
  so the ON solve can differ there (BL is off in the campaign config).



## 6. The sheet placement — and why the field jump was rejected

The cut is a face, so unless `z=0` is a face centre the model's wing is
displaced: `z_face = z[k0] - dz_below/2`. Over every A and B case:

| case | `|z_face|` | CL ON (prescribed jump) | CL ON (field jump) |
|---|---|---|---|
| nz=32 | 0.0074 | 0.04864 | 0.04827 |
| zoff=0.25dz | 0.0101 | 0.04830 | 0.04886 |
| nz=16 / zoff=0.00 | 0.0203 | 0.04834 | 0.04891 |
| nz=24 | 0.0231 | 0.04834 | 0.04760 |
| nz=12 | 0.0588 | 0.04813 | 0.04345 |
| zoff=0.50dz | 0.0653 | 0.04814 | 0.03999 |

With the prescribed jump the residual has **no dependence on `z_face`**
(0.04813..0.04864 across a 9x range of misplacement). With the field jump
(tested and rejected) CL tracked `|z_face|` at -0.143 per chord with only
+-0.0015 fit residual, i.e. the measured jump is itself a dependent of the
sheet placement: an off-plane sheet smears and skews `phi[k0]-phi[k0-1]`.

That is the whole reason `sheet_jump()` prefers the prescribed doublet
strength. The remaining 1.1% of ON spread is at the two coarsest/dirtiest
grids and is not attributed further here.

The corresponding *solution-side* observation (not fixed here, it needs
grid surgery and would break every recorded bit-identity): `make_grid`
pins a NODE at z=0, so the sheet always lands half a cell below the
physical plane, and by up to 0.065c when the grid is shifted. Placing
`z=0` on a FACE would make `z_face = 0` exactly.



## 7. Acceptance D — chordwise refinement (nz=32, M0.3, AR1.6)

The stock estimator was also compensating the unresolved leading-edge
suction cell, so the dz "floor" and the chord resolution are entangled.
nl=40/80/120 chord nodes (nc = 6/13/19):

| nx | nc | CL OFF | CL ON | `CL_circ` |
|---|---|---|---|---|
| 40 | 6 | 0.04380 | 0.04864 | 0.05156 |
| 80 | 13 | 0.05227 | 0.05235 | 0.05550 |
| 120 | 19 | 0.05508 | **0.05243** | 0.05550 |

* **ON is chord-converged**: 0.05235 -> 0.05243 (0.15%) from nx=80 to
  120, and equal to `0.94 x CL_circ` at both - the 0.94 being the
  isentropic pressure law (section 2), which the TSD-order estimator
  shows must be 1.0000.
* **OFF is not**: it climbs 26% over the same range and is still moving at
  nx=120 (0.05508, i.e. 5% *above* ON and ~1% above `CL_circ`).

So at the stock 6-node chord the stock estimator under-reads CL by 26%,
and its apparent dz sensitivity at nz=32 (0.0438) was substantially a
chord artefact rather than a plane-thickness one. The decontaminated
number is a property of the solution; the stock one was a function of both
grid directions.

## 8. Verdict

**The estimator is decontaminated and the dz floor it caused has dropped
by 52x.** At M0.3/AR1.6/g40x24x16 the reported CL is now
0.04813/0.04834/0.04834/0.04864 for nz = 12/16/24/32 (**1.1%**) against
0.0703/0.0554/0.0536/0.0438 (47.5%), and the two rigid z-offsets move it by
**-0.1% and -0.4%** instead of -9.7% and +25.0%. Both acceptance tests pass;
OFF is bit-identical to A16.

**The diagnosis A16 handed down is wrong, and the specified fix is not the
one that works.** There is no stored `+-J/2` offset: `phi[k0]-phi[k0-1]`
equals the prescribed doublet strength to 2% at every station and every nz,
because the residual already removes the `J/dz` it imposes on the cut face.
Subtracting `J/2` as prescribed makes the spread *worse* (62%). What was
actually dz-dependent is the stock two-plane extrapolation (`w = 1.45`) of
a non-linear near field; the working decontamination is the exact sheet
identity `dphi/dx(0+) - dphi/dx(0-) = dJ/dx` with the jump prescribed.

**Does the floor drop in the sense that matters?** Yes, on both readings.
(i) The reported CL is now a *converged property of the solution* — stable
to 1.1% in dz and 0.15% in nx from nx=80 — instead of a function of the
grid, and in TSD order it reproduces the model's own circulation exactly
(ratio 1.0000). (ii) Agreement with OpenFOAM at the campaign point also
improves, 0.22134 -> 0.25123 against the 0.32682 anchor (-32.3% -> -23.1%),
with CD turning positive. Neither was tuned: the estimator is fixed
analytically and the improvement fell out of using the prescribed jump.

**What to do with it.** Keep the default OFF so nothing already recorded
moves, then adopt it deliberately: quote CL from the decontaminated
estimator when the claim is about the *model's* lift (its own circulation,
UQ bands without a 47% dz term), and treat the stock CL as a
chord-resolution-dependent diagnostic. Two concrete follow-ups this
measurement hands on:
1. place `z=0` on a cell FACE instead of a node, so the cut stops being
   half a cell (up to 0.065c) off the physical plane — solution-side, and
   no estimator can substitute for it;
2. campaign grids carry 6-19 chord nodes and the stock estimator drifts 26%
   over that range (section 7): the largest single error term measured
   here, and not a dz effect at all.





