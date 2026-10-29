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

and `J` is read from the field (`phi[k0]-phi[k0-1]`), not re-derived from
the closure mode, so the estimator is mode-independent.  A plain nodal
finite difference of `J` instead of the cell-integrated form loses the
leading-edge cell (44% of `Gamma` on the stock 6-node chord): CL 0.0333 vs
0.0483.

## 3. Acceptance A — dz sweep (g40x24x16, M0.3, AR1.6)

| nz | dz below/above | z[k0] | z_face | CL OFF | CL ON | A16 OFF |
|---|---|---|---|---|---|---|
| 12 | 0.0783/0.0393 | -0.0196 | -0.0588 | 0.07029 | 0.04345 | 0.07029 |
| 16 | 0.0406/0.0406 | 0.0000 | -0.0203 | 0.05541 | 0.04891 | 0.05541 |
| 24 | 0.0263/0.0198 | -0.0099 | -0.0231 | 0.05358 | 0.04760 | 0.05358 |
| 32 | 0.0147/0.0147 | 0.0000 | -0.0074 | 0.04380 | 0.04827 | 0.04380 |

| | max-min spread | range |
|---|---|---|
| OFF | **47.5%** | 0.04380 .. 0.07029 |
| ON | **11.6%** | 0.04345 .. 0.04891 |

Reduction **4.9x**, and OFF is **bit-identical** to A16's recorded
0.07029/0.05541/0.05358/0.04380 (all four agree to <5e-6).

The 11.6% is not one number: it is entirely the nz=12 point. Splitting by
the model's sheet placement `z_face` (section 6) gives 2.7% for the five
grids whose sheet is within half a cell of z=0 (nz=16/24/32, zoff 0/0.25dz)
against 22.5% for the same five OFF — an **8.3x** reduction there.


`max|dphi(ON-OFF)| = 0.0` at every nz: the decontamination is pure
post-processing and does not touch the solve.

## 4. Acceptance B — rigid z offset (nz=16)

| offset | dz below/above | z_face | CL OFF | dCL OFF | CL ON | dCL ON |
|---|---|---|---|---|---|---|
| 0.00 dz | 0.0406/0.0406 | -0.0203 | 0.05541 | 0.0% | 0.04891 | 0.0% |
| 0.25 dz | 0.0406/0.0406 | -0.0101 | 0.05003 | **-9.7%** | 0.04886 | **-0.1%** |
| 0.50 dz | 0.0900/0.0406 | -0.0653 | 0.06926 | **+25.0%** | 0.03999 | **-18.2%** |

0.25 dz passes (<5%). 0.50 dz does **not**: it is still -18.2%. The cause is
measured, not mysterious: `make_grid`'s cosine remap means a half-cell rigid
shift **doubles** `dz` below the wing (0.0406 -> 0.0900), which puts the
model's sheet at `z_face = -0.0653` — 0.065c off the physical plane, the
worst placement in the study. `CL_ON` follows `z_face` at -0.143 per chord
(section 6), so -18.2% is what that placement produces, not an estimator
artefact. The offset that leaves `dz` unchanged (0.25 dz) is -0.1%.

## 5. Acceptance C — campaign regression

<!--C-->

## 6. Acceptance D — chord refinement

<!--D-->

## 7. Verdict

<!--V-->
