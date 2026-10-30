# Free-tip AR1.6 anchors — **INVISCID** (compressible Euler quality)

**Status: these are NOT viscous RANS results.** Read the caveats before using
them for anything.

Collected 2026-10-30 on GCP `cfd-worker-1` (e2-standard-8), OpenFOAM 2406,
`solver = rhoSimpleFoam`, free-tip box (wing y=±0.8, domain y=±4, all outer
faces farfield, **no symmetry planes**).

## What is and is not modelled

| property | status |
|---|---|
| compressible | **YES** — `hePsiThermo` + `perfectGas` + `sensibleInternalEnergy`, N2 JANAF; p/rho/T coupled with energy |
| viscous at the wall | **NO** — wall BC is `slip`; `nSurfaceLayers 0` |
| turbulent | **NO** — `simulationType laminar` |
| rotational | bulk only — no wall vorticity, no boundary layer, no friction drag, no shock/BL interaction |

So: compressible **Euler**. The CD column below is (numerical) wave drag only;
there is no skin friction in it.

## Anchors (converged to 3000 steady iterations unless noted)

| geometry | M | cells | CL | CD (inviscid) | status |
|---|---|---|---|---|---|
| dfvlr_r4 | 0.850 | 136,158 | 0.271861 | 0.054288 | converged |
| dfvlr_r4 | 0.925 | 136,158 | 0.262529 | 0.073529 | converged |
| dfvlr_r4 | 1.000 | 136,158 | 0.226938 | 0.092855 | converged |
| rae2822 | 0.850 | 138,035 | 0.171431 | 0.035337 | converged |
| rae2822 | 0.925 | 138,035 | 0.171657 | 0.051245 | converged |
| rae2822 | 1.000 | 138,035 | 0.155069 | 0.068163 | converged |
| naca0012 | 0.800 | 276,876 (ny30) | 0.048087 | 0.028344 | converged |
| naca0012 | 0.800 | 517,186 (ny60) | 0.050222 | 0.031382 | **still running** (2265/3000) |

Convergence evidence: identical CL on the last two sampled iterations for every
converged row (see `log.simple` in each case).

## Band

`±5.7%` as measured by B1 on these meshes: span ±5.2% (ny30 -> ny60, here
+4.4% on naca: 0.0481 -> 0.0502), domain ±2.1%, time ±1.0%.
Total O(±6%). **This band does NOT include the missing viscous physics**, which
is believed to be larger than the band (see VISCOUS_HANDOFF.md).

## Do not use these to claim TSD accuracy

The TSD solver is inviscid, so comparing it to these is an
inviscid-vs-inviscid comparison: internally consistent, but it validates the
solver against a reference that is missing the dominant real-flow physics at
M0.85-1.00 (shock/boundary-layer interaction, displacement thickness, friction).
Any "TSD agrees to N%" claim based on this file must carry that caveat.
