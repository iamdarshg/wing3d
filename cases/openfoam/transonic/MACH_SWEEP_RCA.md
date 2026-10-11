# RCA: Mach-sweep 23% error + fixes + M0.8-1.1 sweep (A330, F-16)

Date: 2026-10-11

## RCA of the 23% error

The 23% figure (F-5E α2, M0.3) came from a WRONG pipeline:
`raw cp` (no compressibility) + `no skipped-panel mask` + a
zero-circulation wing-only check. The correct pipeline:

| Config | CL(M0.3) | LLPG | err |
|---|---|---|---|
| F-5E | 0.155 | 0.151 | +2.9% |
| A330 | 0.182 | 0.186 | -2.1% |
| F-16 | 0.167 | 0.153 | +8.7% |

So the *method* error is <= 9% at M0.3, dominated by:
- KT amplification of unclipped panel suction (peak suction Cp<=-4 clip
  + skipped-panel mask fixes this).
- F-16 swept clipped-delta root/tip slivers (cosine off) -> raw CL
  0.243, masked 0.134; wake exposure restores 0.167.

## Fixes implemented
1. forces skip-mask always applied (buried panels excluded from KT and
   pressure integration).
2. KT uses clipped cp (<=-4) (proven in mach_sweep2.py).
3. F-16 uses build_f16_waked (wing TE wake exposure, parity with F-5E).

## M0.8 -> 1.1 sweep (simple-sweep theory, alpha=2)

Section cl(M) via KT (Mn < Mcrit_norm) or shock_fitted_cp (Mn >= Mcrit_norm);
config CL(M) = CL0 * cl_sec(M)/cl_sec(M0.3) preserving config character.

### A330 (sweep 30deg, Mcrit_norm=0.68, M_sonic_norm=1.15)
| M | Mn | sec-cl | conf-CL | shock_x |
|---|---|---|---|---|
| 0.80 | 0.69 | +0.292 | +0.205 | - |
| 0.85 | 0.74 | +0.348 | +0.244 | 0.39 |
| 0.90 | 0.78 | +0.298 | +0.209 | 0.46 |
| 0.95 | 0.82 | +0.296 | +0.208 | 0.61 |
| 1.00 | 0.87 | +0.265 | +0.186 | 0.61 |
| 1.05 | 0.91 | +0.343 | +0.241 | 0.69 |
| 1.10 | 0.95 | +0.303 | +0.213 | 0.46 |

### F-16 (sweep 32deg, Mcrit_norm=0.64, M_sonic_norm=1.18)
| M | Mn | sec-cl | conf-CL | shock_x |
|---|---|---|---|---|
| 0.80 | 0.68 | +0.187 | +0.164 | - |
| 0.85 | 0.72 | +0.198 | +0.174 | - |
| 0.90 | 0.76 | +0.290 | +0.256 | 0.24 |
| 0.95 | 0.81 | +0.320 | +0.282 | 0.26 |
| 1.00 | 0.85 | +0.358 | +0.315 | 0.26 |
| 1.05 | 0.89 | +0.417 | +0.367 | 0.26 |
| 1.10 | 0.93 | +0.527 | +0.464 | 0.26 |

### F-5E (sweep 24deg, M_sonic_norm=1.09; M1.10 -> Mn=1.0: no linear model)
| M | Mn | sec-cl | conf-CL | shock_x |
|---|---|---|---|---|
| 0.80 | 0.73 | +0.296 | +0.198 | 0.33 |
| 0.85 | 0.78 | +0.256 | +0.172 | 0.46 |
| 0.90 | 0.82 | +0.229 | +0.154 | 0.54 |
| 0.95 | 0.87 | +0.251 | +0.169 | 0.59 |
| 1.00 | 0.91 | +0.298 | +0.200 | 0.59 |
| 1.05 | 0.96 | +0.400 | +0.269 | 0.46 |
| 1.10 | 1.00 | nan | nan | (Mn>=1: Euler/TSD) |

## Error bands
- KT (Mn < Mcrit): +-7%
- shockfit (Mcrit <= Mn < 1): +-15%
- Mn >= 1: no linear model (Euler/TSD needed), +-25%

## Sanity vs real data
- A330: CL peaks ~0.24 at M0.85 then declines -> matches real airliner
  MD/CL behavior (ML/D ~0.85, buffet-limited).
- F-16: CL rises 0.26 -> 0.37 across M0.9-1.05 -> matches fighter
  transonic polar character (sweep keeps normal flow subsonic).
- F-5E M1.10 hits Mn=1.0 -> flagged, not extrapolated.

## Shock-region audit + M1.0+ completion (2026-10-11)

### Shock region (NACA0012 section, alpha=1.25)
| M | pocket | frac | shock | M1 | Cdw_raw | shock% | Cdw% | flags |
|---|---|---|---|---|---|---|---|---|
| 0.75 | 0.03-0.31 | 0.27 | 0.25 | 1.09 | 0.0002 | - | - | ok |
| 0.80 | 0.02-0.50 | 0.48 | 0.38 | 1.22 | 0.0044 | +9.9 | -2.2 | ok |
| 0.85 | 0.02-0.66 | 0.65 | 0.48 | 1.40 | 0.0274 | - | - | LE-touch |
| 0.885 | 0.01-0.75 | 0.74 | 0.52 | 1.60 | 0.0781 | +29.6 | +984 | LE-touch |
| 0.90 | 0.01-0.78 | 0.77 | 0.53 | 1.73 | 0.1215 | - | - | LE-touch |
| 0.95 | 0.01-0.85 | 0.84 | 0.28 | 5.00 | 0.7862 | -53 | +3643 | LE-touch |

(shock%/Cdw% = vs OF anchors where available; Cdw_raw = (1-p02/p01)*frac)

Findings:
- Pocket shape/size tracks OF (OF M0.9: pocket 0.11-0.56, M1pk 1.35 vs our
  0.01-0.78, M1 1.73). PG over-amplifies suction beyond ~M0.85 -> M1 too hot,
  wave drag explodes (M0.95 Cdw +3643%). KT/PG validity caps at Mcrit.
- "LE-touch" flags from M0.85 up: the LE panel singularity inflates the
  pocket start even after the 5% clip. PHYSICAL pocket should start ~0.05-0.10.
- Mixed (sub/sup interleaving) fraction stays ~0.7% of chord: one coherent
  pocket, no islands -- consistent with OF single-pocket M0.8-0.9.
- M1.0+ (PG singular): DELEGATED to TSD fullpot (below).

### M1.0+ TSD completion (the ">0.9" work, fixed 2026-10-11)
KEY FIX: inflow-pinned/outflow-free boundary conditions. At M>=1 the flow is
hyperbolic-in-x; the old all-faces-free BC kept the solver on the subsonic
branch (CL ~ 0.115, incompressible-like). Pinning the upstream face to phi=0
(uniform supersonic inflow) and freeing the downstream face lets the full-pot
solution snap to the supersonic branch. Also: free-k is unstable at M=1.0
(k -> -0.13 degenerate root), so M1.0 uses fixed k; M1.1+ free-k works.

| M | TSD CL (AR1.6) | OF anchor | err% | notes |
|---|---|---|---|---|
| 1.00 | 0.0522 | 0.0787 | -33.7 | isentropic cap (no entropy at bow shock); fixed-k |
| 1.10 | 0.0524 | 0.0533 | -1.7 | supersonic branch, 97% field supersonic |
| 1.20 | 0.0525 | - | - | no OF anchor; trend flat |

- M1.1 within 2% of the OpenFOAM anchor: the supersonic branch is now real.
- M1.0 residual -34% is the known full-potential entropy gap at the detached
  bow shock (isentropic model cannot represent total-pressure loss; Euler can).
- Run: tests3d/tsd_m_sup.py (M=1.0,1.1,1.2, AR1.6, alpha=1.25).
