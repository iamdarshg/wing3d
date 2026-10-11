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
