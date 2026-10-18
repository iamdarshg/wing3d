# Error Reduction Scope — TSD vs OpenFOAM

Internal analysis: where the next few percent of error can come from, ranked
by evidence and cost. Companion to `CAMPAIGN_REPORT.md`; table mechanics live
in `campaign/final_table.*` + `tests3d/{final_table,gen_report_section}.py`.

Status (2026-10-18 ~14:30 host): **52/52 pairs, 10/52 within 10% as-run** —
the fills, CLI battery (12/12), Probe A (free-tip OF), and the grid/epoch
study are all complete and folded in below; dfvlr ceiling probes (86k
remesh) and g120/ep16 reconfirmation runs are in flight.

Metric: `err% = (CL_TSD − CL_OF)/CL_OF` at each (geo, M), α=1.25°, OF =
rhoCentralFoam Euler slab, TSD = fixed-k imprint full-potential, span 1.6.

---

## 0. Ranked levers

| # | Lever | Affected cells | Measured evidence | Cost | Expected |
|---|-------|----------------|-------------------|------|----------|
| L1 | **Config pairing**: OF is wall-at-tip, TSD is free-tip | naca/biconvex −20…−40%; rae S-shape | Probe A **done**: OF free-tip (naca M0.85) = 0.0629 vs wall 0.0728 → wall inflates +15.8%; free-vs-free residual +24% at campaign grid, **+9.7% at g120** | done (5.3 h) | ≈16 pp of the +44% naca gap is config; free-vs-free under 10% once grid-matched |
| L2 | **Corpus mesh seams** (mixed resolutions per row) | rae (48k vs 54k); dfvlr (45.6k vs 86k) | fills landed: dfvlr CL deltas across the seams run +1.3/+2.1/+2.3% — smooth through both transitions (§3) | done | **≤ ~1%; closed as immaterial** |
| L3 | **TSD FD grid (+ box, + epochs)** | all (the "few %" for close cells) | study **done** at M0.85: naca +13.3% (g120≈g160 within 0.9%), dfvlr −7.3% (g120≈g160 within 1.5%), rae +4.2% with g160 −13.2% (non-monotone); box inert ≤1.4% | study done; g120ep16 reconfirmations running | **the largest single lever: ±4…13% CL, sign per geo** |
| L4 | **TSD estimator + imprint resolution** | all; dfvlr wobble | ep16 at g80 moves toward the fine-grid answer (naca +4.2%, dfvlr −3.5%); drift evidence §1.3 | folded into L3 study | ≤ ~4%; subsumed by L3 |
| L5 | **OF scheme/refinement at the ceiling** | dfvlr +62…+276%; rae M1.0+ | rae peak 0.2778 @ M0.925; dfvlr peak 0.2791 @ M0.900 (new 86k); ceiling probes **running** (dfvlr M0.925 + M1.000 re-run on the 86k mesh) | ~5 h, in flight | may explain part of the dfvlr "saturation" |
| L6 | Mach-alignment accounting | all ±1.5 pp; dfvlr ±4.5 pp | aligned-M table (§1.2) | doc-only | not a fix; band note only |

---

## 1. Closed items (measured; not levers)

### 1.1 OF force-history sampling noise — CLOSED

Full force histories re-audited on both corpus eras (new-era runs from
container `postProcessing/`; old era from `campaign/container_coeffs/*.coeff`,
438 KB dumps that survived the postProcessing purge):

- new era: oscillation (max−min)/mean = **0.01–0.11%** (dfvlr 0.75/0.775
  sd 2–3e-5 on CL≈0.25; rae 0.85 0.01%; naca 0.04–0.11%).
- old era: ≤ 0.71% (naca 0.875), drift |last − window-mean| ≤ 0.19%
  everywhere, except known-defunct dfvlr 0.900 (tmax 0.009) and
  1.075/1.1 (tmax 0.076) — already excluded from the table.

Conclusion: the corpus values are converged signals to ≪1 pp of CL; no
averaging/noise reduction is available.

### 1.2 Mach-label convention — quantified band, not fixable in-corpus

OF rows run at thermodynamic Mach = **1.01683 × label** (corpus a_ref =
353.06 m/s vs a_thermo(T=300 K) = 347.22 m/s). Re-evaluating TSD CL at the
aligned Mach (`tmp_machscope.py`) moves the error column by:

- all cells: mean |Δ| **1.51 pp**, max 10.3 pp;
- rae: ≤ 0.7 pp except M1.00 **−6.8 pp** (31.9% → 25.2%);
- dfvlr: mean |Δ| **4.47 pp**, wobbling both signs (−9.6 pp at 0.875,
  +10.3 pp at 1.05 — alignment makes the worst cell *worse*).

Corpus freestream convention is uniform across all OF rows; changing the
labeling would not reduce errors systematically. Keep the §9 note.

### 1.3 TSD epoch convergence — small cleanup only

From `tests3d/camp_*.txt` epoch logs, drift |last − mean(last 4 epochs)|:

| geo | mean | max | worst cell |
|-----|------|-----|-----------|
| rae2822 | 0.18% | 1.0% | M0.95 (last dipped) |
| naca0012 | 0.34% | 2.1% | M0.90 |
| biconvex | 0.60% | 2.3% | M0.875 |
| dfvlr_r4 | 1.3% | **6.2%** | M0.90 (last-jump 6.9%) |

Estimator switch (window mean) or +epochs moves symmetric cells ≤ ~1–2 pp
and dfvlr a few pp inside +62…+276% errors. Folded into L3's ep16 config —
which confirmed epoch effects are second-order vs grid (§4); the
converged-increment stop rule remains the cheap fix for any follow-on sweep.

**CLI-battery corroboration (12 re-runs of the shipped path):** epoch
trajectories identify the two battery outliers as termination position, not
extraction bias: dfvlr M0.85's last epoch drops 0.5020 → 0.4558 (−9.2%,
landing −9.5% below campaign) while its 4-epoch window mean is −2.5% from
campaign; naca M1.05 drifts +2.7% across 8 epochs (0.0517 → 0.0531), the
campaign value sitting at the CLI's epochs 1–2 (+1.8% final). rae / biconvex
/ naca M0.75 agree ≤1%. A converged-increment stop rule (or window-mean
reporting) is the cheap systematic fix for any follow-on sweep — worth more
than the ep16 point alone.

---

## 2. L1 — Configuration pairing (the big structural lever)

**Geometry facts (all four meshes, verified from `/meshCache/*/boundary` and
case `checkMesh`):** box y = ±0.8 exactly; wing span 1.6 with tips *on* the
`symm` (type symmetry, 600 faces) y-planes. The OF "AR1.6 slab" is therefore
a **wall-at-tip cascade**, while the TSD campaign runs **free-tip AR1.6**.

**What the pairing costs, measured within-family (M0.85):**

| geo | OF (wall) | TSD (free) | OF/TSD − 1 |
|-----|-----------|------------|------------|
| naca0012 | 0.0728 | 0.0506 | **+43.9%** |
| biconvex | 0.0729 | 0.0484 | **+50.6%** |
| rae2822 | 0.2510 | 0.2491 | +0.8% |
| dfvlr_r4 | 0.2674 | 0.5036 | −46.9% (saturation, see L2/L5) |

Symmetric sections carry nearly all of their tip sensitivity (free AR1.6 to
wall-at-tip), cambered sections with strong shocks much less — the symmetric
family's −20…−40% was the prime config suspect:

**Probe A (done): OF free-tip case.**
`cases/openfoam/transonic/campaign_freetip/naca0012_M0.850_a1.25` — identical
to the corpus case except y ±2.4 (ny=9, Δy = 0.5333 matched) and the two
y-faces moved into the `farfield` patch (freestream BCs, true free tips at
1.6c clearance). Mesh: 116,192 cells, wall patch 4,708 faces vs the corpus
4,080 (matched/slightly finer; checkMesh OK, min volume 3.3e-6). Result:
**CL = 0.06287** at t=0.25 (12,525 steps, deltaT 1e-5; last-40-step wobble
±0.14% — converged; force history committed as
`tests3d/probe_naca_freetip_M085_coefficient.dat`).

**Decomposition (naca M0.85):**
- wall-tip inflation: 0.0728/0.0629 − 1 = **+15.8%** — config is ≈45% of the
  +43.9% gap (consistent with the earlier AR2.2 interpolation estimate);
- free-vs-free residual: 0.0629/0.0506 − 1 = **+24.2%**, most of which is
  TSD grid under-resolution (§4): at g120 (0.0573) the gap is **+9.7%**.

Verdict: config confirmed as a real but partial driver; config + grid
together take the naca M0.85 pairing under 10%. Cheap path is the
grid-corrected TSD numbers plus this decomposition (now in
`CAMPAIGN_REPORT.md` §7) — not a physical re-pair.

**Probe B (optional dev; NOT required by the Probe A verdict): TSD wall-mode
via method of images.** The span loading Γ(y) and the imprint come from the panel solve
(`wing3d/solver.py: assemble/solve`), with the TSD inheriting it through
`wing_jump` + `wake_gamma` (`tests3d/tsd_coupled.py: build_case`). The wall
at y = ±0.8 = mirror planes: add reflected copies of wing + wake to the
influence assembly and load/velocity evaluations, solving with mirror BCs
(period-1.6 cascade). Panel system is small (n_chord=24, n_span=12 build
inputs; 1152 panels), images multiply assembly cost only. Validation: image-
count convergence; large-gap limit must reproduce the free-tip result; touch
limit vs the known +44% (naca) / +0.8% (rae) OF deltas. Probe A already
settled the L1 question without this build — keep on the shelf unless a
physically wall-paired TSD sweep is ever wanted (e.g. a uniform-pairing
corpus rebuild).

---

## 3. L2 — Corpus mesh heterogeneity

Per-case polyMesh `nCells` (authoritative — from each case's own mesh):

| geo | M 0.75–0.825 | 0.850 | 0.875 | 0.900 | 0.925–1.10 |
|-----|--------------|-------|-------|-------|------------|
| naca0012 | 49,734 | 49,734 | 49,734 | 49,734 | 49,734 |
| biconvex | 54,822 | 54,780 | 54,780 | 54,822 | 54,780 |
| rae2822 | 48,114 | 48,114 | **53,946** | **48,114** | **53,946** |
| dfvlr_r4 | **86,145** (fills) | **45,606** | **45,606** | **86,145** (new) | **45,606** |

- naca uniform; biconvex uniform to 42 cells (0.08% — snappy rerun jitter).
- rae: seams at M0.875 (jump to 53,946), M0.900 (back to 48,114), M0.925+
  (53,946 again). Three resolution regimes in one row.
- dfvlr: the fills (0.75, 0.775, +filling 0.80/0.825/0.90) run on the
  86,145-cell mesh; the old band (0.85, 0.875, 0.925–1.10) on 45,606 —
  a **1.9× resolution seam** mid-row.

**Verdict (fills landed 2026-10-18): CLOSED — immaterial.** The dfvlr CL
curve reads 0.2549 / 0.2600 / 0.2641 / 0.2674 / 0.2729 / 0.2791 at M =
0.775 / 0.80 / 0.825 / 0.85 / 0.875 / 0.90 — consecutive deltas
+2.0/+1.6/+1.25/+2.1/+2.3%, i.e. smooth through both mesh transitions (no
step ≳1% off the local trend; the new 86k M0.90 value even becomes the row
peak). rae's error row shows the same smoothness across its seams (−0.8 →
−4.7 → −8.6 across 0.85→0.90). No same-M pair exists (the old M0.90 was
defunct), but mesh heterogeneity contributes ≤ ~1% CL — well inside the
L1/L3 levers and NOT worth the 24 h re-run program. The dfvlr "ceiling"
(0.278–0.279) and rae peak stay live suspects only for the *high-M
collapse*, which the L5 probes (§5) address.

---

## 4. L3/L4 — TSD discreteness (the "next few percent" on close cells)

Campaign TSD solve: FD full-potential on **80×32×32** (build_case default),
box xlim −8…10, y/z ±7, lateral faces free (zero residual), inflow pinned
φ=0. Study **done** (`tests3d/tsd_convergence.py`, results committed as
`tests3d/tsd_convergence.json`) — g80/g120/g160 + enlarged box (−12…16, ±11)
+ ep16, at rae / naca / dfvlr M0.85, 8 epochs each (16 for ep16):

| geo | g80 (= campaign; harness check) | g120 | g160 | box | ep16 |
|-----|----------------------------------|------|------|-----|------|
| naca0012 | 0.050603 (bit-exact) | 0.057318 (+13.3%) | 0.056814 (+12.3%) | +0.02% | 0.052715 (+4.2%) |
| rae2822 | 0.247272 (−0.7%) | 0.259523 (+4.2%) | 0.216096 (−13.2%) | −0.6% | 0.247046 (−0.8%) |
| dfvlr_r4 | 0.493878 (−1.9%) | 0.466888 (−7.3%) | 0.460139 (−8.6%) | +1.4% | 0.486212 (−3.5%) |

(±% vs the campaign value in `tsd_*.json`; naca g80 is bit-exact — the
solver is deterministic given identical inputs, so rae/dfvlr's ±0.7–1.9%
scatter comes from their STL-derived section inputs.)

**Findings:** (i) grid sensitivity is REAL and large — up to +13% (naca),
−7…−9% (dfvlr), non-monotone +4/−13% (rae); (ii) naca and dfvlr agree to
≈1% between g120 and g160 → converged values ≈ 0.057 (naca) and ≈ 0.46
(dfvlr); (iii) the enlarged box is inert (≤1.4%) — farfield truncation is
NOT a lever; (iv) ep16 at g80 moves in the same direction as grid
refinement (naca +4.2%, dfvlr −3.5%) — the campaign's 8 epochs were also
under-iterated. Consequences: campaign TSD values are not grid-converged
(carried as an as-run caveat in `CAMPAIGN_REPORT.md` §9); at M0.85 the
corrected pairings are naca: OF-free 0.0629 vs TSD-g120 0.0573 = **+9.7%**
(down from +43.9% as-run); dfvlr: 0.2674 vs ≈0.46 → −42% (the saturation
deepens — that error is NOT grid). A g120+ep16 reconfirmation (rae/naca/
dfvlr M0.85, plus naca M0.95/1.05 spot-checks) is in flight; a full
table-wide g120 re-sweep (≈40+ h host) is justified only if converged
values are needed for every cell.

L4 status: ep16 confirms the epoch sensitivity is real but second-order vs
grid (≤ ~4%); the 48×24 panel-imprint variant (generalize the strip
indexing in `build_case`) remains the untested L4 add-on.

---

## 5. L5 — OF scheme/refinement probes (dfvlr + rae high-M)

Hypothesis: the dfvlr row's +62…+276% and the rae M1.0+ collapse are
partly OF-side numerical behavior — central-upwind dissipation smearing
the shock on coarse meshes (45.6k cells dfvlr), not pure Euler physics.
Two single-case probes, staged after the current fills:

1. **Mesh refinement at the ceiling — RUNNING**: dfvlr M0.925 (the 0.2780
   peak) AND M1.000 (the collapse, 0.2126) re-run on the 86k cache, two
   container lanes launched 2026-10-18 (~5 h) — direct Δs vs the 45.6k
   values. (The same-M seam re-run at M0.85 was dropped: L2 closed as
   immaterial from the fills, §3.)
2. **Scheme variant**: same case with `fluxScheme Kurganov` (or equivalent
   less-dissipative central scheme) — isolates the numerics at fixed mesh;
   ~4–6 h.

If (1) moves CL ≥ ~3%, the mesh seam (L2) and the ceiling story merge into
one "resolution" program; if both are flat, the OF saturation is genuine
Euler behavior under this box/BC, and the dfvlr row is a config/physics
documentation case rather than an error to reduce.

---

## 6. Sequence and housekeeping

1. dfvlr fills → 52/52 — **done 2026-10-18**; §6/§7/§9 refreshed and
   committed (`ebec6d7`).
2. Probe A → L1 — **done** (CL = 0.06287; config ≈ 45% of the naca gap, §2).
3. TSD convergence study → L3/L4 — **done** (numbers in §4); g120/ep16
   reconfirmations + naca M0.95/1.05 spot-checks running on host.
4. Seam verdict from fills — **done, closed** (§3): mesh heterogeneity
   ≤ ~1%.
5. L5 ceiling probes (dfvlr 0.925 + 1.000 on the 86k cache) — **running**
   in the container; the scheme variant (L5-2) and rae high-M refinement
   are decided from their result.
6. CLI battery addendum to §2 — **done** (`4cce6c6`; 10/12 pairs ≤1%, two
   flags closed as termination position, §1.3).

Any TSD-side solver change (L1 wall-mode, L3 grid) automatically propagates
to `wing3d/tsdcli.py` (shared `run_sup` path) and would require regenerating
the campaign TSD sweep + the §6 table.

Evidence artifacts: `tests3d/scope_mach_align.py` (alignment band),
`tests3d/tsd_epoch_drift.py` (epoch drift), `tests3d/mesh_map.sh` (per-case
mesh map; runs inside the container), `tests3d/tsd_convergence.py` +
`tsd_convergence.json` (L3 study), `tests3d/tsd_conv2.py` (g120ep16
reconfirmations), `tests3d/make_freetip_probe.py` (Probe A case generator),
`tests3d/probe_naca_freetip_M085_*` (probe result + full force history). The
force-history audit is reproducible from the container `postProcessing/`
(all fills) and `campaign/container_coeffs/*.coeff` (old era).
