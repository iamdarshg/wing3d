# Error Reduction Scope — TSD vs OpenFOAM

Internal analysis: where the next few percent of error can come from, ranked
by evidence and cost. Companion to `CAMPAIGN_REPORT.md`; table mechanics live
in `campaign/final_table.*` + `tests3d/{final_table,gen_report_section}.py`.

Status at write time (2026-10-17 23:30 host): **49/52 pairs, 10/49 within
10%**; dfvlr M0.800/0.825/0.900 fills running; CLI battery finishing.

Metric: `err% = (CL_TSD − CL_OF)/CL_OF` at each (geo, M), α=1.25°, OF =
rhoCentralFoam Euler slab, TSD = fixed-k imprint full-potential, span 1.6.

---

## 0. Ranked levers

| # | Lever | Affected cells | Measured evidence | Cost | Expected |
|---|-------|----------------|-------------------|------|----------|
| L1 | **Config pairing**: OF is wall-at-tip, TSD is free-tip | naca/biconvex −20…−40%; rae S-shape | mesh boundaries + tips-on-wall geometry; OF/TSD ratio +44% (naca), +51% (biconvex) vs +0.8% (rae) | OF free-tip probe **running** (≈4–6 h); optional TSD wall-mode dev ≈0.5–1 d | reclassify symmetric family; potentially within ±10% under matched pairing |
| L2 | **Corpus mesh seams** (mixed resolutions per row) | rae (48k vs 54k); dfvlr (45.6k vs 86k) | per-case polyMesh notes (table below) | curve test from running fills; 1 same-M re-run ≈4 h | ±2–5% (TBD) |
| L3 | **TSD FD grid + box convergence** | all (the "few %" for close cells) | campaign grid 80×32×32; box xlim −8…10, y/z ±7; no refinement study yet | study queued (≈3 h host) | unknown; ±0.5–2% |
| L4 | **TSD estimator + imprint resolution** | all; dfvlr wobble | epoch drift ≤1.0% rae / ≤2.3% biconvex / 6.2% dfvlr | ep16 in L3 study; imprint variant small add-on | ±0.3–1 pp (dfvlr more) |
| L5 | **OF scheme/refinement at the ceiling** | dfvlr +70…+276%; rae M1.0+ | both rows peak at 0.2778–0.2780 at M0.925 on their coarse meshes | ≈4–6 h per variant (staged) | may explain part of the "saturation" |
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
and dfvlr a few pp inside +70…+276% errors. Worth folding into L3's ep16
config; not a standalone program.

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
wall-at-tip), cambered sections with strong shocks much less. The symmetric
family's −20…−40% may therefore be mostly **config, not physics** — but the
split cannot be proven from analysis alone; two probes decompose it:

**Probe A (running): OF free-tip case.**
`cases/openfoam/transonic/campaign_freetip/naca0012_M0.850_a1.25` — identical
to the corpus case except y ±2.4 (ny=9, Δy = 0.5333 matched), the two
y-faces moved into the `farfield` patch (freestream BCs, true free tips at
1.6c clearance), snappy caps 80k/95k, force constants copied verbatim
(magUInf 300.099, rhoInf 1.1231, Aref 1.6). Launched at nice 10 so the
dfvlr fills keep priority. Decision tree:
- OF_free(naca 0.85) ≈ 0.050–0.056 (within ±10% of TSD free) → **config
  confirmed** for the symmetric family; then either (a) build the TSD
  wall-mode and report wall-paired numbers, or (b) extend the corpus
  convention with free-tip OF rows for naca/biconvex (26 cases — only if
  the pairing is to be made uniform), else (c) document the mismatch.
- OF_free stays ≈ 0.07 → wall effect is NOT the driver; the −30% moves to
  OF-side investigation (L3/L5 analogues on the OF side).

**Probe B (optional dev, only if path (a)): TSD wall-mode via method of
images.** The span loading Γ(y) and the imprint come from the panel solve
(`wing3d/solver.py: assemble/solve`), with the TSD inheriting it through
`wing_jump` + `wake_gamma` (`tests3d/tsd_coupled.py: build_case`). The wall
at y = ±0.8 = mirror planes: add reflected copies of wing + wake to the
influence assembly and load/velocity evaluations, solving with mirror BCs
(period-1.6 cascade). Panel system is small (n_chord=24, n_span=12 build
inputs; 1152 panels), images multiply assembly cost only. Validation: image-
count convergence; large-gap limit must reproduce the free-tip result; touch
limit vs the known +44% (naca) / +0.8% (rae) OF deltas.

---

## 3. L2 — Corpus mesh heterogeneity

Per-case polyMesh `nCells` (authoritative — from each case's own mesh):

| geo | M 0.75–0.825 | 0.850 | 0.875 | 0.900 | 0.925–1.10 |
|-----|--------------|-------|-------|-------|------------|
| naca0012 | 49,734 | 49,734 | 49,734 | 49,734 | 49,734 |
| biconvex | 54,822 | 54,780 | 54,780 | 54,822 | 54,780 |
| rae2822 | 48,114 | 48,114 | **53,946** | **48,114** | **53,946** |
| dfvlr_r4 | **86,145** (fills) | **45,606** | **45,606** | **86,145** (filling) | **45,606** |

- naca uniform; biconvex uniform to 42 cells (0.08% — snappy rerun jitter).
- rae: seams at M0.875 (jump to 53,946), M0.900 (back to 48,114), M0.925+
  (53,946 again). Three resolution regimes in one row.
- dfvlr: the fills (0.75, 0.775, +filling 0.80/0.825/0.90) run on the
  86,145-cell mesh; the old band (0.85, 0.875, 0.925–1.10) on 45,606 —
  a **1.9× resolution seam** mid-row.

**Experiment already in flight:** when the 3 fills land, the dfvlr curve
0.75/0.775 (86k) → 0.80/0.825 (86k) → 0.85 (45.6k) → 0.90 (86k) shows
whether the old mesh under/over-resolves relative to the new one. If the
seam is material (≥2–3%), the clean fix is a stage-2 mini-campaign
re-running the 45.6k-era dfvlr cells on the 86k cache (6 cases ≈ 24 h
container) and/or the three 48.1k-era rae cells (0.875 already 54k; 0.90
and the 0.75–0.85 group) — decide from the seam delta. Note the dfvlr
"circulation ceiling" (0.278 @ M0.925) and the rae peak (0.2778 @ M0.925)
both sit on the coarse meshes — mesh dissipation is a live suspect (L5).

---

## 4. L3/L4 — TSD discreteness (the "next few percent" on close cells)

Campaign TSD solve: FD full-potential on **80×32×32** (build_case default),
box xlim −8…10, y/z ±7, lateral faces free (zero residual), inflow pinned
φ=0. Study script ready: `tests3d/tsd_convergence.py` — g80/g120/g160 +
enlarged box (−12…16, ±11) + ep16, at rae / naca / dfvlr M0.85, ~3 h host
(queued behind the CLI battery for host RAM). Reads the campaign value from
`tsd_*.json` and prints the delta, so the campaign replication (g80) doubles
as a harness check.

If ΔCL(g120−g80) or ΔCL(g160−g80) is material (≥1%), the campaign sweep can
be re-run on the converged grid (13 M × 4 geo × ~10–15 min at g120 ≈ 8–13 h
host, staged), or the delta applied as a documented correction basis. If
≪0.5%, the residual mid-band rae ±9% wiggle is not FD discretization —
attention moves to the imprint shape (L4) and config (L1).

L4 add-ons: ep16 (folded into the study) tests whether the dfvlr wobble
persists with more epochs; a 48×24 panel-imprint variant (generalize the
strip indexing in `build_case`) quantifies how much the coarse 24-station
Γ(y) imprint moves CL.

---

## 5. L5 — OF scheme/refinement probes (dfvlr + rae high-M)

Hypothesis: the dfvlr row's +70…+276% and the rae M1.0+ collapse are
partly OF-side numerical behavior — central-upwind dissipation smearing
the shock on coarse meshes (45.6k cells dfvlr), not pure Euler physics.
Two single-case probes, staged after the current fills:

1. **Mesh refinement at the ceiling**: dfvlr M0.925 (the 0.2780 peak) on
   the 86k cache — direct Δ vs the 45.6k value; ~4 h. (Doubles as the L2
   same-M seam measurement if run at M0.85 instead; prefer 0.925 to test
   the ceiling, or 0.85 for the seam — one case, pick by what the fills
   show.)
2. **Scheme variant**: same case with `fluxScheme Kurganov` (or equivalent
   less-dissipative central scheme) — isolates the numerics at fixed mesh;
   ~4–6 h.

If (1) moves CL ≥ ~3%, the mesh seam (L2) and the ceiling story merge into
one "resolution" program; if both are flat, the OF saturation is genuine
Euler behavior under this box/BC, and the dfvlr row is a config/physics
documentation case rather than an error to reduce.

---

## 6. Sequence and housekeeping

1. dfvlr fills (running) → 52/52; refresh §6/§7 via `final_table.py` +
   `gen_report_section.py` as they land.
2. Probe A result → L1 decision (report the config split honestly either way).
3. TSD convergence study (queued) → L3/L4 numbers; adopt or document.
4. Seam verdict from fills → optional dfvlr same-M re-run (L2/L5-1).
5. Optional scheme probe (L5-2) and rae high-M refinement if (4) motivates.
6. CLI battery addendum to §2 when its shell lands (`tests3d/cli_battery.json`).

Any TSD-side solver change (L1 wall-mode, L3 grid) automatically propagates
to `wing3d/tsdcli.py` (shared `run_sup` path) and would require regenerating
the campaign TSD sweep + the §6 table.

Evidence artifacts: `tests3d/scope_mach_align.py` (alignment band),
`tests3d/tsd_epoch_drift.py` (epoch drift), `tests3d/mesh_map.sh` (per-case
mesh map; runs inside the container), `tests3d/tsd_convergence.py` (L3
study), `tests3d/make_freetip_probe.py` (Probe A case generator). The
force-history audit is reproducible from the container `postProcessing/`
(all fills) and `campaign/container_coeffs/*.coeff` (old era).
