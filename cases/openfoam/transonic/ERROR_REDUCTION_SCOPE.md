# Error Reduction Scope — TSD vs OpenFOAM

Internal analysis: where the next few percent of error can come from, ranked
by evidence and cost. Companion to `CAMPAIGN_REPORT.md`; table mechanics live
in `campaign/final_table.*` + `tests3d/{final_table,gen_report_section}.py`.

Status (2026-10-20 ~22:15 host): **52/52 pairs, 10/52 within 10% as-run** —
fills, CLI battery (12/12), Probe A, the grid/epoch study, the g120ep16
confirmation set (rae +4.0%, naca +12.5%, biconvex +12.8%, naca M0.95
+6.7%, naca M1.05 +12.9%, dfvlr −18.2% drifting) and both g160ep16 runs
(dfvlr −9.95%, still −1.0% in its last epoch; rae −17.45%, settled at
0.2054–0.2056) are complete and folded in below. No lanes running. Host
lanes were terminated by a machine restart at 2026-10-19 21:19 (no
campaign data lost); restarted 2026-10-20 13:48 under a 2-lane cap. L5: the first probe pair reproduced the corpus ceiling values
bit-exactly on the same (refined) mesh — a determinism check, not a mesh
test — and the preserved coarse mesh's checkMesh failure is now
**root-caused**: the `constant` points file carries 49 stray points
(59,360 vs the required 59,311); the real snapped points were in
`0/polyMesh/points` all along. Topology + good points assembled → lanes
verify "Mesh OK" → the coarse re-probe pair (M0.925/M1.000) **completed**
(resumed from t=0.15 after the restart): CL = 0.3626993 / 0.2529717 vs
refined 0.2779525 / 0.2126116 → **+30.5% / +19.0% → resolution matters**
(§5). The dfvlr transonic rows are strongly mesh-sensitive; the corpus
anchors are refined-mesh values but are not demonstrably mesh-converged.

Metric: `err% = (CL_TSD − CL_OF)/CL_OF` at each (geo, M), α=1.25°, OF =
rhoCentralFoam Euler slab, TSD = fixed-k imprint full-potential, span 1.6.

---

## 0. Ranked levers

| # | Lever | Affected cells | Measured evidence | Cost | Expected |
|---|-------|----------------|-------------------|------|----------|
| L1 | **Config pairing**: OF is wall-at-tip, TSD is free-tip | naca/biconvex −20…−40%; rae S-shape | Probe A **done**: OF free-tip (naca M0.85) = 0.0629 vs wall 0.0728 → wall inflates +15.8%; free-vs-free residual +24.2% at campaign grid, **+10.4% at g120ep16** | done (5.3 h) | ≈16 pp of the +44% naca gap is config; free-vs-free ≈ +10.4% once grid-converged |
| L2 | **Corpus mesh seams** (mixed resolutions per row) | rae (48k vs 54k); dfvlr (45.6k vs 86k) | fills landed: dfvlr CL deltas across the seams run +1.3/+2.1/+2.3% — smooth through both transitions (§3) | done | **≤ ~1%; closed as immaterial** |
| L3 | **TSD FD grid (+ box, + epochs)** | all (the "few %" for close cells) | study + g120ep16 confirmations **done** (5 geos, 6 cells): rae +4.0% (stable), naca +12.5%, biconvex +12.8%, naca M0.95 +6.7%, naca M1.05 +12.9%, dfvlr −18.2% g120ep16 (was −7.3% at ep8 — the ep8 g120/g160 "agreement" was under-iteration), dfvlr g160ep16 −9.95% (still −1.0% in the last epoch), rae g160ep16 −17.45% (settled, 20.6% below g120 = shock branch); box inert ≤1.4% | no lanes running | **the largest single lever: −18…+13% CL, sign per geo** |
| L4 | **TSD estimator + imprint resolution** | all; dfvlr wobble | ep16 at g80 moves toward the fine-grid answer (naca +4.2%, dfvlr −3.5%); drift evidence §1.3 | folded into L3 study | ≤ ~4%; subsumed by L3 |
| L5 | **OF scheme/refinement at the ceiling** | dfvlr +62…+276%; rae M1.0+ | dfvlr peak 0.2791 @ M0.900; M0.925/M1.000 probe pair reproduced the corpus values **bit-exactly** → the corpus ceiling values are already refined-mesh (86k) runs (solver deterministic); the preserved 45.6k mesh's checkMesh failure is **root-caused** (constant points = 59,360 with 49 strays; real snapped 59,311 points in `0/polyMesh`) → lanes fixed + coarse pair **complete to t=0.25**: CL 0.3626993 / 0.2529717 vs refined 0.2779525 / 0.2126116 = **+30.5% / +19.0%** | done | **resolution matters** — dfvlr anchors are mesh-sensitive, not mesh-converged; scheme variant (L5-2) unstaged |
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
| dfvlr_r4 | 0.2674 | 0.5036 | −46.9% (as-run; L5: mesh-sensitive, §5) |

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
  TSD grid under-resolution (§4): at the converged g120ep16 value (0.0569)
  the gap is **+10.4%**.

Verdict: config confirmed as a real but partial driver; config + grid
take the naca M0.85 pairing from −30.5% as-run to ≈**+10%** free-vs-free
(**+10.4%** at g120ep16) — just over the line, the remainder being
free-vs-free numerics. Cheap path is the grid-corrected TSD numbers plus
this decomposition (now in `CAMPAIGN_REPORT.md` §7) — not a physical
re-pair.

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
  a **1.9× resolution seam** mid-row. *(2026-10-19 caveat: this per-value
  mesh attribution comes from the case mesh files at audit time; for the
  ceiling cells M0.925/M1.000 it is now known to be wrong — those values
  are 86k-mesh runs, §5 — and the dfvlr runner re-ran all 11 cases on the
  86k mesh; the seam-smoothness verdict below is unaffected.)*

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
≈+4% (rae, below), negative and growing with epochs for dfvlr; (ii) naca
agrees to ≈1% between g120 and g160 → converged ≈ 0.057 (naca); the
dfvlr g120↔g160 ep8 "agreement" (≈0.46) was **under-iteration** — at ep16
g120 dropped a further −11.8% to 0.4119 and kept drifting (nsup 714→783
through ep15); (iii) the enlarged box is inert (≤1.4%) — farfield
truncation is NOT a lever; (iv) ep16 at g80 moves in the same direction as
grid refinement (naca +4.2%, dfvlr −3.5%) — the campaign's 8 epochs were
also under-iterated.

**The rae g160 "−13%" is under-iteration, not grid behavior.** Its g160
epoch trajectory never settles: CL falls monotonically through ep7
(0.2375 → 0.2426 → 0.2388 → 0.2294 → 0.2224 → 0.2196 → 0.2182 → 0.2161)
while the supersonic count keeps climbing (nsup 593 → 829 → 1015). 8
epochs are not enough at 160×64×64; rae's fine-grid answer is its g120
family (≈ +4%, below) and a g160 ep16 run followed (launched 2026-10-20
15:25).

**Second-stage confirmations (g120 + ep16, `tests3d/tsd_conv2.py`):**

| run | campaign | g120ep16 | vs campaign | vs g120ep8 | status |
|-----|----------|----------|-------------|------------|--------|
| rae2822 M0.85 | 0.24907 | **0.25906** | **+4.01%** | −0.18% | done — stable (0.258–0.261 from ep2) |
| naca0012 M0.85 | 0.05060 | **0.05693** | **+12.49%** | −0.69% | done — 0.0568–0.0573 family (±1% wobble) |
| biconvex M0.85 | 0.04839 | **0.05460** | **+12.81%** | — | done — settled at 0.0543–0.0546 from ep9 |
| naca0012 M0.95 | 0.05157 | **0.05503** | **+6.71%** | — | done — window mean 0.0548 (M0.95 wobbles ±3%) |
| dfvlr_r4 M0.85 | 0.50364 | **0.41187** | **−18.22%** | **−11.78%** | done — NOT settled (0.5188→0.4119; −1.3/−1.6/+0.4% at ep12–15) |
| naca0012 M1.05 | 0.05220 | **0.05892** | **+12.89%** | — | done — flat 0.0589 from ep5 (nsup ~271k) |
| dfvlr_r4 M0.85 g160ep16 | 0.50364 | **0.45354** | **−9.95%** | — | done — 0.4875→0.4535; still −1.0% in the last epoch (NOT settled) |
| rae2822 M0.85 g160ep16 | 0.24907 | **0.20561** | **−17.45%** | — | done — settled 0.2054–0.2056 (ep14–15); 20.6% below g120 (grid-dependent state) |

Consequences: campaign TSD values are not grid-converged (as-run caveat,
`CAMPAIGN_REPORT.md` §9). Corrected M0.85 pairings:
- **naca**: OF-free 0.0629 vs TSD-converged 0.0569 → **+10.4%** (from
  −30.5% as-run wall-paired; from +43.9% free-vs-free) — just over the
  10% line after both levers.
- **rae** (wall-paired all along): OF 0.2510 vs TSD 0.2591 → **+3.2%**
  (as-run −0.8%; the g120 grid correction moves rae *up*). **But the g160
  grid settles at 0.2056 (−17.5% vs campaign) — 20.6% below g120**, a
  grid-dependent shock state (branch vs fine-grid artifact unresolved),
  so the finest-grid rae pairing is
  OF/TSD = **+22.1%**, not +3.2%. The rae +3.2% is a g120-grid result,
  and the rae under-10% mid-band claim is a coarse-grid result: rae is
  one of the cases where grid convergence changes the verdict.
- **dfvlr**: 0.2674 vs 0.4535 (g160ep16, still −1.0% in its last epoch) →
  OF/TSD−1 = **−41.0%** (as-run −46.9%; the g120ep16 −35.1% was
  grid/iteration-under-converged, not convergence). The TSD value keeps
  moving with grid/epochs (0.5036 as-run → 0.4119 g120ep16 → 0.4535
  g160ep16, last still falling), so the dfvlr gap carries a ±10% band,
  not a settled number. Not an OF-side grid artifact (L5, §5).

A full table-wide g120 re-sweep (≈40+ h host) is justified only if
converged values are needed for every cell; the M0.85 four-geometry block
is the representative set. The shift is not a low-M artifact: naca M1.05
(top of the band) moved **+12.89%** at g120ep16 (0.05892 vs 0.05220),
flat from ep5.

L4 status: ep16 confirms the epoch sensitivity is real but second-order vs
grid (≤ ~4%); the 48×24 panel-imprint variant (generalize the strip
indexing in `build_case`) remains the untested L4 add-on.

---

## 5. L5 — OF scheme/refinement probes (dfvlr + rae high-M)

Hypothesis: the dfvlr row's +62…+276% and the rae M1.0+ collapse are
partly OF-side numerical behavior — central-upwind dissipation smearing
the shock on coarse meshes, not pure Euler physics.

**Attempt 1 (2026-10-18) — invalid as a mesh test, kept as a determinism
check.** dfvlr M0.925 + M1.000 re-run in the container on the refined
cache (`/meshCache/dfvlr_r4`) returned CL = 0.2779525 / 0.2126116 —
**bit-identical to the corpus `container_coeffs/*.coeff` values** (all 12
columns). So: (i) rhoCentralFoam is run-to-run bit-reproducible here — a
useful determinism result; (ii) the corpus ceiling values were *already*
refined-mesh runs — the "45.6k vs 86k" pairing this probe was meant to
test does not exist in the corpus. (Parsing note: the 12-column file is
time, Cd, Cd(f), Cd(r), Cl, Cl(f), Cl(r), Cm, ... — CL is column 5; the
compare script's CD column-7 read hits Cl(rear), ignored.)

**Attempt 2 (2026-10-19/20): coarse rebuild — ROOT CAUSE FOUND, pair COMPLETE.**
The preserved coarse mesh's checkMesh failure was a **points/topology
mismatch, not geometry corruption**: `constant/polyMesh/points` carries
59,360 points — the good 59,311 plus **49 stray points** (md5
`b7adaff7e649fbd5951936cf0bd0dd27`) that shift every index and produce the
11,606 "negative-volume" cells / 72,913 "wrongly oriented" faces in
checkMesh. The correctly snapped points (59,311) were in
`0/polyMesh/points` all along; `checkMesh` at `/dfvlr_coarse` had silently
combined the time-dir points over constant topology (hence the confusing
"Mesh OK" there while every copy/consumer read the corrupt file and the
solver died at Time ≈ 4e-05). Fix (`tmp_lane_fix.sh`): assemble the good
topology (`constant/polyMesh`) + the good points (`0/polyMesh/points`) →
both lanes now verify **"Mesh OK" (59,311 pts / 45,606 cells)** and
`rhoCentralFoam` ran M0.925/M1.000 to t=0.25 (same protocol as the
corpus — past the old crash point, CL in `coeff_last.dat`-equivalent
`coefficient.dat`). Good-points md5 `ad7b845b4ae4da6e3cce44f936379a3a`;
faces/owner/neighbour/boundary were identical everywhere.

**Result (2026-10-20): resolution matters.** Both lanes reached `End` at
t=0.25 (the 2026-10-19 21:19 machine restart killed the first attempt at
t≈0.165; resumed from the 0.15 write-time dir after the reboot). Final
CL: **0.3626993 (M0.925)** and **0.2529717 (M1.000)** vs the refined-mesh
anchors 0.2779525 / 0.2126116 → **ΔCL = +30.5% / +19.0%**, far above the
≥3% "resolution matters" threshold. The coarse mesh *over-predicts* CL
(shock/LE-suction smearing), and the flow was settled long before t=0.25
(CL flat from t≈0.16; the pre-shutdown read at t≈0.1656 differs by
<0.05%). Consequence: the dfvlr OF anchors are refined-mesh values, but
the 30% coarse→refined swing shows they are **not demonstrably mesh-
converged** — the dfvlr high-M discrepancy cannot be attributed cleanly
without a refinement program, and the scheme variant (L5-2) is now
unstaged rather than merely pending. Evidence:
`tests3d/l5_coarse_dfvlr_M0.925_45k.dat`, `...M1.000_45k.dat`; driver
`tmp_l5_verdict.sh`.

Decision rule (unchanged): if the refined-vs-coarse ΔCL is < ~1%, the OF
saturation is genuine Euler behavior under this box/BC and the dfvlr row
is a config/physics documentation case; if ≥ ~3%, resolution matters and
the ceiling story merges with L2. **Applied: +30.5% / +19.0% → resolution
matters.** The L5 continuation is therefore an OF mesh-refinement program
(coarse → refined → finer on the dfvlr high-M points), not the scheme
swap; the scheme variant (same case, less-dissipative flux scheme, fixed
mesh; ~4–6 h) is unstaged.

---

## 6. Sequence and housekeeping

1. dfvlr fills → 52/52 — **done 2026-10-18**; §6/§7/§9 refreshed and
   committed (`ebec6d7`).
2. Probe A → L1 — **done** (CL = 0.06287; config ≈ 45% of the naca gap, §2).
3. TSD convergence study → L3/L4 — **done**; the full g120ep16 set landed
   (rae +4.01%, naca +12.49%, biconvex +12.81%, naca M0.95 +6.71%, naca
   M1.05 +12.89%, dfvlr −18.22% still drifting). Host lanes were
   terminated by the 2026-10-19 21:19 machine restart; restarted 2026-10-20
   13:48 under a 2-lane cap — dfvlr g160ep16 done (−9.95%, not settled);
   rae g160ep16 done (−17.45%, settled at 0.2054–0.2056, 20.6% below g120).
   No lanes running.
4. Seam verdict from fills — **done, closed** (§3): mesh heterogeneity
   ≤ ~1%.
5. L5 ceiling probes — **done** (§5). Attempt 1 (86k re-run) disproved the
   intended pairing by bit-exact reproduction; the coarse checkMesh failure
   root-caused as a points/topology mismatch (49 stray points in
   `constant/polyMesh/points`; good points in `0/polyMesh/points`); lanes
   fixed (`tmp_lane_fix.sh`), coarse pair completed to t=0.25 →
   **ΔCL = +30.5% / +19.0% → resolution matters**. Follow-on = OF mesh
   refinement program; the scheme variant (L5-2) is unstaged.
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
