# Error Reduction Scope — TSD vs OpenFOAM

Internal analysis: where the next few percent of error can come from, ranked
by evidence and cost. Companion to `CAMPAIGN_REPORT.md`; table mechanics live
in `campaign/final_table.*` + `tests3d/{final_table,gen_report_section}.py`.

Status (2026-10-20 ~22:15 host): **52/52 pairs, 10/52 within 10% as-run** —
fills, CLI battery (12/12), Probe A, the grid/epoch study, the g120ep16
confirmation set (rae +4.0%, naca +12.5%, biconvex +12.8%, naca M0.95
+6.7%, naca M1.05 +12.9%, dfvlr −18.2% drifting) and both g160ep16 runs
(dfvlr −9.95%, still −1.0% in its last epoch; rae −17.45%, settled at
0.2054–0.2056) are complete and folded in below, plus the rae g140ep16
discriminator: CL 0.19683 (−20.97%), *below* both g120 (0.25906) and
g160 (0.20561) — non-monotone → supersonic-mask state selection, not a
grid trend (§4). One run queued behind a RAM gate: g120 ep24
metastability test (host at 0.41 GB free 2026-10-21 15:25; auto-launches
when free > 2.5 GB). Host
lanes were terminated by a machine restart at 2026-10-19 21:19 (no
campaign data lost); restarted 2026-10-20 13:48 under a 2-lane cap. L5: the first probe pair reproduced the corpus ceiling values
bit-exactly on the same (refined) mesh — a determinism check, not a mesh
test — and a preserved 45.6k copy (`/dfvlr_coarse`, not a campaign case)
failed checkMesh with the points/topology mismatch now **root-caused**:
the `constant` points file carries 49 stray points (59,360 vs the
required 59,311); the good points were in that scratch case's
`0/polyMesh/points`. The corpus itself is uniform 86k (§3), so the coarse
pair is a **reconstructed** sensitivity probe. Topology + good points
assembled → lanes verify "Mesh OK" → the coarse re-probe pair
(M0.925/M1.000) **completed**
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
| L2 | **Corpus mesh seams** (mixed resolutions per row) | rae (48k vs 54k); dfvlr none (uniform 86k) | fills landed: **dfvlr is uniform 86k at every M** — the 45.6k case dirs are failed-retry debris, no delivered value came from them (§3); rae's real seams are smooth | done | **closed: no dfvlr seam; mesh sensitivity itself is real but is L5's lever, not a corpus seam** |
| L3 | **TSD FD grid (+ box, + epochs)** | all (the "few %" for close cells) | study + g120ep16/g160ep16 confirmations **done** (4 geos, 8 runs: 6 g120ep16 + 2 g160ep16): rae +4.0% (stable), naca +12.5%, biconvex +12.8%, naca M0.95 +6.7%, naca M1.05 +12.9%, dfvlr −18.2% g120ep16 (was −7.3% at ep8 — the ep8 g120/g160 "agreement" was under-iteration), dfvlr g160ep16 −9.95% (still −1.0% in the last epoch), rae g160ep16 −17.45% (settled, 20.6% below g120), rae g140ep16 −20.97% (below both → supersonic-mask state selection, §4); box inert ≤1.4% | g120 ep24 RAM-gated (host RAM-starved) | **the largest single lever: −21…+13% CL, sign per geo** |
| L4 | **TSD estimator + imprint resolution** | all; dfvlr wobble | ep16 at g80 moves toward the fine-grid answer (naca +4.2%, dfvlr −3.5%); drift evidence §1.3 | folded into L3 study | ≤ ~4%; subsumed by L3 |
| L5 | **OF scheme/refinement at the ceiling** | dfvlr +62…+276%; rae M1.0+ | refined-mesh dfvlr peak 0.2791 @ M0.900 (reconstructed 45.6k probe: 0.3627 @ M0.925); M0.925/M1.000 probe pair reproduced the corpus values **bit-exactly** → the corpus ceiling values are already refined-mesh (86k) runs (solver deterministic); a preserved 45.6k copy (`/dfvlr_coarse`, **not a campaign case**) failed checkMesh (**root-caused**: constant points = 59,360 with 49 strays; good 59,311 points in that scratch case's `0/polyMesh`) → lanes fixed + **reconstructed** coarse pair **complete to t=0.25**: CL 0.3626993 / 0.2529717 vs refined 0.2779525 / 0.2126116 = **+30.5% / +19.0%** | done | **resolution matters** — dfvlr anchors are mesh-sensitive, not mesh-converged; corpus never used the 45.6k mesh; scheme variant (L5-2) unstaged |
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

Delivered-run `nCells` (reconstructed from per-case provenance; the dfvlr
case directories are stale — see the bullet below):

| geo | M 0.75–0.825 | 0.850 | 0.875 | 0.900 | 0.925–1.05 |
|-----|--------------|-------|-------|-------|------------|
| naca0012 | 49,734 | 49,734 | 49,734 | 49,734 | 49,734 |
| biconvex | 54,822 | 54,780 | 54,780 | 54,822 | 54,780 |
| rae2822 | 48,114 | 48,114 | **53,946** | **48,114** | **53,946** |
| dfvlr_r4 | **86,145** | **86,145** | **86,145** | **86,145** | **86,145** |

- naca uniform; biconvex uniform to 42 cells (0.08% — snappy rerun jitter).
- rae: seams at M0.875 (jump to 53,946), M0.900 (back to 48,114), M0.925+
  (53,946 again). Three resolution regimes in one row.
- dfvlr: **uniform 86,145 cells at every M (0.750–1.100)**. The 45,606-cell
  meshes now sitting in the 0.850/0.875/0.925–1.100 case directories are
  debris from a failed 2026-10-15 retry (`solver_df2.sh`,
  `mesh=dfvlr_r4.old`), not the solve: their
  `constant/polyMesh/points` is the 59,360-point corrupt file (49 strays),
  every `log.run` starts at the points-file mtime and ends without forces
  (`results2/*.dat` zero-byte; `batchv2.log`: `NOFORCES`), and the corrupt
  mesh dies at t≈4e-5 — yet every delivered history reaches t ≥ 0.075
  (most 0.25). The old-band results came from `dfvlr_runner.sh` ("all 11
  dfvlr cases ... with the new 86k mesh"); M0.925/M1.000 reproduce the 86k
  probe **bit-exactly**. So there is **no dfvlr resolution seam**, and the
  case-directory mesh is not the run mesh. (`log.check` is stale everywhere
  — 21,000 cells from 2026-10-12 — and is unusable for attribution.) Full
  evidence: `tests3d/dfvlr_mesh_audit.md`.

**Verdict (fills landed 2026-10-18; provenance-corrected 2026-10-20):
CLOSED for the delivered data — there is no dfvlr seam.** The dfvlr CL
curve reads 0.2549 / 0.2600 / 0.2641 / 0.2674 / 0.2729 / 0.2791 at M =
0.775 / 0.80 / 0.825 / 0.85 / 0.875 / 0.90 — consecutive deltas
+2.0/+1.6/+1.25/+2.1/+2.3%, i.e. smooth (the 86k M0.90 value is the row
peak). rae's error row shows the same smoothness across its real seams
(−0.8 → −4.7 → −8.6 across 0.85→0.90). Because dfvlr is uniform 86k, the
dfvlr "seam check" is a smoothness check on one mesh, not a mesh test: it
does **not** bound mesh sensitivity. L5's reconstructed 45.6k mesh gives
+30.5%/+19.0% vs the 86k anchors at M0.925/1.000 (§5), so the OF anchors
are mesh-sensitive even though the corpus never used that mesh. The dfvlr
"ceiling" (0.278–0.279) and rae peak stay live suspects only for the
*high-M collapse*, which the L5 probes (§5) address.

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

**Findings:** (i) grid sensitivity is REAL and large — up to +13% (naca);
rae moves +4.2% at g120ep8 but its settled g160ep16 value is 0.20561
(−17.45% vs campaign, 20.6% below g120); dfvlr is negative and grows with
epochs; (ii) naca
agrees to ≈1% between g120 and g160 → converged ≈ 0.057 (naca); the
dfvlr g120↔g160 ep8 "agreement" (≈0.46) was **under-iteration** — at ep16
g120 dropped a further −11.8% to 0.4119 and kept drifting (nsup 714→783
through ep15); (iii) the enlarged box is inert (≤1.4%) — farfield
truncation is NOT a lever; (iv) ep16 at g80 moves in the same direction as
grid refinement (naca +4.2%, dfvlr −3.5%) — the campaign's 8 epochs were
also under-iterated.

**The rae g160 "−13%" is under-iterated at ep8 and grid-dependent at
ep16.** At ep8 the g160 trajectory never settles: CL falls overall
through ep7 after an ep1→ep2 rise (0.2375 → 0.2426 → 0.2388 → 0.2294 →
0.2224 → 0.2196 → 0.2182 → 0.2161) while the supersonic count keeps
climbing (nsup 593 → 829 → 1015). The g160 ep16 follow-up completed and
settled at CL 0.20561 (−17.45% vs campaign) — 20.6% below the g120ep16
value (0.25906, +4.01%) — so rae's g160 state is grid-dependent, not the
g120 family.

**g140 discriminator (complete 2026-10-21 00:29): the grid trend is NOT
monotone — it is a supersonic-mask state selection.** rae M0.85 at
140×56×56, ep16: **CL 0.19683 (−20.97% vs campaign)** — *below* both the
g120 anchor (0.25906) and the g160 anchor (0.20561). The three-grid
sequence 0.2591 → 0.1968 → 0.2056 kills the smooth-refinement reading.
The trajectories separate into two states: g120 locks immediately into a
small-mask state (nsup 212→353, CL 0.258–0.261 from ep2), while g140 and
g160 both drift into a large-mask state (nsup → ~1226/1235) with CL
falling for 10+ epochs to 0.19–0.21. Within the large-mask state the
g140→g160 refinement moves CL by +4.4% (0.1968→0.2056) — a normal
discretization effect; the g120→g140 step (−24%) is the state switch.
Mechanism: each epoch freezes the supersonic mask and face coefficients
from the current solution and re-solves (`tsd_m_sup.py` `run_sup` loop) —
a lagged-coefficient fixed-point iteration whose basin is set by the ep0
seed, and the ep0 seed is grid-dependent (nsup at ep0: 212 / 1357 / 685).
Neither branch is demonstrably converged (nsup still creeping at ep15).
Consequence: rae's +3.2% (g120) and −18.1% (g160) are two states of an
unconverged nonlinear iteration, not a grid-convergence sequence; the
under-10% mid-band claim rests on the small-mask state. Next test
(queued): g120 extended to ep24 — if the small-mask state is a true fixed
point it stays at ≈0.259; if it drifts, all grids converge large-mask and
the "true" TSD rae answer is ≈0.20.

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
| rae2822 M0.85 g140ep16 | 0.24907 | **0.19683** | **−20.97%** | — | done — large-mask state (nsup→1226); below both g120 and g160 → non-monotone (§4) |

Consequences: campaign TSD values are not grid-converged (as-run caveat,
`CAMPAIGN_REPORT.md` §9). Corrected M0.85 pairings:
- **naca**: OF-free 0.0629 vs TSD-converged 0.0569 → **+10.4%** (from
  −30.5% as-run wall-paired; from +43.9% free-vs-free) — just over the
  10% line after both levers.
- **rae** (wall-paired all along): OF 0.2510 vs TSD 0.2591 → **+3.2%**
  (as-run −0.8%; the g120 grid correction moves rae *up*). **But the g160
  grid settles at 0.2056 (−17.5% vs campaign) — 20.6% below g120**, and
  the g140 discriminator lands even lower (0.19683, −21.0%): the
  three-grid sequence is non-monotone, so this is a supersonic-mask
  state selection (small-mask g120 vs large-mask g140/g160), not a grid
  trend (§4). The two state pairings: small-mask err% = **+3.2%**
  (OF/TSD − 1 = −3.1%), large-mask err% = **−18.1%** at g160 (−21.6% at
  g140). The rae under-10% mid-band claim rests on the small-mask state:
  rae is the case where solver state — not grid convergence — changes
  the verdict.
- **dfvlr**: 0.2674 vs 0.4535 (g160ep16, still −1.0% in its last epoch) →
  OF/TSD−1 = **−41.0%** (as-run −46.9%; the g120ep16 −35.1% was
  grid/iteration-under-converged, not convergence). The TSD value keeps
  moving with grid/epochs (0.5036 as-run → 0.4119 g120ep16 → 0.4535
  g160ep16, last still falling), so the dfvlr gap carries a ±10% band,
  not a settled number. The OF side is also mesh-sensitive (L5: a
  reconstructed coarse 45.6k probe gives +30.5%/+19.0% vs the refined 86k
  anchors at M0.925/1.000), so neither side is a settled reference.

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
test does not exist among the delivered values (the 45.6k meshes now in
the case dirs are later debris, §3). (Parsing note: the 12-column file is
time, Cd, Cd(f), Cd(r), Cl, Cl(f), Cl(r), Cm, ... — CL is column 5; the
compare script's CD column-7 read hits Cl(rear), ignored.)

**Attempt 2 (2026-10-19/20): coarse rebuild — ROOT CAUSE FOUND, pair COMPLETE.**
A preserved 45.6k copy (`/dfvlr_coarse`; **no campaign case has this
mesh** — the delivered dfvlr row is uniform 86k, §3) failed checkMesh
with a **points/topology mismatch, not geometry corruption**:
`constant/polyMesh/points` carries 59,360 points — the good 59,311 plus
**49 stray points** (md5 `b7adaff7e649fbd5951936cf0bd0dd27`) that shift
every index and produce the 11,606 "negative-volume" cells / 72,913
"wrongly oriented" faces in checkMesh. The good points (59,311) were in
that scratch case's `0/polyMesh/points`; `checkMesh` there had silently
combined the time-dir points over constant topology (hence the confusing
"Mesh OK" while every copy/consumer read the corrupt file and the solver
died at Time ≈ 4e-05). Fix (`tmp_lane_fix.sh`): assemble the good
topology (`constant/polyMesh`) + the good points (`0/polyMesh/points`) →
both lanes now verify **"Mesh OK" (59,311 pts / 45,606 cells)** and
`rhoCentralFoam` ran M0.925/M1.000 to t=0.25 (same protocol as the
corpus — past the old crash point, CL in `coeff_last.dat`-equivalent
`coefficient.dat`). Good-points md5 `ad7b845b4ae4da6e3cce44f936379a3a`;
faces/owner/neighbour/boundary were identical everywhere.

**Result (2026-10-20): resolution matters.** This is a **reconstructed**
sensitivity probe: the corpus never used the 45.6k mesh (all delivered
dfvlr values are 86k, §3). Both lanes reached `End` at
t=0.25 (the 2026-10-19 21:19 machine restart killed the first attempt at
t≈0.165; resumed from the 0.15 write-time dir after the reboot). Final
CL: **0.3626993 (M0.925)** and **0.2529717 (M1.000)** vs the refined-mesh
anchors 0.2779525 / 0.2126116 → **ΔCL = +30.5% / +19.0%**, far above the
≥3% "resolution matters" threshold. The reconstructed coarse mesh
*over-predicts* CL
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
the ceiling story merges with mesh sensitivity (L5), not with an L2 seam.
**Applied: +30.5% / +19.0% → resolution matters.** The L5 continuation is therefore an OF mesh-refinement program
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
   rae g140ep16 done (−20.97%, below both g120 and g160 → mask-state
   selection, §4). g120 ep24 metastability test queued behind a RAM gate
   (host at 0.41 GB free 2026-10-21 15:25).
4. Seam verdict from fills — **done, provenance-corrected**: dfvlr is
   uniform 86k at every M (no seam; the 45.6k dirs are failed-retry
   debris), rae's seams are smooth. Mesh sensitivity itself is L5:
   reconstructed coarse-vs-refined = +30.5%/+19.0%.
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
