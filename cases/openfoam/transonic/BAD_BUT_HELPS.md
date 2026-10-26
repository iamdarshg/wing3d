# BAD BUT HELPS — TSD ablation ledger (wing3d vs 3D RANS truth)

Program rule (north star, 2026-10-22): physically unsupportable fixes that
*improve* a score are logged here with their measured cost; they must never
be silently kept. Removals that *lower* error are kept permanently.
Companion: cases/openfoam/transonic/AUDIT_2026-10-22.md (forensic verdicts)
and ERROR_REDUCTION_SCOPE.md §7 (targets).

Grading context: the TSD leg is inviscid full-potential by construction,
graded against 3D compressible VISCOUS (RANS) truth. The inviscid-viscous
gap belongs in the reported band, never tuned away to fake a PASS.

- Ablation protocol: per change, run ctest flow suite (probe targets +
  M6 diagnostic) with the knock-out knob, compare vs the anchored
  baseline, keep or BAD-BUT-HELPS.
- Targets are gradeable only once A5 publishes cl_pub/status in
  published_refs.json (3D RANS anchors); until then entries below are
  PROVISIONAL (direction only), not passes.

---

## TSD unsupportable-item ledger (audit items, commit e1309c4 lineage)

| # | item | status | measured cost | kept? |
|---|---|---|---|---|
| 1 | fixed-k imprint (`free_k=False`, circulation frozen at incompressible panel value; Kutta rows computed then discarded, tsd_m_sup.py:63-67) | UNSUPPORTABLE as a transonic CL model | (pending ablation) | audit: to test free-k live |
| 2 | rank-deficient BC rows (`R[-1]=0`, `R[:,0,:]=0` etc., ~14% DOFs with no equation; values set by Krylov path) | UNSUPPORTABLE | (pending ablation; knob WING3D_BC_ROWS=1) | audit: to test real farfield rows |
| 3 | savgol section smoothing in tsdcli (tuned until anchor CLs matched) | SHAKY / tuned | (pending ablation; knob WING3D_RAW_SECTION=1) | audit: to test raw slopes |
| 4 | TE-gap forcing to 0.002c in geometry.build_wing | SHAKY | (pending ablation; knob WING3D_TE_GAP=0) | audit: to test no morph |
| 5 | mask-freeze per epoch without fixed-point/stability test (state selection grid/path dependent; rae g120/g140/g160 three branches) | SHAKY implementation | (pending: state band, not a point) | audit: keep as band |

## Measurement log

(appended as ablations complete — each run: date, tag, knob state,
measured CL per probe, delta vs baseline, keep/BBH verdict)
## Ladder round 1 (span 1.6, g80x32x32 ep8, dfvlr probes) — 2026-10-25/26

Baseline (fixed-k, physical branch): CL 0.4655493 / 0.4918337 / 0.5054811
  err vs known anchors (fine125k / fine125k / corpus86k): +45.00% / +35.52% / +137.75%

| ablation | verdict | M0.85 | M0.925 | M1.00 | action |
|---|---|---|---|---|---|
| U1 free_k=1 | CATASTROPHIC | CL 0.004037 (-98.74%, **-143.74pp**) | 0.004879 (-98.66%, -134.18pp) | 0.005466 (-97.43%, -235.18pp) | REVERTED. k collapses to 0.020 -> lift dies. |
| U2 bc_rows=1 (real farfield rows) | MIXED | 0.438152 (+36.47%, **-8.53pp**) | 0.497836 (+37.17%, +1.65pp) | 0.508256 (+139.05%, +1.31pp) | NOT adopted: net win only at M0.85, small losses elsewhere. Kept as an option, not default. |
| U4 te_gap=0 | CATASTROPHIC | -1.373063 (-527.65%, **-572.65pp**) | (running) | (running) | REVERTED. Removing the 0.002c TE morph flips the solver to the degenerate negative-CL branch. |

BAD BUT HELPS entries (physically unsupportable, but removing them destroys the score):
1. **fixed-k imprint (k=1)** — on crystal-clear physics grounds (circulation never
   reselected; Kutta rows discarded) it is unsupportable. But U1 proves removing it
   (free_k live) collapses CL to ~0.005 at every Mach (-98.7%), because k -> 0.020.
   MEASURED COST OF REMOVAL: -143.7 / -134.2 / -235.2 pp of CL error.
2. **TE-gap 0.002c morph** — phys ically it edits the true (open-TE) R4 section.
   But U4 (te_gap=0) drives CL to -1.373 at M0.85 (sign flip to the degenerate
   branch). MEASURED COST OF REMOVAL: -572.7 pp. It is load-bearing for branch
   selection, not for the physics it claims.
3. **span 6.0** (from the earlier run) — the solver sits on the negative-CL branch
   (CL ~ -0.85). Not a knob we use; recorded so nobody re-enables it.

NOTE: the two catastrophes (U1, U4) and the span-6.0 disease are all the SAME
underlying defect: state selection (mask / branch) is unstable, so physically
correct changes flip the solution to a degenerate branch. That defect, not the
individual knobs, is what the solver work must fix.
