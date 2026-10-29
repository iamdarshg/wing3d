# A6 final report (draft skeleton, fill as ladder completes)

Date: 2026-10-22. Agent A6 (TSD solver leg), wing3d TSD-vs-OpenFOAM.

## 0. Before-state (verified by ctest, this session)

### Geometry suite (instant)
- naca0012_geom  PASS  (t/c=0.1200 vs 0.120, camber 0.0000 vs 0.000, TE 0.0000)
- rae2822_geom   PASS  (t/c=0.1211 vs 0.121, camber 0.0126 vs 0.013, TE 0.0000)
- dfvlr_r4_geom  PASS  (t/c=0.1338 vs 0.135, camber 0.0210 vs 0.021, TE 0.0051)
- onera_m6_geom  PENDING (OF-side geometry; no TSD 2D proxy by design)
Result: 3/4 PASS, 1 PENDING -> results/a6_geom_baseline.json (exit 0).

### Flow suite
All 4 targets PENDING (not graded):
- onera_m6_case2308       -> PENDING (of_case_pending)  [diagnostic ceiling]
- dfvlr_r4_of3drans_probe -> PENDING (of_anchor_pending)
- rae2822_of3drans_probe  -> PENDING (of_anchor_pending)
- naca0012_of3drans_probe -> PENDING (of_anchor_pending)
Result: 0/4 PASS, 4 PENDING -> results/a6_flow_baseline.json (exit 0).
No flow anchor is gradeable until A5 publishes 3D RANS anchors.

### Solver state (HEAD, commit lineage e1309c4..778c60a)
- fixed-k imprint (free_k=False), rank-deficient BC rows (WING3D_BC_ROWS
  off), savgol STL smoothing in tsdcli, te_gap=0.002c morph.
- Container: dfvlr_fine band in flight (M0.975 solving; 1.000/1.025/1.050
  remain), RAM-gated 1.8 GB, ~5h/solve. My ladder gate 2.5 GB, yields.

## 1. Changes made (commits)
- e872b25  ablate prep: tsd_campaign geometry aliases (dfvlr_r4_true),
  onera_m6_sharp_te clean diagnostic raise, tsd_ablate.py ladder runner
  (ctest-protocol mirror, MEASURED status for PENDING), tsdcli
  WING3D_RAW_SECTION knob, tsd_conv2 WING3D_FREE_K env.
- 13925f8  BAD_BUT_HELPS.md ledger created (audit-derived items).
- 778c60a  ablation plan doc (tests3d/ablation_plan_a6.md).

## 2. Ablation ladder (fill from ablate_*.json)
- baseline: span 6.0 (ctest default), g80x32x32 ep8, probes MEASURED
- U1 free_k=1, U2 bc_rows=1, U4 te_gap=0, U3 raw section (tsdcli path)

## 3. BAD BUT HELPS (append as measured)
(fill)

## 4. Error distribution vs anchors
No 3D RANS anchor gradeable yet. Provisional direction only:
- dfvlr fine-branch (Euler 125k, in flight): M0.925 CL=0.3629209.
  TSD baseline at g80 (span1.6): ~0.4776 (campaign) -> +31.6% vs fine.
  (MUST re-measure at ctest span 6.0 and mark config caveat.)

## 5. Honest band
(fill after ladder)
