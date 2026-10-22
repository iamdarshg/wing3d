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