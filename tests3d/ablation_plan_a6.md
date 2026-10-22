# A6 ablation ladder plan (2026-10-22) — TSD unsupportable/shaky items

Referee: tests3d/validation/ctest_validate.py (read + run, never edit).
Reference DB: published_refs.json (A5 owns; read-only here).

## Flow targets (from published_refs)
- onera_m6_case2308 — M6 (M=0.84, alpha=3.06). PUBLISHED 3D RANS SA-neg
  CL=0.268/CD=0.0083, but a SWEPT 3D WING. Our TSD is 2D-section
  (span-extruded). INVISCID. NEVER force-fit to 0.268; keep
  diagnostic/PENDING. A5 keeps status=of_case_pending.
- dfvlr_r4_of3drans_probe (M=0.925 a=1.25) — PENDING until A5 anchor.
- rae2822_of3drans_probe (M=0.85 a=1.25) — PENDING until A5 anchor.
- naca0012_of3drans_probe (M=0.85 a=1.25) — PENDING until A5 anchor.

Until anchors land: targets are MEASURED at (M, alpha), status MEASURED
(not PASS/FAIL). Ladder records CL per knob so the error is one
subtraction when anchors land.

## Ablation order (user protocol: unsupportable first, then shaky)

### U1 free-k closure (unsupportable) — knob WING3D_FREE_K=1
Fixed-k (free_k=False) freezes circulation at the incompressible panel
value; Kutta rows are computed then discarded (tsd_m_sup.py:63-67).
Audit verdict: UNSUPPORTABLE as a transonic CL model. Test free-k live.
History: campaign-era diagnostics claimed free-k non-viable (k collapse /
spurious low-circulation solution ~0.494, tangency residual grows) — the
measurement must REDO this claim on the current fast-epoch solver.

### U2 rank-deficient BC rows (unsupportable) — knob WING3D_BC_ROWS=1
~14% of DOFs have no equation (R[-1]=0, R[:,0,:]=0 etc.); values set by
the Krylov path. Replaced by real farfield rows (extrapolation outflow,
Dirichlet lateral/top/bottom).

### U3 savgol section smoothing (shaky/tuned) — knob WING3D_RAW_SECTION=1
tsdcli slices STLs then savgol-smooths (w=41) until anchor CLs matched
("rae +0.04%"). Test raw slopes (the audit measured raw slopes shifted
the branch CL 0.25 -> 0.16 at M0.85 on rae — repro on current code).

### U4 TE-gap forcing (shaky) — knob WING3D_TE_GAP=0
build_wing morphs sections to te_gap=0.002c; dfvlr (true TE 0.0051c)
gets CLOSED DOWN. Test no morph (raw TE).

### S5 mask-freeze per epoch (shaky implementation)
State selection grid/path-dependent (rae g120/g140/g160 three branches).
Not a knob per se: the state band {small-mask, large-mask} is REPORTED,
part of the uncertainty band, not tuned. The fixed-point/stability test
is the g120 ep24 metastability result (small-mask is a stable fixed
point, CL 0.25845). This is a band, not a point: reported.

## Grading decisions (user protocol)
- If error goes down after removing the item -> keep removal permanently.
- Else -> BAD BUT HELPS entry with measured cost (append
  cases/openfoam/transonic/BAD_BUT_HELPS.md), never silently kept.
- Baseline = fixed-k / zero-BC-rows / savgol-on / te_gap=0.002 (the
  as-committed solver). Each ablation changes exactly one knob.
- Since no 3D RANS anchors are gradeable yet, use the provisional dfvlr
  fine-branch (dfvlr_fine program, in flight) for DIRECTION only, and
  the campaign-Euler-slab history for context (NOT gradeable). Final
  verdicts are recorded pending A5's anchors, and the measured CL/knob
  table is what A5 grades against.

## RAM discipline
- Container: dfvlr_fine band sequential, gate 1.8 GB, ~5h/solve.
  Active: M0.975 (started 11:44 UTC); remaining: 1.000, 1.025, 1.050.
- My ladder: gate 2.5 GB, one solve at a time, waits up to 2h per solve.
  Yields to the container.
