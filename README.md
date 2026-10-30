# wing3d — transonic TSD / full-potential wing solver, validated against OpenFOAM

A 3D transonic **small-disturbance / full-potential** wing solver, benchmarked
against compressible OpenFOAM. The goal is a solver that is **3–5× faster than
OpenFOAM at matched accuracy**, with a CL uncertainty band derived from its own
hyperparameters that lies **within** OpenFOAM's band.

**Start here: [`HANDOFF.md`](HANDOFF.md)** — it carries the verified state, the
open problems in priority order, the run commands, and every trap already paid
for. If you are resuming cold, read that first; it is the source of truth.

## Layout

| path | what |
|---|---|
| `wing3d/tsd.py` | the solver: conservative residual, transpiration BC, entropy correction, coupled-Kutta modes, surface-Cp decontamination |
| `wing3d/ibl.py` | integral boundary layer (Thwaites laminar / Head turbulent) |
| `wing3d/geometry.py`, `wing3d/solver.py` | wing/section geometry, 3D Morino panel solver (used only as an imprint/initial guess) |
| `wing3d/shock.py`, `wing3d/transonic.py` | shock/estimator layer (diagnostic; part of it is anchor-calibrated — see the audit) |
| `tests3d/` | runners, sweeps, ablation ladders, and the `a*_report.md` development trail |
| `tests3d/validation/` | ctest-style harness (`ctest_validate.py`), `published_refs.json`, `known_anchors.json` |
| `cases/openfoam/` | OpenFOAM cases: historic campaign corpus, and the **free-tip AR1.6 anchors** |
| `tools/of_mesh/` | body-fitted mesher pathway (NASA TMR `hcf_wing`), VM runbook, scripts |
| `HANDOFF.md` | **master handoff** |
| `cases/openfoam/transonic/VISCOUS_HANDOFF.md` | what must change to make the anchors viscous + turbulent |
| `cases/openfoam/transonic/BAD_BUT_HELPS.md` | ledger of unsupportable-but-helpful fixes with their measured cost |
| `cases/openfoam/transonic/ERROR_REDUCTION_SCOPE.md` | the campaign scope + targets |
| `VALIDATION.md` | earlier external-validation summary |

## Install

```bash
python -m pip install -r requirements.txt
```

Python 3.12 is what we use. `numpy` and `scipy` are required; `matplotlib` is
only for diagnostics.

OpenFOAM (only needed for the reference side): **ESI OpenFOAM 2406**
(`/usr/lib/openfoam/openfoam2406`), installable from
<https://dl.openfoam.com>. See `tools/of_mesh/vm/README.md` for the exact
setup used on the GCP worker.

## Run the solver

```python
import sys; sys.path.insert(0, '.'); sys.path.insert(0, 'tests3d')
from tsd_campaign import build_profile, build_section
from tsd_m_sup import run_sup

prof, sec = build_profile('dfvlr_r4'), build_section('dfvlr_r4')
L, xc, ix, j0 = run_sup(0.85, span=1.6, epochs=8, free_k=False,
                        profile=prof, section=sec)
print(L['CL'])
```

The validated configuration is `WING3D_IMPRINT=ll` plus `WING3D_ENTROPY=1`;
`WING3D_BL=1` adds the boundary layer, `WING3D_CPDECONTAM=1` the decontaminated
surface-Cp estimator (see below).

## Run the validation harness

```bash
python -X utf8 tests3d/validation/ctest_validate.py --suite geometry
python -X utf8 tests3d/validation/ctest_validate.py --suite flow --targets dfvlr_r4_of3drans_probe
```

## Honest status (2026-10-30)

* Subcritical agreement is good: dfvlr M0.85 **−3.8%** versus the free-tip
  anchor (band ±5.7%).
* **The solver does not yet predict transonic lift loss**: its CL is flat
  across Mach while the truth rises then collapses. This is the #1 open
  physics problem — see `HANDOFF.md` §3 P2.
* **The anchors are inviscid** (compressible Euler: slip wall, laminar, no
  layers). Making them genuinely viscous + turbulent is priority 1 — see
  `cases/openfoam/transonic/VISCOUS_HANDOFF.md`.
* The 3–5× speed target has not been measured.

## Harness environment variables

| var | effect |
|---|---|
| `WING3D_IMPRINT=ll` | Prandtl lifting-line imprint (replaces the broken panel imprint) — **required** for current results |
| `WING3D_ENTROPY=1` | non-isentropic (shock entropy) density correction |
| `WING3D_BL=1` | viscous–inviscid coupling: integral BL displacement as transpiration |
| `WING3D_CPDECONTAM=1` | decontaminated surface-Cp estimator (removes a 10–27% dz sensitivity) |
| `WING3D_KUTTA=6/7` | coupled bound-vortex circulation (experimental; Mach response still weak) |
| `WING3D_NONPLANAR=1` | per-(i,j) surface height (experimental; dihedral still fails) |
| `WING3D_ALPHA`, `WING3D_TE_GAP`, `WING3D_BC_ROWS`, `WING3D_FREE_K` | ablation knobs |
