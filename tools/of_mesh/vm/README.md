# GCP worker runbook (cfd-worker-1)

Set up 2026-10-29 because the local Docker/WSL container was saturating the
user's machine (WSL held 1.8 GB; the box had 0.7 GB free of 15.2 GB).

## Instance

| field | value |
|---|---|
| name | `cfd-worker-1` |
| zone | `asia-south1-c` |
| machine | `e2-standard-8` — 8 vCPU / 32 GB |
| disk | 80 GB pd-balanced |
| image | ubuntu-2204-lts |
| IP | 34.14.174.135 |
| cost | ~$0.27/hr, **on-demand** (preemptible CPU quota is 0 in this project) |

Created with a startup script installing: ESI **OpenFOAM 2406**
(`/usr/lib/openfoam/openfoam2406`, from dl.openfoam.com), `gfortran`,
`python3` + numpy/scipy, git.

**Stop it when idle** — `gcloud compute instances stop cfd-worker-1
--zone=asia-south1-c` — it bills while RUNNING only.

## What runs here

The heavy OpenFOAM anchor work (the part that was crushing the local box).
The Python TSD side is light (~300 MB) and can stay local.

## Reproduce from scratch

```bash
gcloud compute ssh cfd-worker-1 --zone=asia-south1-c
git clone https://github.com/iamdarshg/wing3d.git ~/cfd
bash ~/vm_build_all2.sh     # blockMesh + snappyHexMesh + checkMesh for the 8 cases
bash ~/vm_restart.sh        # patch fvSchemes/controlDict and launch all 8 solves
bash ~/vm_status.sh         # progress table
```

## Gotchas found the hard way (all fixed in vm_restart.sh)

1. **`~` is not expanded by `scp`** through gcloud/pscp — use
   `/home/iamdarshg/...` absolute paths.
2. **Meshes and STLs are not in the repo** (only dicts). Each case needs
   `constant/triSurface/<name>.stl` copied from
   `cases/openfoam/transonic/stl/` before `snappyHexMesh`. The STL name is
   read from the case's own `snappyHexMeshDict`.
3. **Background jobs die on SSH exit** even with `nohup` — use
   `setsid nohup ... < /dev/null &`.
4. **The committed `fvSchemes` is the rhoSimpleFoam variant**
   (`ddtSchemes { default steadyState; }`). rhoCentralFoam has
   `"rho.*" { solver diagonal; }` in fvSolution, and with steadyState the
   rho equation's diagonal is **zero** → `Floating point exception` in
   `diagonalSolver::solve` on step 1.
   Fix: `ddtSchemes { default localEuler; }` plus
   `fluxScheme Kurganov;` and `reconstruct(rho|U|T) vanLeer;`.
5. **`controlDict` had `deltaT 1`** → the time loop takes zero steps
   (`Starting time loop` then `End`). Set `deltaT 1e-07`,
   `adjustTimeStep yes`, `maxCo 0.1`, `maxDeltaT 5e-05`.
6. **rhoCentralFoam cannot start from a uniform freestream** on a sharp
   airfoil — it dies on the janaf thermo inversion
   (`Max iterations exceeded ... old T:10000 f:2e8`). Two known ways out:
   (a) gentle start (maxCo 0.1 / tiny dt, used here), or (b) what B1 did —
   converge rhoSimpleFoam first, then `mapFields` into the rhoCentralFoam
   case.
7. **`rm -rf [0-9]*/` deletes `0/`** (the initial fields) — it matches the
   `0` directory. Restore with `git checkout -- cases/.../freetip/`.

## Cases

`cases/openfoam/transonic/freetip/` — free-tip AR1.6 (wing at y=±0.8, box
at y=±4, all outer faces farfield, **no symmetry**). This is the
config-matched reference for the TSD (free tip), as opposed to the historic
campaign slabs whose tips were clamped by pressure boundaries.

Rows: naca0012 M0.80 (ny30, ny60), dfvlr_r4 M0.85/0.925/1.00,
rae2822 M0.85/0.925/1.00, plus a cascade control.
