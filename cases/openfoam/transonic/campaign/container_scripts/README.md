# OpenFOAM Campaign Lane Architecture (container `oftrans`)

Snapshot of the container's `/tmp` worker scripts that produced the 52-case
OF corpus (original 36-case programme + the M0.75-0.825 fills + the
dfvlr_r4 M0.900 replacement). All lanes: `rhoCentralFoam`, `localEuler`
(LTS), `deltaT 2e-5`, `maxCo 0.5`, `maxDeltaT 1e-4`, `endTime 0.25`,
`writeInterval 0.03`, `purgeWrite 2`, forces on the wing patch (Aref 1.6,
lRef 1.0, CofR (0.25,0,0)).

## Lanes

| script | coverage | notes |
|---|---|---|
| `solver4.sh W NP` | rae2822 all M + dfvlr_r4 all M + biconvex/naca lows, W/NP partition (ran as W0/W1) | skip: any `results2/<case>.done` |
| `solver_df2.sh` | dfvlr + biconvex lows on the `.old` mesh cache | skip: `results2/<case>.done2` |
| `solver_df3b.sh` | dfvlr_r4 {0.750..0.825, 0.900} on the golden 86k mesh | skip: `.done2` or `.done`=="done" |
| `solver5.sh` | naca0012 {0.775, 0.800, 0.825}; claims at start | skip: `.done2` |
| `solver_v4/v5/v6.sh` | dfvlr_r4 {0.900}/{0.800}/{0.825}, one case each | launched by `launch_v4/v5/v6.sh` wrappers that wait for a lane slot to free |

## Marker protocol (results2/)

| marker content | meaning |
|---|---|
| `.done` = `done` | real completed run (`.dat` has the final coefficient row) |
| `.done` = `claimed-<lane>` | reserved for that lane (stops other lanes' arrival races) |
| `.done2` = `done` / `no-forces` | completion record for the `.done2`-checking lanes; `no-forces` = solver produced no forces (rerun separately) |
| `.done2` = `claimed-<lane>` | pre-claim so slow lanes skip a case another lane owns |
| empty `.dat` | ignored by `tests3d/final_table.py` (falls back to `container_coeffs/` corpus) |

## Case directory layout

Each `cases/openfoam/transonic/campaign/<geo>_M<mach>_a1.25/` holds the
case dicts (`system/`), initial fields (`0/`), and the triSurface STL;
`constant/polyMesh` is copied per run from `/meshCache/<geo>` (golden
snappyHexMesh caches) or `<geo>.old` (biconvex/dfvlr proven caches).
Logs: `log.run` (per case), `batch*.log` (per lane).
