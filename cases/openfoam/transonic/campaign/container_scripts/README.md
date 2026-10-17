# OpenFOAM Campaign Lane Architecture (container `oftrans`)

Snapshot of the container's `/tmp` worker scripts that produced the 52-case
OF corpus (original 36-case programme + the M0.75-0.825 fills + the
dfvlr_r4 M0.900 replacement). All lanes: `rhoCentralFoam`, `localEuler`
(LTS), `deltaT 2e-5` (`1e-5` on the dfvlr_r4 cases), `maxCo 0.5`,
`maxDeltaT 1e-4`, `endTime 0.25`,
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

## dfvlr_r4 controlDict caveat (2026-10-17)

The dfvlr_r4 case family stores several entries per line
(`startFrom startTime; startTime 0; stopAt endTime; endTime 0.25;`).
Line-oriented `sed` edits that rewrite entire lines delete the trailing
entries silently — this broke the first v4/v5/v6 attempt
(`startTime`/`stopAt`/`endTime` vanished => FATAL IO at startup). The
dfvlr_r4 lanes therefore run the case controlDicts as generated:
`solver_v4/v5/v6.sh` carry no `sed` (rewritten versions kept here). The
M0.900 case's `0/U` had additionally been seeded from the old M0.85 copy;
it was regenerated to the correct freestream `(317.677 0 6.932)`
(U = 0.9 x 353.059 m/s, alpha 1.25 deg) before the 2026-10-17 relaunch.

## Case directory layout

Each `cases/openfoam/transonic/campaign/<geo>_M<mach>_a1.25/` holds the
case dicts (`system/`), initial fields (`0/`), and the triSurface STL;
`constant/polyMesh` is copied per run from `/meshCache/<geo>` (golden
snappyHexMesh caches) or `<geo>.old` (biconvex/dfvlr proven caches).
Logs: `log.run` (per case), `batch*.log` (per lane).
