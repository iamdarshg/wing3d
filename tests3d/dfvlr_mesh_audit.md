# dfvlr_r4 campaign mesh audit — what each case directory actually contains

**Date:** 2026-10-20 · **Container:** `oftrans` (started for this audit, left running)
**Scope:** all 15 `dfvlr_r4_M*_a1.25` cases in `/caseCampaign/campaign/`
**Reproduce:** `tests3d/dfvlr_mesh_audit.sh` (copied into the container and run as
`bash /tmp/dfvlr_mesh_audit.sh`) → raw output committed as
`tests3d/dfvlr_mesh_audit_raw.txt`.

The question is narrower than "which mesh produced the coefficients":

1. what mesh sits in the case directory **now** (this table), and
2. what mesh produced the **delivered** coefficient (second table) — because for
   this row the two answers are different.

---

## 1. Directory truth — mesh contents of each campaign case

Authoritative source: the `note "nPoints:… nCells:… nFaces:…"` line in
`constant/polyMesh/owner` (OpenFOAM writes it at mesh-construction time), cross-checked
against the actual point count in `constant/polyMesh/points`.

| M | nCells (owner note) | nPoints (note → actual file) | `0/polyMesh` present? | `constant/polyMesh/points` md5 | points mtime | verdict (dir contents) |
|---|---|---|---|---|---|---|
| 0.750 | 86,145 | 110,001 → 110,001 | **no** | `e1fa0119701082004979a363708b6267` | 2026-10-15 17:29 | **86k** |
| 0.775 | 86,145 | 110,001 → 110,001 | **no** | `e1fa0119701082004979a363708b6267` | 2026-10-15 22:38 | **86k** |
| 0.800 | 86,145 | 110,001 → 110,001 | **no** | `e1fa0119701082004979a363708b6267` | 2026-10-17 16:45 | **86k** |
| 0.825 | 86,145 | 110,001 → 110,001 | **no** | `e1fa0119701082004979a363708b6267` | 2026-10-17 16:45 | **86k** |
| 0.850 | **45,606** | 59,311 → **59,360** | **no** | `b7adaff7e649fbd5951936cf0bd0dd27` | 2026-10-15 16:43 | **45.6k (corrupt points)** |
| 0.875 | **45,606** | 59,311 → **59,360** | **no** | `b7adaff7e649fbd5951936cf0bd0dd27` | 2026-10-15 16:43 | **45.6k (corrupt points)** |
| 0.900 | 86,145 | 110,001 → 110,001 | **no** | `e1fa0119701082004979a363708b6267` | 2026-10-17 16:45 | **86k** |
| 0.925 | **45,606** | 59,311 → **59,360** | **no** | `b7adaff7e649fbd5951936cf0bd0dd27` | 2026-10-15 16:43 | **45.6k (corrupt points)** |
| 0.950 | **45,606** | 59,311 → **59,360** | **no** | `b7adaff7e649fbd5951936cf0bd0dd27` | 2026-10-15 16:43 | **45.6k (corrupt points)** |
| 0.975 | **45,606** | 59,311 → **59,360** | **no** | `b7adaff7e649fbd5951936cf0bd0dd27` | 2026-10-15 16:43 | **45.6k (corrupt points)** |
| 1.000 | **45,606** | 59,311 → **59,360** | **no** | `b7adaff7e649fbd5951936cf0bd0dd27` | 2026-10-15 16:43 | **45.6k (corrupt points)** |
| 1.025 | **45,606** | 59,311 → **59,360** | **no** | `b7adaff7e649fbd5951936cf0bd0dd27` | 2026-10-15 16:43 | **45.6k (corrupt points)** |
| 1.050 | **45,606** | 59,311 → **59,360** | **no** | `b7adaff7e649fbd5951936cf0bd0dd27` | 2026-10-15 16:43 | **45.6k (corrupt points)** |
| 1.075 | **45,606** | 59,311 → **59,360** | **no** | `b7adaff7e649fbd5951936cf0bd0dd27` | 2026-10-15 16:50 | **45.6k (corrupt points)** |
| 1.100 | **45,606** | 59,311 → **59,360** | **no** | `b7adaff7e649fbd5951936cf0bd0dd27` | 2026-10-15 16:50 | **45.6k (corrupt points)** |

Notes:

- **No campaign case has a time-directory `polyMesh` at all** — `polyMesh_dirs` for every case
  is exactly `<case>/constant/polyMesh`. `0/polyMesh/points` exists only in the scratch case
  `/dfvlr_coarse` (md5 `ad7b845b4ae4da6e3cce44f936379a3a`, 59,311 pts), not in any `campaign/` case.
- The 45,606 cases are internally inconsistent: the `owner` note says 59,311 points while the
  `points` file carries **59,360** (md5 `b7adaff…`) — the 49-stray-point corruption root-caused in
  `ERROR_REDUCTION_SCOPE.md` §5. Byte-identical to `/meshCache/dfvlr_r4.old`,
  `/meshCache/dfvlr_r4.stale45606`, `/dfvlr_coarse/constant/polyMesh/points` and the
  `…_stale59k` probe mesh.
- The 86k cases are internally consistent (110,001 in both places, md5 `e1fa0119…`), byte-identical
  to `/meshCache/dfvlr_r4` and the `…_86k` probe meshes.
- Reference md5s: good coarse points (59,311) = `ad7b845b4ae4da6e3cce44f936379a3a`;
  corrupt coarse points (59,360) = `b7adaff7e649fbd5951936cf0bd0dd27`; refined points (110,001) =
  `e1fa0119701082004979a363708b6267`.

### Stale `log.check` — do not use these for provenance

| M | `log.check` mtime | points / cells it reports |
|---|---|---|
| 0.850, 0.875, 0.925, 0.950, 0.975, 1.000, 1.050, 1.075, 1.100 | 2026-10-12 10:20–10:32 | 28,684 / **21,000** |
| 0.900 | 2026-10-15 16:39 | 59,311 / **45,606** |
| 0.750–0.825 | — | no `log.check` |

Every `log.check` predates the mesh now on disk (21k-cell leftovers from the 2026-10-12 meshing
attempt; the 0.900 one is from the failed 2026-10-15 coarse retry, before that case was rebuilt at
86k on 2026-10-17). They are worthless as mesh attribution. (`/caseCampaign/mesh.log`,
2026-10-12 11:16, also records a superseded 54,510-cell dfvlr meshing attempt.)

---

## 2. Run provenance — which mesh produced the delivered coefficients

The directory contents above are **not** the mesh that produced the delivered values. The old-band
case dirs were overwritten during a failed retry, which is the origin of the contradiction.

| M | delivered value source | t_end | run mesh | proof strength |
|---|---|---|---|---|
| 0.750 | `results2/dfvlr_r4_M0.750_a1.25.dat` | 0.25 | 86k | case dir rebuilt from golden cache (`solver_df3b.sh`) |
| 0.775 | `results2/…M0.775….dat` | 0.25 | 86k | same |
| 0.800 | `results2/…M0.800….dat` | 0.25 | 86k | rebuilt 2026-10-17 (`solver_v5`) |
| 0.825 | `results2/…M0.825….dat` | 0.25 | 86k | rebuilt 2026-10-17 (`solver_v6`) |
| 0.850 | `results/…M0.850….dat` (old era) | 0.25 | **86k** | documented (`dfvlr_runner.sh` header, `dfchk.sh`); not bit-tested |
| 0.875 | `results/…M0.875….dat` (old era) | 0.25 | **86k** | documented; not bit-tested |
| 0.900 | `results2/…M0.900….dat` | 0.25 | 86k | rebuilt 2026-10-18 (`solver_v4`) |
| 0.925 | `results/…M0.925….dat` (old era) | 0.25 | **86k** | **proven bit-exact** vs 86k probe2 |
| 0.950 | `results/…M0.950….dat` (old era) | 0.2014 | **86k** | documented + unreachable on corrupt coarse mesh |
| 0.975 | `results/…M0.975….dat` (old era) | 0.25 | **86k** | documented |
| 1.000 | `results/…M1.000….dat` (old era) | 0.25 | **86k** | **proven bit-exact** vs 86k probe2 |
| 1.025 | `results/…M1.025….dat` (old era) | 0.25 | **86k** | documented |
| 1.050 | `results/…M1.050….dat` (old era) | 0.25 | **86k** | documented |
| 1.075 | `container_coeffs/…M1.075….coeff` only (not in report table) | 0.0764 | **86k** | partial run; corrupt coarse dies at t≈4e-5 |
| 1.100 | `container_coeffs/…M1.100….coeff` only (not in report table) | 0.0754 | **86k** | same |

Evidence for the run-mesh column:

1. **The old-band runner says 86k.** `campaign/container_coeffs/dfvlr_runner.sh` (committed
   2026-10-14 in `dd6e57b`) opens with: *“Dedicated DFVLR runner: runs all 11 dfvlr cases
   sequentially with the **new 86k mesh** + hardened numerics.”* Its sibling check script
   `dfchk.sh` opens with *“=== dfvlr M0.85 mesh check (should be **86145**) ===”*. `batch.log`
   records the successful old-era completions: `M0.850 done t=0.25 CL=2.674053e-01`,
   `M0.875 … CL=2.729032e-01`, `M0.925 … CL=2.779525e-01`, `M0.950 … t=0.2014`.
2. **Two of those values are proven 86k, bit-exactly.** Probe2 re-ran M0.925/M1.000 on
   `/meshCache/dfvlr_r4` (110,001 pts, `e1fa0119…`): CL = 2.779525e-01 / 2.126116e-01, and the
   final coefficient rows are byte-identical to the corpus dumps
   (`…/probes/dfvlr_r4_M0.925_86k/postProcessing/forces/0/coefficient.dat` vs
   `campaign/container_coeffs/dfvlr_r4_M0.925_a1.25.coeff`).
   The repaired coarse mesh gives CL = 0.3626993 / 0.2529717 at the same points
   (`ERROR_REDUCTION_SCOPE.md` §5) — nowhere near the corpus values.
3. **The 45,606 mesh physically cannot produce any delivered value.** Its `constant/polyMesh/points`
   is the 59,360-point corrupt file; as documented in §5, every consumer died at
   t ≈ 4e-5. The only repaired 45.6k mesh was assembled on 2026-10-19/20 (`tmp_lane_fix.sh`).
   Yet every delivered dfvlr history reaches t ≥ 0.075 (most 0.25) — so none can be a
   45.6k result.
4. **The coarse copies are dated to the failed 2026-10-15 retry, not the solve.** Every
   45.6k case has `log.run` headed `Date: Oct 15 2026`, starting at exactly the points-file
   mtime (e.g. M0.850 `16:43:11`, M1.075 `16:50:04`), and every one ends without forces
   (`End` immediately, or an FPE stack trace). `results2/*.dat` for 0.850, 0.875, 0.925–1.050
   are **0 bytes**, and 1.075/1.100 have no `.dat` at all (`claimed-v3` marker only).
   `solver_df2.sh` logs `start … (mesh=dfvlr_r4.old)` and `NOFORCES` for the 16:49–16:50 copies;
   the 16:43 copies are from the same retry window, while `/meshCache/dfvlr_r4` was still the
   coarse cache (it was replaced by the 86k copy at 16:52).
5. **The later rebuilds are dated correctly.** 0.750/0.775 points mtimes are 2026-10-15 17:29 /
   22:38; 0.800/0.825/0.900 are 2026-10-17 16:45 — all after the 86k cache promotion, and all
   carry the `e1fa0119…` (86k) points. No old-band case was rebuilt.

---

## 3. Raw evidence

Command:

```bash
docker start oftrans
docker cp tests3d/dfvlr_mesh_audit.sh oftrans:/tmp/dfvlr_mesh_audit.sh
docker exec oftrans bash /tmp/dfvlr_mesh_audit.sh
```

Verbatim header note + points line, per case
(`grep -m1 note <case>/constant/polyMesh/owner`; `awk` first integer line of `points`):

```text
dfvlr_r4_M0.750_a1.25  note "nPoints:110001  nCells:86145  nFaces:281816  nInternalFaces:252127";  points=110001  md5=e1fa0119701082004979a363708b6267
dfvlr_r4_M0.775_a1.25  note "nPoints:110001  nCells:86145  nFaces:281816  nInternalFaces:252127";  points=110001  md5=e1fa0119701082004979a363708b6267
dfvlr_r4_M0.800_a1.25  note "nPoints:110001  nCells:86145  nFaces:281816  nInternalFaces:252127";  points=110001  md5=e1fa0119701082004979a363708b6267
dfvlr_r4_M0.825_a1.25  note "nPoints:110001  nCells:86145  nFaces:281816  nInternalFaces:252127";  points=110001  md5=e1fa0119701082004979a363708b6267
dfvlr_r4_M0.850_a1.25  note "nPoints:59311   nCells:45606  nFaces:150169  nInternalFaces:129548";  points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
dfvlr_r4_M0.875_a1.25  note "nPoints:59311   nCells:45606  nFaces:150169  nInternalFaces:129548";  points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
dfvlr_r4_M0.900_a1.25  note "nPoints:110001  nCells:86145  nFaces:281816  nInternalFaces:252127";  points=110001  md5=e1fa0119701082004979a363708b6267
dfvlr_r4_M0.925_a1.25  note "nPoints:59311   nCells:45606  nFaces:150169  nInternalFaces:129548";  points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
dfvlr_r4_M0.950_a1.25  note "nPoints:59311   nCells:45606  nFaces:150169  nInternalFaces:129548";  points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
dfvlr_r4_M0.975_a1.25  note "nPoints:59311   nCells:45606  nFaces:150169  nInternalFaces:129548";  points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
dfvlr_r4_M1.000_a1.25  note "nPoints:59311   nCells:45606  nFaces:150169  nInternalFaces:129548";  points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
dfvlr_r4_M1.025_a1.25  note "nPoints:59311   nCells:45606  nFaces:150169  nInternalFaces:129548";  points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
dfvlr_r4_M1.050_a1.25  note "nPoints:59311   nCells:45606  nFaces:150169  nInternalFaces:129548";  points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
dfvlr_r4_M1.075_a1.25  note "nPoints:59311   nCells:45606  nFaces:150169  nInternalFaces:129548";  points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
dfvlr_r4_M1.100_a1.25  note "nPoints:59311   nCells:45606  nFaces:150169  nInternalFaces:129548";  points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
```

Reference meshes:

```text
/meshCache/dfvlr_r4               points=110001  md5=e1fa0119701082004979a363708b6267
/meshCache/dfvlr_r4.old           points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
/meshCache/dfvlr_r4.stale45606    points=59360   md5=b7adaff7e649fbd5951936cf0bd0dd27
/dfvlr_coarse/constant/polyMesh/points  points=59360  md5=b7adaff7e649fbd5951936cf0bd0dd27   (unrunnable)
/dfvlr_coarse/0/polyMesh/points         points=59311  md5=ad7b845b4ae4da6e3cce44f936379a3a   (good points, scratch case only)
/caseCampaign/probes/dfvlr_r4_M0.925_coarse45k/constant/polyMesh/points  md5=ad7b845b4ae4da6e3cce44f936379a3a  (repaired 45.6k, 2026-10-19)
/caseCampaign/probes/dfvlr_r4_M0.925_86k/constant/polyMesh/points        md5=e1fa0119701082004979a363708b6267
```

Bit-exact reproduction (final coefficient row, all 12 columns):

```text
probe  dfvlr_r4_M1.000_86k : 0.25  9.980133e-02 4.995648e-02 4.984485e-02 2.126116e-01 … 5.712804e-06 -6.530434e-06
corpus M1.000_a1.25.coeff  : 0.25  9.980133e-02 4.995648e-02 4.984485e-02 2.126116e-01 … 5.712804e-06 -6.530434e-06
```

(CL is column 5 = 2.126116e-01; probe2 log: `M0.925 DONE t=0.25 CL=2.779525e-01`,
`M1.000 DONE t=0.25 CL=2.126116e-01`.)

Failed-retry provenance:

```text
case-dir `log.run` start time vs points-file mtime (each run starts the same second as its mesh copy):
  M0.850 16:43:11 / 16:43:11   M0.925 16:43:15 / 16:43:15   M1.000 16:43:17 / 16:43:17   M1.075 16:50:04 / 16:50:04
  M0.875 16:43:14 / 16:43:14   M0.950 16:43:16 / 16:43:16   M1.025 16:43:18 / 16:43:18   M1.100 16:50:07 / 16:50:07
  M0.975 16:43:17 / 16:43:16   M1.050 16:43:19 / 16:43:18

/caseCampaign/campaign/dfvlr_r4_M0.850_a1.25/log.run : Date: Oct 15 2026 / Time: 16:43:11 → "End" (no forces)
/caseCampaign/campaign/dfvlr_r4_M1.075_a1.25/log.run: Date: Oct 15 2026 / Time: 16:50:04 → FPE stack trace
results2/dfvlr_r4_M0.850_a1.25.dat = 0 bytes ; …M0.925….dat = 0 bytes ; …M1.050….dat = 0 bytes
batchv2.log: 16:50:04 start dfvlr_r4_M1.075_a1.25 (mesh=dfvlr_r4.old) → NOFORCES
container_scripts/solver_df2.sh: CPATH="/meshCache/$g.old" ; rm -rf "$c/constant/polyMesh" ; cp -r "$CPATH" …
```

---

## 4. Conclusion

**The case directories and the delivered coefficients disagree for the dfvlr old band, and both
docs read the directories as if they were run provenance.** The directory audit is unambiguous:
0.750–0.825 and 0.900 contain the **86,145-cell** mesh; 0.850, 0.875 and 0.925–1.100 contain the
**45,606-cell topology with the corrupt 59,360-point file** and no time-directory `polyMesh`. But
those 45.6k copies are dated 2026-10-15 16:43–16:50 and are accompanied by `log.run` files from that
same minute that end with no forces, zero-byte `results2/*.dat`, and `batchv2.log` entries
`(mesh=dfvlr_r4.old) → NOFORCES`. They are the debris of a failed retry, not the solve. The delivered
0.850–1.050 values were produced on 2026-10-13/14/15 by `dfvlr_runner.sh` — whose header states it
runs all 11 dfvlr cases **“with the new 86k mesh”**, and whose companion `dfchk.sh` expects
**“should be 86145”** at M0.85. M0.925 and M1.000 are proven 86k **bit-exactly** against the probe2
re-run; the repaired 45.6k mesh gives CL 0.3626993/0.2529717 versus the corpus 0.2779525/0.2126116.
And the corrupt 45.6k mesh cannot run past t≈4e-5, whereas every delivered history reaches
t ≥ 0.075 — so **no delivered dfvlr coefficient can have come from the 45,606-cell mesh. The dfvlr
row is uniform at 86,145 cells; there is no 45.6k/86k seam.**

### Exactly what is wrong

**`CAMPAIGN_REPORT.md` §9, lines 250–256** — the dfvlr clause of the mesh-resolution bullet:

> “dfvlr_r4 M0.75–0.825 and 0.90 ran on 86,145 cells while 0.85/0.875/0.925–1.10 ran on 45,606”

The entire second half is wrong (0.85/0.875/0.925–1.10 did **not** run on 45,606), as is the
methodology parenthetical on line 250 (“per-case polyMesh notes”), because those notes describe
post-hoc over-written directories. Correct replacement sentence:

> **“dfvlr_r4 ran on the 86,145-cell mesh at every M (0.750–1.100): the fills/replacement cases
> (0.750–0.825, 0.900) were rebuilt from the golden 86k cache on 2026-10-15/17/18, and the original
> 0.850–1.100 band was run on the same 86k mesh by `dfvlr_runner.sh` on 2026-10-13/14/15 — proven
> bit-exactly at M0.925 and M1.000 against an 86k re-run. The 45,606-cell mesh now sitting in the
> 0.850/0.875/0.925–1.100 case directories is a stale copy from a failed 2026-10-15 retry
> (`solver_df2.sh`, `meshCache/dfvlr_r4.old`) whose `constant/polyMesh/points` carries 49 stray
> points (59,360 vs 59,311); those directories produced no forces, so no delivered dfvlr value comes
> from the 45.6k mesh. dfvlr is therefore uniform-resolution, not a two-mesh row.”**

**`CAMPAIGN_REPORT.md` §9, lines 257–260** — “the correctly snapped points live in
`0/polyMesh/points`”: true only of the `/dfvlr_coarse` scratch case. **No campaign case contains a
`0/polyMesh`**; the good 59,311-point file exists solely at `/dfvlr_coarse/0/polyMesh/points`
(and, post-fix, in the `…_coarse45k` probe). Reword to name `/dfvlr_coarse/0/polyMesh/points`.

**`ERROR_REDUCTION_SCOPE.md` §3, table row on line 172** — three cells are wrong:

> `| dfvlr_r4 | **86,145** (fills) | **45,606** | **45,606** | **86,145** (new) | **45,606** |`

The `0.850`, `0.875` and `0.925–1.10` columns must all read **86,145** (the row becomes
86,145 everywhere). The “per-case polyMesh nCells (authoritative)” framing on line 165 is the root
error: authoritative for what is on disk today, not for what was solved.

**`ERROR_REDUCTION_SCOPE.md` §3, bullet on lines 177–183** — “the old band (0.85, 0.875,
0.925–1.10) on 45,606 — a 1.9× resolution seam mid-row” is wrong: there is no dfvlr resolution seam
and no 1.9× transition. The 2026-10-19 caveat in that bullet is directionally right but mis-frames
the fix as “the dfvlr runner re-ran all 11 cases on the 86k mesh” — no such successful re-run
happened; the old band was **originally** 86k, and the 45.6k directories are the failed-retry
debris. Correct wording: *“dfvlr_r4 is uniform at 86,145 cells for all M; the 45,606 figures
previously tabulated here were read off case directories overwritten by a failed 2026-10-15
retry, not off the meshes that produced the results.”*

**`ERROR_REDUCTION_SCOPE.md` §0, L2 row (line 36)** — the parenthetical “dfvlr (45.6k vs 86k)” and
the “smooth through both transitions” evidence are moot: dfvlr has one mesh, so it contributes no
seam. The “≤ ~1%; closed as immaterial” verdict still stands, trivially.

What is **not** affected: the rae2822 seam attribution (48,114 ↔ 53,946) is a separate row and was
not audited here; §3’s overall “mesh heterogeneity ≤ ~1%” conclusion remains valid. `CAMPAIGN_REPORT.md`
§4.1’s “dfvlr 86,145” is correct. The L5 result itself is unaffected — but note that the “coarse”
0.3627/0.2530 numbers are a **reconstructed** 45.6k mesh (2026-10-19/20), not a corpus run, so L5
measures mesh sensitivity of a mesh the corpus never used for these cells.
