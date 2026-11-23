# Report: roofer for better LoD2.2 roofs of German buildings (r2c)

Date: 2026-11-17. Fork `overflw/roofer`; base is `main` at `5a6c794` (includes the Eigen fix from #193/#194).

## Update 2026-11-18: measured in r2c (38 German tiles: 30 ALS, 8 bDOM)

r2c ran its production pipeline, scored on hold-out points against r2c main. The sandbox results below did **not** carry over:

- **#2 data-cost cap: rejected.**
  - The original version also capped the ground label. A roof face whose best plane missed it by more than the cap then cost as much as ground and turned into a hole in the roof. This was fixed in e9c915f (only roof planes are capped).
  - Cap 1.0 with ground uncapped: usable +888 m², but the 95 % CI [−381, 2092] spans 0; bad_m2 −683, off_m2 +751, missed_m2 +168, 14 buildings lose roof area.
  - Holes still occur by another path (e.g. bb `DEBBAL660002yYzuBL`, one gable slope lost).
  - Losses where one plane leans across two (sawtooth roofs, hip ends) and where small dormer or ledge faces are dropped.
- **#1 coplanar merge: rejected.**
  - Usable −1600 m² [−3553, −158], off_m2 +1194; sawtooth roofs −826 m².
  - The RMS-ratio rule joins faces that r2c's per-face rule keeps apart. That rule requires each face to keep its own fit on the merged plane.
  - It also segfaulted on buildings with too few points for region growing. Fixed in 722fdef.
  - Worth trying instead: a merge using r2c's per-face rule that also snaps the neighbours onto the merged plane.
- **#3 and #4 (bug fixes): merged.**

**Lesson.** The US QL1 test areas and roofprint outlines did not predict the German result. Neither the sawtooth/hip-end failure nor the roof holes showed up there. The "Recommended r2c trial" below is superseded. The graph-cut robustness item (ranked item 3) remains open.

## TL;DR

- roofer builds from source with the Nix dev shell. Eigen is linked, so CGAL does not fall back to Diagonalize_traits. 147 of 153 tests pass; the other 6 need test data from `data.3dgi.xyz`, which the sandbox proxy blocks.
- No German data could be downloaded: every German open-data host returns 403 from the sandbox. The evaluation therefore uses **USGS 3DEP QL1 ALS** (Champaign IL, 18–24 pts/m², class 6) with outlines derived from class-6 points. That gives:
  - **campus**: 12 large blocks, 42,400 m²;
  - **Urbana**: 91 houses, 13,300 m²;
  - **bDOM-like**: Urbana with synthetic image-matching noise, so models built from noisy input can be scored against the clean cloud as a real reference.
- 4 draft PRs:

| PR | Change | Default | Evidence |
|---|---|---|---|
| [#2](https://github.com/overflw/roofer/pull/2) | `arrangement-optimiser.max-data-distance`: truncated L1 data cost in the graph cut | off | **Strongest result.** Campus: bad_m2 1705 → 352 (cap 1.0) / 278 (cap 0.5); slivers 65 → 35/14; bad faces 50 → 26/12; level faces 303 → 266/224. Cap 1.0 keeps missed_m2 at baseline; cap 0.5 nearly triples it (43 → 119). |
| [#1](https://github.com/overflw/roofer/pull/1) | `plane-detector.merge-coplanar`: merge adjacent regions that fit one plane | off | Splits −19 to −35 % on campus and bDOM-like (crease_m2 1583 → 1185, 546 → 400). Against the clean reference on bDOM-like: usable +23 m², bad −75 m². Neutral on clean ALS houses. |
| [#3](https://github.com/overflw/roofer/pull/3) | Angle-clustering wrap-around at π (bug) | fix | Near-N-S lines were never clustered. Small, mixed effect on the metrics; this is a correctness fix. |
| [#4](https://github.com/overflw/roofer/pull/4) | `snap-threshold` was compared against a squared distance (bug) | fix | No effect at the default of 1 m. Matters for tuned values and for `unit-scale`. |

- **Recommended r2c trial:**

  ```toml
  [reconstruction.arrangement-optimiser]
  max-data-distance = 1.0

  [reconstruction.plane-detector]
  merge-coplanar = true
  ```

  On campus this gives bad_m2 1705 → 282, splits 69 → 46, faces 950 → 856 and missed_m2 43 → 65, with usable_m2 unchanged. Try it in both `coarse` and `fine`. It may let `fine` work on large roofs without r2c's post-merge.

- **Finding:** the graph cut is fragile. A millimetre-level refit of the planes flipped whole faces: a pyramid facet on one house, a 1,223 m² face on a campus block. This matches r2c's "a large main face flips to a wrong plane". The cap in #2 makes the large-face case robust. The root cause (see ranked item 3) is still open.

## 1. What was done

1. A plan (`PLAN.md`) and four parallel research tracks: code review, scorer, extra data, literature/upstream. Their findings are summarised in sections 4–6.
2. Build: `nix develop` failed at first because the LAStools source tarball from GitHub got a 403. Fixed locally by cloning the tag over git and `nix-store --add-fixed` (same store path). Tests are configured with `FETCHCONTENT_SOURCE_DIR_*` pointing at empty directories, because the test data is blocked.
3. Evaluation tooling on this branch:
   - `tools/score_lod22.py`: r2c metrics (rms_in, outliers, missed_m2, bad_m2, splits/crease_m2, usable_m2 at OFF 0.15, face counts, level faces, slivers, faces without points), with even/odd hold-out, a `compare` mode and pytest self-tests.
   - `tools/eval/`: EPT fetch, outline derivation, bDOM-like noise, runner, table script, plan-view renders, and an inner-wall counter (eave-kink proxy).
   - `tools/extradata/`: prototypes for vegetation filtering (multi-return, NDVI, planarity guard) and DOP NDVI sampling.
4. Each variant was run on the three test sets and scored, and renders were checked: the top movers, plus b10, b11, b64 and b72.

## 2. Results

r2c base config (plane-epsilon 0.2, min-plane-points 25, extension 10, complexity-factor 0.95, clip-terrain false). Scores are on hold-out points.

**Campus: 12 buildings**

| metric | base | merge | cap 0.5 | cap 1.0 | merge + cap 1.0 |
|---|---:|---:|---:|---:|---:|
| faces | 950 | 935 | 746 | 872 | 856 |
| level faces | 303 | 289 | 224 | 266 | 256 |
| slivers | 65 | 64 | 14 | 35 | 34 |
| faces without points | 7 | 9 | 1 | 2 | 2 |
| bad faces | 50 | 49 | 12 | 26 | 27 |
| rms_in p90 | 0.163 | 0.163 | 0.140 | 0.152 | 0.153 |
| bad_m2 | 1705 | 1628 | 278 | 352 | 282 |
| missed_m2 | 43 | 62 | 119 | 47 | 65 |
| splits / crease_m2 | 69 / 1583 | 56 / 1185 | 44 / 1217 | 55 / 1282 | 46 / 993 |
| usable_m2 | 23,360 | 23,360 | 23,500 | 23,450 | 23,430 |
| off_m2 | 2140 | 2206 | 2538 | 2216 | 2274 |

**Urbana: 91 houses**

All variants are within noise: usable 7,535–7,553 m², bad_m2 66–160, splits 5–8. Only 4 of 629 planes get merged. A complexity-factor of 0.90 vs 0.98 shows the known trade-off (fewer splits but more off_m2, and the reverse), with or without the cap.

**bDOM-like, scored against clean ALS**

| | base | merge | merge + cap 1.0 |
|---|---:|---:|---:|
| usable_m2 | 4,998 | 5,021 | 5,025 |
| bad_m2 | 5,167 | 5,092 | 5,118 |
| splits | 228 | 209 | 203 |

The self-referential score of the same models shows usable −31 m² for merge. Split faces follow noise blobs, so the per-face metric rewards them. This is a concrete case of the "metrics are self-referential" warning.

**Vegetation gate** (`tools/extradata/filter_roof_veg.py`: multi-return non-last and not locally planar). It removes 2.9 % (Urbana) and 2.3 % (campus) of class-6 points. Urbana, scored against the unfiltered cloud: faces 585 → 575, bad_m2 138 → 74, usable 7,541 → 7,521. Usable drops because the removed points still count as off.

**Runtime:** no measurable change from any option. The campus run (12 large blocks) takes 67–136 s depending on the number of threads.

## 3. Ranked improvement list

Impact is for r2c's goals; effort is for implementing it in roofer.

| # | Improvement | Impact | Effort | Risk | Status |
|---|---|---|---|---|---|
| 1 | Truncated data cost in the graph cut (#2) | high (large roofs, slivers, level faces) | done | medium: changes what complexity-factor means; cap 0.5 loses superstructures | draft PR, needs r2c ablation |
| 2 | Adjacency-aware coplanar merge (#1) | medium-high (splits, bDOM) | done | low (off by default; step and angle guards) | draft PR |
| 3 | **Graph-cut robustness for no-data and big faces.** Faces with no heightfield cells cost 0 for every plane (`rasterise_polygon(...,false)`). Holes are ignored (the outer ring is rasterised only). Labels start at 0 (`preset_labels=false`). Tiny plane changes flip whole faces. Fix: subtract the inner rings; give no-data faces a nearest-neighbour fill or the neighbours' cost; preset labels to the plane covering each face's centroid. | high (wrong planes) | 1–2 days | medium | not started; reproduction in the #1 description (Urbana b10) |
| 4 | **Eave kinks.** Hip and ridge lines are clipped to the point overlap (`PlaneIntersector.cpp:100-109`). Gaps are bridged by `chain()` segments, which then go into the regulariser as real lines and get extended by `extension`. The footprint never reaches the regulariser. Fix (a): pass the footprint to `LineRegulariser::compute` and snap the end of a priority-2 line to a footprint vertex within about 1 m of its supporting line. Fix (b): flag bridge segments in `chain()` and do not extend them. | medium-high on hip roofs | ~1 day | low | not reproduced here: hips meet the corners when the outlines are derived from roof points. Likely tied to ALKIS ground footprints and overhangs, so it needs r2c data. |
| 5 | Sliver / small-face relabel after the graph cut. A face with area < A_min or thinness < w_min takes the neighbour plane with the lowest step (`arr_dissolve_step_edges_naive` already computes steps). | medium | ~0.5 day | low | the cap already removes half of the slivers |
| 6 | **Multi-class and multi-return input in `StreamCropper`.** `building_classes = [6, 20]`, plus `filter-multireturn` and `ndvi-max` read from LAS returns/NIR at `StreamCropper.cpp:794`. Removes a LAZ pre-pass in r2c and fixes Berlin roofs filed under vegetation classes. | medium (vegetation) | ~0.5–1 day | low | prototype as pre-processing in `tools/extradata` |
| 7 | Plane-detector normals for bDOM: `normal-neighbour-count` 10–15. Region growing that tests distance only once a region is large. A second dormer pass on the residuals of large regions. | medium (bDOM shattering, dormers) | config: none; code: 1–2 days | low-medium | not tested |
| 8 | GoCoPP-style plane-configuration refinement (merge/split/add/transfer; the maintainer has a fork). Subsumes #2 and part of #7. | high | weeks | medium | research |
| 9 | Orthophoto edges as candidate lines in the arrangement (upstream issue 3DBAG/roofer#154; Gao et al. 2023 roofline detection) | high for bDOM and dormers | weeks | medium | research |
| 10 | `clip_roof_faces_to_terrain` runs even with `clip-terrain=false` (`ArrangementDissolver.cpp:128`). With r2c's constant `h_ground` it only bites when a face has a wrong plane that dips below ground. | low | trivial | low | not done |

## 4. Code review: key findings (file:line on `main`)

- **Graph cut** (`ArrangementOptimiser.cpp:188-256, 292`): uncapped L1 volume vs `(1−cf)·length`. Scored on the max-of-alpha-triangles heightfield. No-data cells are skipped; holes are ignored; labels start at 0. → #2 and item 3.
- **Region growing** (`RegionGrower.hpp:61-96`, `PlaneDetectorBase.hpp:123-130`): greedy growth with a refit every `refit-interval` inliers and no final refit. `adjacencies` is collected but never used. The only merge is CGAL `regularize_planes` (global, 25°, needs exact equality, and prints to `std::cout`). → #1.
- **Normals**: k=5 PCA plus a 41° gate (`normal-angle-threshold` default acos 0.75). Noisy on bDOM. The gate can also let a 15–25° dormer be absorbed by a 45° main plane. → item 7.
- **Line regularisation**: angle wrap bug → #3. `snap_threshold` compared against a squared distance → #4. Bridge segments from `chain()` are extended like real lines → item 4.
- **No min-area / sliver filter** anywhere. `arr_filter_biggest_face` is unused. → item 5.

## 5. Extra data: assessment

Attribute flow:
- Only XYZ, classification and GPS time (for the acquisition year) reach reconstruction (`StreamCropper.cpp:771-801`, `add_point` :374-443).
- Returns, intensity, RGB and NIR are dropped.
- Only one building class is accepted (`bld-class`).

| Source | Value | Where | Risks |
|---|---|---|---|
| Multi-return gate | high for vegetation in ALS (measured in r2c: 88–99 % of off-building vegetation vs 5–13 % of roof points) | r2c pre-processing now; roofer `StreamCropper` later (item 6) | roof edges and glass are also multi-return → combine with a planarity guard (as in the prototype) |
| Per-point NIR (ST bDOM) | high | pre-processing, or item 6 | green roofs (high NDVI), shadows |
| DOP RGBI NDVI raster | high where there is no per-point NIR | r2c pre-processing (`sample_dop_ndvi.py`, PDAL `filters.colorization`) | lean in non-true orthos; date mismatch with leaf-off ALS; shadows. In a synthetic test with 1 m lean, 13 % of a gable roof was removed until ridges were protected by a two-plane / normal-consistency check |
| DOP edges / ridges | medium-high for bDOM | new line input to the arrangement (item 9) | shadow and texture edges |
| ALKIS `dachform` / official LoD2 `roofType` | medium, as a prior and plausibility check | r2c post-check via `id-attribute`; re-run with other parameters on a mismatch | sparse and outdated; never use as a hard constraint |

No DOP endpoint (NI STAC, NRW, BB, BE, RP) was reachable from the sandbox.

## 6. Upstream / literature (short)

**Upstream issues and PRs:**
- 3DBAG/roofer#172 (open): planes extend beyond the data. The maintainer blames missing boundary lines. Relates to items 3 and 4.
- #154 (open): external break lines. The hook for orthophoto lines.
- #169: `unit-scale` for sparse data.
- #103: low-resolution DSM.
- #193/#194: the Eigen fix (already in this base).
- No upstream work on plane merging, dormers or vegetation.

**Literature, transferable ideas by gain per effort:**
1. Adjacency-aware merge of planes by fit (Khoshelham 2005; Yan, Shan & Jiang 2014). Implemented in #1.
2. GoCoPP (Yu & Lafarge 2022).
3. Orthophoto roofline priors (Gao, Peters & Stoter 2023).
4. NDVI and multi-return vegetation removal.
5. Graph-cut label costs and compactness, City3D-style.
6. Kinetic shape reconstruction or PolyFit, which would replace the core (too big a change).
7. Learned wireframes (PBWR, Point2Roof) as line candidates. RoofDiffusion as a possible bDOM denoiser. Both are research items.

## 7. What was not verified (needs r2c's local pipeline)

- **German data:**
  - ALS of Berlin, NRW and NI;
  - real bDOM (20 cm image matching), here replaced by a synthetic noise model;
  - **ALKIS ground footprints**: the outlines used here are roofprints, so overhang effects and eave kinks (item 4) are not exercised.
- r2c's full chain: stamp → roofer → CityJSONSeq → 3DCityDB, the `coarse`/`fine`/`refine` choice and the post-merge. Run `ablate.py` on the 24 fresh tiles with:
  - `max-data-distance` ∈ {0.5, 1.0, 1.5};
  - `merge-coplanar` on/off;
  - each of these in the `coarse` and `fine` configs.
- An eye check against r2c's 135 judged buildings. My renders cover about a dozen buildings only.
- The dormer trade-off of the cap on German small houses: missed_m2 with cap 0.5.
- The 6 data-dependent upstream tests (wippolder, issue-64, issue-71, underpasses).

## 8. Next steps for orthophoto integration (ordered)

1. **Now, in r2c pre-processing.** Add a vegetation gate before stamping: multi-return plus a planarity guard on all ALS, and per-point NIR NDVI on ST bDOM. Use `tools/extradata/filter_roof_veg.py` or the PDAL pipeline. Measure on Berlin.
2. **DOP NDVI where no per-point NIR exists.** Prefer TrueDOP (BE) and true-ortho DOPs (NI, NRW). Use a VRT over the COG/JP2 tiles, a median-smoothed NDVI and a shadow mask. Skip tiles where the DOP and ALS dates are more than 2 years apart. Protect ridges and dormer edges with a two-plane or normal-consistency test.
3. **Berlin class handling.** Add multi-class building input and the return/NIR gates in `StreamCropper` (item 6). That avoids a LAZ re-write per tile.
4. **Orthophoto line priors.** Detect rooflines on TrueDOP (e.g. Gao et al. 2023) and add them as an optional extra segment input in `ArrangementBuilder`. Let the graph cut choose whether to use them (upstream #154). This is the most promising route for bDOM, whose alpha-shape boundaries are poor.
5. Use ALKIS `dachform` or the official LoD2 `roofType` as a post-check, with re-runs for mismatches.
