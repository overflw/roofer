# Plan: higher-quality LoD2.2 roof polygons for German buildings (r2c)

Goal: one face per real roof plane, dormers/superstructures caught, faces fit the
points, no spurious level faces or slivers. Only roof polygons + footprints matter.

## Tracks (run in parallel)

a. **Architecture/code review** of the reconstruction chain
   (PlaneDetector -> AlphaShaper -> LineDetector -> LineRegulariser ->
   ArrangementBuilder -> ArrangementOptimiser -> ArrangementDissolver ->
   ArrangementSnapper -> ArrangementExtruder, driven by
   `apps/roofer-app/reconstruct_building.hpp`). Output: weak spots tied to r2c
   symptoms (over-segmentation, missed dormers, wrong planes on large roofs,
   eave kinks, slivers, level faces), each with file:line and a fix idea.
b. **Evaluation tooling**: `tools/score_lod22.py`, a reproducible scorer for
   roofer CityJSONSeq output vs. the input LAS/LAZ implementing r2c metrics
   (rms_in, outliers, missed_m2, bad_m2, splits/crease_m2, usable_m2 @ OFF 0.15,
   face counts, faces w/o points, level faces, slivers) with even/odd hold-out.
c. **Extra data**: orthophoto RGB/NIR (NDVI vegetation mask, edge/ridge priors),
   multi-return/intensity gating, ALKIS roof attributes. Feasibility inside
   roofer vs. r2c pre-processing; cheap prototype (point filter by returns).
d. **Literature/upstream scan**: open 3DBAG/roofer issues/PRs, City3D, PolyFit,
   KSR, learned primitive fitting -> ideas that fit roofer.

## Build & evaluate
- Build roofer with Nix dev shell (CGAL_EIGEN3_ENABLED checked), run ctest.
- Test data: data.3dgi.xyz is blocked from the sandbox (403); try public
  German ALS (Berlin/NRW open data) + ALKIS/OSM footprints; otherwise synthesize
  roofs with known planes (gable+dormer, hip, large flat-with-superstructures,
  noise levels of ALS and bDOM) so ground truth is known.
- Baseline (r2c base config) vs. change; report metrics; note what couldn't be
  verified.

## Implementation candidates (pick after a-d)
- Plane merging inside roofer (coplanar adjacent regions after region growing),
  replacing r2c's post merge.
- Robust plane refit / outlier-aware plane choice on large roofs.
- Eave-kink fix in regularisation/arrangement (extend hip lines to footprint
  corners).
- Sliver / tiny level face removal.
- Optional per-point mask attribute (vegetation) honoured by roofer.

One focused draft PR per change against overflw/roofer main; REPORT.md on this
branch as a draft PR at the end.
