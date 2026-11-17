# Evaluation helpers (r2c-style)

Scripts used to compare roofer variants. `../score_lod22.py` is the scorer.

- `ept_fetch.py DATASET lon0 lat0 lon1 lat1 out.laz EPSG`: fetch a bbox from a USGS 3DEP EPT dataset (public S3).
- `make_footprints.py in.laz out.gpkg EPSG`: outlines from class-6 points (stand-in for ALKIS; these are roofprints).
- `bdom_noise.py in.laz out.laz`: add image-matching-like noise (smooth 2-4 m blobs, 0.22 m std, plus 0.05 m white). Score the result against the clean cloud for a true reference.
- `run.sh BIN NAME AREA [extra.toml]`: run `$EVAL_DIR/bin/BIN` with `base.toml` (r2c base config) plus the overrides in `extra.toml`.
- `score_all.sh AREA run...`: score the runs and print a table.
- `render.py out.png pc.laz BID rundir...`: plan view per run: faces coloured by aspect/slope, plus point residuals.
- `render_walls.py`, `innerwalls.py`: inner (step) walls, a proxy for eave kinks.

Layout under `$EVAL_DIR`: `data/AREA.laz`, `data/AREA_fp.gpkg`, `bin/<binaries>`, `runs/`, `venv/` (pip install -r ../requirements-score.txt matplotlib).
