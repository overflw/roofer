"""render.py OUT.png PC.laz BID RUNDIR [RUNDIR...]  -- plan view per run: faces coloured by
aspect/slope (top), points coloured by signed residual to their face (bottom)."""
import sys, colorsys, numpy as np
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import score_lod22 as S
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import shapely
out, pc, bid, runs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
blds = [{b.bid: b for b in S.read_buildings([r + "/out"], "identificatie")}.get(bid) for r in runs]
fp = next(b for b in blds if b is not None)
allp = shapely.unary_union([f.poly for b in blds if b for f in b.faces])
x0, y0, x1, y1 = allp.buffer(2).bounds
P = S.read_points([pc], [6], bbox=(x0, y0, x1, y1), use_all=True)
fig, ax = plt.subplots(2, len(runs), figsize=(5 * len(runs), 10), squeeze=False)
for j, (r, b) in enumerate(zip(runs, blds)):
    a, a2 = ax[0, j], ax[1, j]
    res = np.full(len(P.x), np.nan)
    for f in (b.faces if b else []):
        n = f.normal; asp = (np.degrees(np.arctan2(n[0], n[1])) % 360) / 360
        slope = np.degrees(np.arccos(min(1, n[2])))
        col = (0.75, 0.75, 0.75) if slope < 3 else colorsys.hsv_to_rgb(asp, 0.25 + 0.75 * min(slope, 50) / 50, 0.95)
        for g in getattr(f.poly, "geoms", [f.poly]):
            xs, ys = g.exterior.xy; a.fill(xs, ys, color=col, ec="k", lw=0.6)
            for h in g.interiors: hx, hy = h.xy; a.fill(hx, hy, color="w", ec="k", lw=0.6)
        m = shapely.contains_xy(f.poly, P.x, P.y)
        res[m] = P.z[m] - f.plane_z(P.x[m], P.y[m])
    nf = len(b.faces) if b else 0
    a.set_title(f"{r.rstrip('/').split('/')[-1]}: {nf} faces"); a.set_aspect("equal"); a.set_xlim(x0, x1); a.set_ylim(y0, y1)
    sc = a2.scatter(P.x, P.y, c=np.clip(res, -0.5, 0.5), s=2, cmap="RdBu_r", vmin=-0.5, vmax=0.5)
    nan = np.isnan(res); a2.scatter(P.x[nan], P.y[nan], s=2, c="k")
    for f in (b.faces if b else []):
        for g in getattr(f.poly, "geoms", [f.poly]):
            xs, ys = g.exterior.xy; a2.plot(xs, ys, "k-", lw=0.4)
    off = np.nanmean(np.abs(res) > 0.15)
    a2.set_title(f"residual (±0.5 m), off>0.15: {off:.0%}"); a2.set_aspect("equal"); a2.set_xlim(x0, x1); a2.set_ylim(y0, y1)
fig.suptitle(bid); fig.tight_layout(); fig.savefig(out, dpi=70)
