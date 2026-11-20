#!/usr/bin/env python3
"""Score roofer LoD2.2 roof models against the point cloud they were made from.

The scorer reads roofer CityJSONSeq output (``*.city.jsonl``: the first line is
a CityJSON header carrying ``transform``, every further line a
``CityJSONFeature`` with integer vertices) and the LAS/LAZ point cloud(s) used
for reconstruction. Every ``RoofSurface`` of every ``lod: "2.2"`` geometry
(Solid, MultiSurface, CompositeSurface, MultiSolid, CompositeSolid) of a
``Building`` or its ``BuildingPart`` children is scored against the points
whose XY lie inside the face polygon (2D projection, holes respected).

Usage
-----
    python tools/score_lod22.py [score] --city OUT.city.jsonl [more ...] \\
        --pc points.laz [more ...] [--classes 6] [--all] \\
        [--out-json s.json] [--out-csv buildings.csv] [--out-face-csv faces.csv]

    python tools/score_lod22.py compare a.json b.json

``--city`` accepts files or directories (directories are searched recursively
for ``*.city.jsonl``, ``*.terrain.city.jsonl`` is skipped; a ``metadata.json``
next to/above split features supplies the transform, as written by roofer
with ``--split-cjseq``). Install dependencies with
``pip install -r tools/requirements-score.txt``.

Point sets ("hold" mode, the default)
-------------------------------------
Points are filtered by classification (``--classes``, default 6). A point's
parity is its index within its LAS file (before class filtering), so the split
is stable and independent of the filter. *Decision* points (even index) are
only used to fit new models (the merged plane of the split test). *Score*
points (odd index) are used for every reported number. ``--all`` uses all
points for both.

Metric definitions (all areas are 2D plan areas in m2 unless named *_3d)
------------------------------------------------------------------------
Face plane: Newell normal of the face's exterior ring (oriented upwards) through
the ring's vertex mean. Residual r = z_point - z_plane(x_point, y_point)
(vertical distance). Faces with |n_z| < 0.1 (near-vertical) are counted as
``steep`` and get no point metrics.

* n_pts: score points inside the face. ``nodata`` face: n_pts == 0.
* rms_in: sqrt(mean(r^2)) over points with |r| <= --inlier-dist (0.5).
* outlier_share: share of the face's points with |r| > --inlier-dist.
* bad face: n_pts > 0 and (rms_in > --rms-max (0.2) or no inliers at all).
* missed_m2: score points with --missed-min (0.5) < r <= --missed-max (4.0)
  are rasterised on a global --cell (0.5 m) grid; enclosed empty cells are
  filled (binary_fill_holes); cells whose centre lies outside the face shrunk
  by --missed-margin (1.0 m) are dropped; 8-connected blobs with an area in
  [--blob-min, --blob-max] = [4, 30] m2 are unmodelled dormers. missed_m2 is
  their summed cell area.
* bad_m2 (per face) = face area if the face is bad, else min(missed_m2, area)
  (no double counting); summed per building / dataset.
* splits / crease_m2: two roof faces of one building are adjacent if the part
  of one boundary within --adj-dist (0.05 m) of the other boundary is longer
  than --adj-len (0.5 m). A plane z = ax + by + c is least-squares fitted to the
  decision points of both faces that lie within --inlier-dist of their own
  face plane. With rms_m the merged plane's rms_in on a face's score points,
  the pair is a split if for both faces rms_m <= max(own rms_in) * 1.1 + 0.02,
  rms_m <= --rms-max and the merged plane keeps >= 90 % of the face's own
  inlier count (guard). Each face needs >= --split-min-pts decision and score
  points. crease_m2 = summed area of the smaller face of each split pair,
  every face counted at most once.
* usable_m2 / off_m2 / empty_m2: cells of the --cell grid whose centre lies in
  the face shrunk by --usable-margin (0.3 m) are eligible (eligible_m2). A cell
  with score points is usable if the median |r| of its points <= --off (0.15),
  OFF if more than 50 % of its points have |r| > --off, empty if it has no
  points.
* level face: |n_z| > cos(--level-deg (3 deg)).
* sliver: area_3d < --sliver-area (1.0 m2) or area_3d / perimeter_3d^2 <
  --sliver-thinness (0.01; a square is 0.0625, a 1:25 rectangle ~0.0092).
* rms_in p50/p90 over faces with inliers, plain and area-weighted (2D area).
* fallback building: a Building with no LoD2.2 RoofSurface at all (its score
  points inside the lod 0 footprint are still counted). ``rf_extrusion_mode``
  values are tallied when present.

Outputs: summary table on stdout, ``--out-json`` (params, inputs, summary),
``--out-csv`` (per building), ``--out-face-csv`` (per face).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from dataclasses import dataclass, field

import numpy as np

SCHEMA_VERSION = 1

# ----------------------------------------------------------------------------
# CityJSON parsing
# ----------------------------------------------------------------------------


@dataclass
class Face:
    bid: str
    part: str
    k: int  # running surface index within the part's lod2.2 geometries
    poly: object  # shapely Polygon / MultiPolygon (2D)
    normal: np.ndarray  # unit, n_z >= 0
    centre: np.ndarray  # point on plane
    area2d: float
    area3d: float
    perim3d: float

    def plane_z(self, x, y):
        n = self.normal
        return self.centre[2] - (n[0] * (x - self.centre[0]) + n[1] * (y - self.centre[1])) / n[2]


@dataclass
class Building:
    bid: str
    attributes: dict
    footprint: object = None
    parts: list = field(default_factory=list)
    faces: list = field(default_factory=list)


def _newell(ring):
    c = ring.mean(axis=0)
    p = ring - c
    q = np.roll(p, -1, axis=0)
    n = np.array([
        np.sum((p[:, 1] - q[:, 1]) * (p[:, 2] + q[:, 2])),
        np.sum((p[:, 2] - q[:, 2]) * (p[:, 0] + q[:, 0])),
        np.sum((p[:, 0] - q[:, 0]) * (p[:, 1] + q[:, 1])),
    ])
    return n, c


def _perimeter(ring):
    return float(np.linalg.norm(ring - np.roll(ring, -1, axis=0), axis=1).sum())


def _make_poly(rings2d):
    import shapely

    if len(rings2d[0]) < 3:
        return shapely.Polygon()
    poly = shapely.Polygon(rings2d[0], [r for r in rings2d[1:] if len(r) >= 3])
    if not poly.is_valid:
        poly = shapely.make_valid(poly)
        polys = [g for g in shapely.get_parts(poly) if g.geom_type == "Polygon"]
        poly = shapely.MultiPolygon(polys) if len(polys) > 1 else (polys[0] if polys else shapely.Polygon())
    return poly


def _iter_semantic_surfaces(geom):
    """Yield (surface_boundary, semantic_type) for every surface of a geometry."""
    t = geom.get("type")
    b = geom.get("boundaries", [])
    sem = geom.get("semantics") or {}
    surfs = sem.get("surfaces") or []
    vals = sem.get("values")

    def stype(v):
        if v is None or not (0 <= v < len(surfs)):
            return None
        return surfs[v].get("type")

    def at(path):
        v = vals
        for i in path:
            if v is None or i >= len(v):
                return None
            v = v[i]
        return v

    if t in ("MultiSurface", "CompositeSurface"):
        for i, s in enumerate(b):
            yield s, stype(at([i]))
    elif t == "Solid":
        for si, shell in enumerate(b):
            for i, s in enumerate(shell):
                yield s, stype(at([si, i]))
    elif t in ("MultiSolid", "CompositeSolid"):
        for di, solid in enumerate(b):
            for si, shell in enumerate(solid):
                for i, s in enumerate(shell):
                    yield s, stype(at([di, si, i]))


def _is_lod22(geom):
    lod = geom.get("lod")
    try:
        return abs(float(lod) - 2.2) < 1e-9
    except (TypeError, ValueError):
        return False


def _footprint(obj, verts):
    import shapely

    for g in obj.get("geometry", []) or []:
        if str(g.get("lod")) in ("0", "0.0") and g.get("type") in ("MultiSurface", "CompositeSurface"):
            polys = []
            for s in g["boundaries"]:
                rings = [verts[np.asarray(r, dtype=np.int64)][:, :2] for r in s]
                p = _make_poly(rings)
                if not p.is_empty:
                    polys.append(p)
            if polys:
                return shapely.union_all(polys)
    return None


def parse_feature(feat, transform, id_attribute=None):
    """Return a list of Building records from one CityJSON(Feature) object."""
    v = np.asarray(feat.get("vertices") or np.zeros((0, 3)), dtype=float).reshape(-1, 3)
    if transform:
        v = v * np.asarray(transform["scale"], float) + np.asarray(transform["translate"], float)
    cos = feat.get("CityObjects", {})
    out = []
    for oid, obj in cos.items():
        if obj.get("parents") or not str(obj.get("type", "")).startswith("Building"):
            continue
        attrs = obj.get("attributes") or {}
        bid = oid
        if id_attribute and attrs.get(id_attribute) is not None:
            bid = str(attrs[id_attribute])
        b = Building(bid=bid, attributes=attrs, footprint=_footprint(obj, v))
        # building itself + all descendants (via "children" and "parents" links)
        todo, seen = [oid], []
        while todo:
            cid = todo.pop(0)
            if cid in seen or cid not in cos:
                continue
            seen.append(cid)
            todo.extend(cos[cid].get("children", []) or [])
            todo.extend(k for k, o in cos.items() if cid in (o.get("parents") or []))
        for cid in seen:
            cobj = cos[cid]
            has = False
            k = 0
            for g in cobj.get("geometry", []) or []:
                if not _is_lod22(g):
                    continue
                has = True
                for surf, st in _iter_semantic_surfaces(g):
                    kk = k
                    k += 1
                    if st != "RoofSurface":
                        continue
                    rings = [v[np.asarray(r, dtype=np.int64)] for r in surf if len(r) >= 3]
                    if not rings:
                        continue
                    n, c = _newell(rings[0])
                    nn = np.linalg.norm(n)
                    if nn == 0:
                        continue
                    n = n / nn
                    if n[2] < 0:
                        n = -n
                    poly = _make_poly([r[:, :2] for r in rings])
                    a2 = float(poly.area)
                    a3 = a2 / n[2] if n[2] > 1e-6 else 0.5 * nn
                    b.faces.append(Face(bid=bid, part=cid, k=kk, poly=poly, normal=n, centre=c,
                                        area2d=a2, area3d=float(a3),
                                        perim3d=sum(_perimeter(r) for r in rings)))
            if has:
                b.parts.append(cid)
        out.append(b)
    return out


def _find_metadata(path):
    d = os.path.dirname(os.path.abspath(path))
    for _ in range(6):
        m = os.path.join(d, "metadata.json")
        if os.path.isfile(m):
            return m
        nd = os.path.dirname(d)
        if nd == d:
            break
        d = nd
    return None


def _collect_city_files(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            for root, _, names in os.walk(p):
                for n in sorted(names):
                    if (n.endswith(".city.jsonl") or n.endswith(".city.json")) and not n.endswith(".terrain.city.jsonl"):
                        files.append(os.path.join(root, n))
        else:
            files.append(p)
    return sorted(files)


def read_buildings(paths, id_attribute=None):
    buildings = []
    meta_cache = {}
    for path in _collect_city_files(paths):
        transform = None
        with open(path, "r", encoding="utf-8") as fh:
            for ln, line in enumerate(fh):
                line = line.strip()
                if not line:
                    continue
                obj = json.loads(line) if path.endswith(".jsonl") else json.loads(line + fh.read())
                if obj.get("type") == "CityJSON":
                    transform = obj.get("transform")
                    if obj.get("CityObjects"):
                        buildings += parse_feature(obj, transform, id_attribute)
                    continue
                if transform is None:
                    m = _find_metadata(path)
                    if m and m not in meta_cache:
                        with open(m, "r", encoding="utf-8") as mf:
                            meta_cache[m] = json.load(mf).get("transform")
                    transform = meta_cache.get(m)
                buildings += parse_feature(obj, transform, id_attribute)
    return buildings


# ----------------------------------------------------------------------------
# Point cloud
# ----------------------------------------------------------------------------


@dataclass
class Points:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    score: np.ndarray  # bool
    decide: np.ndarray  # bool


def read_points(paths, classes, bbox=None, use_all=False, chunk_size=5_000_000):
    import laspy

    xs, ys, zs, par = [], [], [], []
    classes = np.asarray(sorted(set(classes)))
    for p in paths:
        with laspy.open(p) as fh:
            offset = 0
            for ch in fh.chunk_iterator(chunk_size):
                n = len(ch)
                cls = np.asarray(ch.classification)
                m = np.isin(cls, classes)
                x = np.asarray(ch.x)[m]
                y = np.asarray(ch.y)[m]
                z = np.asarray(ch.z)[m]
                idx = (np.nonzero(m)[0] + offset).astype(np.int64)
                if bbox is not None:
                    k = (x >= bbox[0]) & (x <= bbox[2]) & (y >= bbox[1]) & (y <= bbox[3])
                    x, y, z, idx = x[k], y[k], z[k], idx[k]
                xs.append(x); ys.append(y); zs.append(z); par.append((idx & 1).astype(np.uint8))
                offset += n
    cat = lambda a, dt: np.concatenate(a).astype(dt) if a else np.zeros(0, dt)
    x, y, z, pa = cat(xs, float), cat(ys, float), cat(zs, float), cat(par, np.uint8)
    if use_all:
        sc = np.ones(len(x), bool)
        de = sc
    else:
        sc = pa == 1
        de = pa == 0
    return Points(x, y, z, sc, de)


class GridIndex:
    """Bucket points on a coarse grid for fast bbox candidate queries."""

    def __init__(self, x, y, cs=10.0):
        self.cs = cs
        self.x, self.y = x, y
        if len(x) == 0:
            self.nx = self.ny = 0
            return
        self.x0, self.y0 = float(x.min()), float(y.min())
        kx = ((x - self.x0) // cs).astype(np.int64)
        ky = ((y - self.y0) // cs).astype(np.int64)
        self.nx, self.ny = int(kx.max()) + 1, int(ky.max()) + 1
        key = kx * self.ny + ky
        self.order = np.argsort(key, kind="stable")
        self.skeys = key[self.order]

    def query(self, b):
        if self.nx == 0:
            return np.zeros(0, np.int64)
        f = lambda v, o, nmax: min(max(int((v - o) // self.cs), 0), nmax - 1)
        kx0, kx1 = f(b[0], self.x0, self.nx), f(b[2], self.x0, self.nx)
        ky0, ky1 = f(b[1], self.y0, self.ny), f(b[3], self.y0, self.ny)
        kxs = np.arange(kx0, kx1 + 1, dtype=np.int64)
        lo = np.searchsorted(self.skeys, kxs * self.ny + ky0, "left")
        hi = np.searchsorted(self.skeys, kxs * self.ny + ky1, "right")
        idx = np.concatenate([self.order[a:c] for a, c in zip(lo, hi)]) if len(kxs) else np.zeros(0, np.int64)
        x, y = self.x[idx], self.y[idx]
        return idx[(x >= b[0]) & (x <= b[2]) & (y >= b[1]) & (y <= b[3])]


# ----------------------------------------------------------------------------
# Scoring
# ----------------------------------------------------------------------------

DEFAULTS = dict(
    classes=[6], inlier_dist=0.5, rms_max=0.2, cell=0.5,
    missed_min=0.5, missed_max=4.0, missed_margin=1.0, blob_min=4.0, blob_max=30.0,
    adj_dist=0.05, adj_len=0.5, split_factor=1.1, split_add=0.02, split_min_pts=5,
    split_keep_inliers=0.9,
    off=0.15, usable_margin=0.3, level_deg=3.0, sliver_area=1.0, sliver_thinness=0.01,
    min_nz=0.1, use_all=False,
)


def _rms_in(r, d):
    inl = np.abs(r) <= d
    n = int(inl.sum())
    return (float(np.sqrt(np.mean(r[inl] ** 2))) if n else float("nan")), n


def _face_grid(face, cell):
    b = face.poly.bounds
    ix0, iy0 = int(math.floor(b[0] / cell)), int(math.floor(b[1] / cell))
    ix1, iy1 = int(math.floor(b[2] / cell)), int(math.floor(b[3] / cell))
    return ix0, iy0, ix1 - ix0 + 1, iy1 - iy0 + 1


def _inner_mask(poly, margin, ix0, iy0, nx, ny, cell):
    import shapely

    inner = poly.buffer(-margin)
    if inner.is_empty:
        return np.zeros((nx, ny), bool)
    cx = (np.arange(ix0, ix0 + nx) + 0.5) * cell
    cy = (np.arange(iy0, iy0 + ny) + 0.5) * cell
    gx, gy = np.meshgrid(cx, cy, indexing="ij")
    shapely.prepare(inner)
    return shapely.contains_xy(inner, gx.ravel(), gy.ravel()).reshape(nx, ny)


def score_face(face, P, idx, p):
    """Point metrics of one face. idx: indices of all class-filtered points inside it."""
    from scipy import ndimage

    res = dict(n_pts=0, n_in=0, rms_in=float("nan"), outlier_share=float("nan"), bad=False,
               missed_m2=0.0, missed_blobs=0, usable_m2=0.0, off_m2=0.0, empty_m2=0.0, eligible_m2=0.0,
               steep=False)
    if face.normal[2] < p["min_nz"]:
        res["steep"] = True
        return res, None
    si = idx[P.score[idx]]
    x, y = P.x[si], P.y[si]
    r = P.z[si] - face.plane_z(x, y)
    n = len(si)
    res["n_pts"] = n
    cell = p["cell"]
    ix0, iy0, nx, ny = _face_grid(face, cell)
    elig = _inner_mask(face.poly, p["usable_margin"], ix0, iy0, nx, ny, cell)
    res["eligible_m2"] = float(elig.sum()) * cell * cell
    if n == 0:
        res["empty_m2"] = res["eligible_m2"]
        return res, None
    rms, n_in = _rms_in(r, p["inlier_dist"])
    res.update(rms_in=rms, n_in=n_in, outlier_share=1.0 - n_in / n)
    res["bad"] = bool(n_in == 0 or rms > p["rms_max"])

    ix = np.clip(np.floor(x / cell).astype(np.int64) - ix0, 0, nx - 1)
    iy = np.clip(np.floor(y / cell).astype(np.int64) - iy0, 0, ny - 1)

    # --- missed (unmodelled structures above the face)
    up = (r > p["missed_min"]) & (r <= p["missed_max"])
    if up.any():
        mask = np.zeros((nx, ny), bool)
        mask[ix[up], iy[up]] = True
        mask = ndimage.binary_fill_holes(mask)
        mask &= _inner_mask(face.poly, p["missed_margin"], ix0, iy0, nx, ny, cell)
        lab, nlab = ndimage.label(mask, structure=np.ones((3, 3), bool))
        if nlab:
            sizes = np.bincount(lab.ravel())[1:] * cell * cell
            keep = (sizes >= p["blob_min"]) & (sizes <= p["blob_max"])
            res["missed_m2"] = float(sizes[keep].sum())
            res["missed_blobs"] = int(keep.sum())

    # --- usable / off cells
    cid = ix * ny + iy
    ok = elig.ravel()[cid]
    if ok.any():
        c, a = cid[ok], np.abs(r[ok])
        o = np.lexsort((a, c))
        c, a = c[o], a[o]
        _, start, cnt = np.unique(c, return_index=True, return_counts=True)
        med = 0.5 * (a[start + (cnt - 1) // 2] + a[start + cnt // 2])
        offfrac = np.add.reduceat((a > p["off"]).astype(np.int64), start) / cnt
        ca = cell * cell
        res["usable_m2"] = float((med <= p["off"]).sum()) * ca
        res["off_m2"] = float((offfrac > 0.5).sum()) * ca
        res["empty_m2"] = res["eligible_m2"] - len(cnt) * ca
    else:
        res["empty_m2"] = res["eligible_m2"]
    return res, r


def find_splits(faces, fres, pidx, P, p):
    """Return list of (i, j, rms_i_merged, rms_j_merged) split pairs."""
    import shapely

    if len(faces) < 2:
        return []
    polys = np.array([f.poly for f in faces], dtype=object)
    tree = shapely.STRtree(polys)
    a_idx, b_idx = tree.query(shapely.buffer(polys, p["adj_dist"]), predicate="intersects")
    pairs = sorted({(int(i), int(j)) for i, j in zip(a_idx, b_idx) if i < j})
    out = []
    d = p["inlier_dist"]
    for i, j in pairs:
        fi, fj = faces[i], faces[j]
        if fres[i]["steep"] or fres[j]["steep"]:
            continue
        shared = shapely.intersection(fi.poly.boundary, shapely.buffer(fj.poly.boundary, p["adj_dist"])).length
        if shared <= p["adj_len"]:
            continue
        ri, rj = fres[i]["rms_in"], fres[j]["rms_in"]
        if not (np.isfinite(ri) and np.isfinite(rj)):
            continue
        sets, ok = [], True
        for f in (fi, fj):
            di = pidx[f_key(f)]
            di = di[P.decide[di]]
            di = di[np.abs(P.z[di] - f.plane_z(P.x[di], P.y[di])) <= d]
            if len(di) < p["split_min_pts"]:
                ok = False
            sets.append(di)
        if not ok:
            continue
        dd = np.concatenate(sets)
        mx, my = P.x[dd].mean(), P.y[dd].mean()
        A = np.c_[P.x[dd] - mx, P.y[dd] - my, np.ones(len(dd))]
        coef = np.linalg.lstsq(A, P.z[dd], rcond=None)[0]
        lim = min(max(ri, rj) * p["split_factor"] + p["split_add"], p["rms_max"])
        rms_m = []
        for k, f in ((i, fi), (j, fj)):
            si = pidx[f_key(f)]
            si = si[P.score[si]]
            if len(si) < p["split_min_pts"]:
                ok = False
                break
            r = P.z[si] - (coef[0] * (P.x[si] - mx) + coef[1] * (P.y[si] - my) + coef[2])
            rm, nin = _rms_in(r, d)
            if not (np.isfinite(rm) and rm <= lim and nin >= p["split_keep_inliers"] * fres[k]["n_in"]):
                ok = False
                break
            rms_m.append(rm)
        if ok:
            out.append((i, j, rms_m[0], rms_m[1]))
    return out


def f_key(f):
    return (f.part, f.k)


def _wquantile(v, w, q):
    if len(v) == 0:
        return float("nan")
    o = np.argsort(v, kind="stable")
    v, cw = v[o], np.cumsum(w[o])
    if cw[-1] <= 0:
        return float("nan")
    return float(v[min(np.searchsorted(cw, q * cw[-1]), len(v) - 1)])


def score(buildings, P, p):
    import shapely

    index = GridIndex(P.x, P.y)
    cos_level = math.cos(math.radians(p["level_deg"]))
    brows, frows = [], []
    for b in buildings:
        faces = b.faces
        br = dict(building_id=b.bid, n_parts=len(b.parts), n_faces=len(faces),
                  extrusion_mode=b.attributes.get("rf_extrusion_mode", ""),
                  rf_roof_planes=b.attributes.get("rf_roof_planes", ""),
                  rf_rmse_lod22=b.attributes.get("rf_rmse_lod22", ""),
                  fallback=int(len(faces) == 0))
        if not faces:
            n = 0
            if b.footprint is not None and not b.footprint.is_empty:
                c = index.query(b.footprint.bounds)
                c = c[P.score[c]]
                shapely.prepare(b.footprint)
                n = int(shapely.contains_xy(b.footprint, P.x[c], P.y[c]).sum())
            br.update(n_pts=n, roof_area_2d=0.0, roof_area_3d=0.0)
            brows.append(br)
            continue
        bounds = np.array([f.poly.bounds for f in faces])
        cand = index.query((bounds[:, 0].min(), bounds[:, 1].min(), bounds[:, 2].max(), bounds[:, 3].max()))
        cx, cy = P.x[cand], P.y[cand]
        pidx, fres, resid = {}, [], []
        for f, fb in zip(faces, bounds):
            m = (cx >= fb[0]) & (cx <= fb[2]) & (cy >= fb[1]) & (cy <= fb[3])
            sub = cand[m]
            if not f.poly.is_empty and len(sub):
                shapely.prepare(f.poly)
                sub = sub[shapely.contains_xy(f.poly, P.x[sub], P.y[sub])]
            else:
                sub = sub[:0]
            pidx[f_key(f)] = sub
            res, r = score_face(f, P, sub, p)
            fres.append(res)
            resid.append(r)
        splits = find_splits(faces, fres, pidx, P, p)
        crease_faces = set()
        partner = {}
        for i, j, _, _ in splits:
            small = i if faces[i].area2d <= faces[j].area2d else j
            crease_faces.add(small)
            partner.setdefault(i, []).append(j)
            partner.setdefault(j, []).append(i)
        all_pts = np.unique(np.concatenate([v[P.score[v]] for v in pidx.values()]))
        sq, nin = 0.0, 0
        for f, res, r, fi in zip(faces, fres, resid, range(len(faces))):
            level = bool(f.normal[2] > cos_level)
            thin = f.area3d / f.perim3d ** 2 if f.perim3d > 0 else 0.0
            sliver = bool(f.area3d < p["sliver_area"] or thin < p["sliver_thinness"])
            bad_m2 = f.area2d if res["bad"] else min(res["missed_m2"], f.area2d)
            if r is not None and res["n_in"]:
                inl = np.abs(r) <= p["inlier_dist"]
                sq += float(np.sum(r[inl] ** 2))
                nin += int(inl.sum())
            frows.append(dict(
                building_id=b.bid, part_id=f.part, face=f.k, area_2d=f.area2d, area_3d=f.area3d,
                slope_deg=math.degrees(math.acos(min(1.0, f.normal[2]))), thinness=thin,
                level=int(level), sliver=int(sliver), steep=int(res["steep"]),
                nodata=int(res["n_pts"] == 0 and not res["steep"]),
                n_pts=res["n_pts"], n_in=res["n_in"], rms_in=res["rms_in"],
                outlier_share=res["outlier_share"], bad=int(res["bad"]), bad_m2=bad_m2,
                missed_m2=res["missed_m2"], missed_blobs=res["missed_blobs"],
                usable_m2=res["usable_m2"], off_m2=res["off_m2"], empty_m2=res["empty_m2"],
                eligible_m2=res["eligible_m2"], crease=int(fi in crease_faces),
                split_with=";".join(f"{faces[q].part}:{faces[q].k}" for q in partner.get(fi, [])),
            ))
        fr = frows[-len(faces):]
        s = lambda k: float(sum(x[k] for x in fr))
        br.update(
            n_pts=int(len(all_pts)), roof_area_2d=s("area_2d"), roof_area_3d=s("area_3d"),
            n_nodata=int(s("nodata")), n_level=int(s("level")), n_sliver=int(s("sliver")),
            n_steep=int(s("steep")), n_bad=int(s("bad")),
            rms_in=(math.sqrt(sq / nin) if nin else float("nan")),
            outlier_share=(1.0 - nin / s("n_pts")) if s("n_pts") else float("nan"),
            bad_m2=s("bad_m2"), bad_rms_m2=float(sum(x["area_2d"] for x in fr if x["bad"])),
            missed_m2=s("missed_m2"), missed_blobs=int(s("missed_blobs")),
            splits=len(splits), crease_m2=float(sum(faces[q].area2d for q in crease_faces)),
            usable_m2=s("usable_m2"), off_m2=s("off_m2"), empty_m2=s("empty_m2"), eligible_m2=s("eligible_m2"),
            _sq=sq, _nin=nin,
        )
        brows.append(br)
    return summarize(brows, frows), brows, frows


def summarize(brows, frows):
    g = lambda rows, k: float(sum(r.get(k, 0) or 0 for r in rows))
    scored = [r for r in brows if not r["fallback"]]
    area = g(frows, "area_2d")
    rms = np.array([f["rms_in"] for f in frows], float)
    w = np.array([f["area_2d"] for f in frows], float)
    ok = np.isfinite(rms)
    nin = g(scored, "_nin")
    npts_faces = g(frows, "n_pts")
    modes = {}
    for r in brows:
        if r.get("extrusion_mode"):
            modes[r["extrusion_mode"]] = modes.get(r["extrusion_mode"], 0) + 1
    div = lambda a, b: a / b if b else float("nan")
    s = dict(
        buildings=len(brows), buildings_scored=len(scored), buildings_fallback=len(brows) - len(scored),
        faces=len(frows), faces_nodata=int(g(frows, "nodata")), faces_level=int(g(frows, "level")),
        faces_sliver=int(g(frows, "sliver")), faces_steep=int(g(frows, "steep")), faces_bad=int(g(frows, "bad")),
        roof_area_2d_m2=area, roof_area_3d_m2=g(frows, "area_3d"),
        points_scored=int(g(scored, "n_pts")), points_fallback_footprints=int(g([r for r in brows if r["fallback"]], "n_pts")),
        rms_in_pooled=math.sqrt(g(scored, "_sq") / nin) if nin else float("nan"),
        outlier_share=div(npts_faces - nin, npts_faces),
        rms_in_face_p50=float(np.median(rms[ok])) if ok.any() else float("nan"),
        rms_in_face_p90=float(np.quantile(rms[ok], 0.9)) if ok.any() else float("nan"),
        rms_in_face_p50_aw=_wquantile(rms[ok], w[ok], 0.5),
        rms_in_face_p90_aw=_wquantile(rms[ok], w[ok], 0.9),
        bad_rms_m2=float(sum(f["area_2d"] for f in frows if f["bad"])),
        missed_m2=g(frows, "missed_m2"), missed_blobs=int(g(frows, "missed_blobs")),
        bad_m2=g(frows, "bad_m2"), bad_frac=div(g(frows, "bad_m2"), area),
        splits=int(g(scored, "splits")), crease_m2=g(scored, "crease_m2"),
        usable_m2=g(frows, "usable_m2"), off_m2=g(frows, "off_m2"), empty_m2=g(frows, "empty_m2"),
        eligible_m2=g(frows, "eligible_m2"),
        usable_frac=div(g(frows, "usable_m2"), area), off_frac=div(g(frows, "off_m2"), area),
        usable_frac_eligible=div(g(frows, "usable_m2"), g(frows, "eligible_m2")),
    )
    for k, v in sorted(modes.items()):
        s[f"extrusion_mode_{k}"] = v
    return s


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------

LOWER_BETTER = {"buildings_fallback", "faces_nodata", "faces_sliver", "faces_bad", "rms_in_pooled",
                "outlier_share", "rms_in_face_p50", "rms_in_face_p90", "rms_in_face_p50_aw",
                "rms_in_face_p90_aw", "bad_rms_m2", "missed_m2", "missed_blobs", "bad_m2", "bad_frac",
                "splits", "crease_m2", "off_m2", "off_frac", "empty_m2"}
HIGHER_BETTER = {"usable_m2", "usable_frac", "usable_frac_eligible"}


def _fmt(v):
    if isinstance(v, float):
        if not math.isfinite(v):
            return "nan"
        return f"{v:.4f}" if abs(v) < 10 else f"{v:.1f}"
    return str(v)


def print_summary(s, out=None):
    out = out or sys.stdout
    w = max(len(k) for k in s)
    for k, v in s.items():
        print(f"{k:<{w}}  {_fmt(v):>14}", file=out)


def _write_csv(path, rows):
    rows = [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, restval="")
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (_fmt(v) if isinstance(v, float) else v) for k, v in r.items()})


def _json_safe(o):
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if isinstance(o, dict):
        return {k: _json_safe(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_json_safe(v) for v in o]
    return o


def run_score(args):
    p = dict(DEFAULTS)
    for k in DEFAULTS:
        if hasattr(args, k):
            p[k] = getattr(args, k)
    t0 = time.time()
    buildings = read_buildings(args.city, args.id_attribute)
    bb = None
    geoms = [f.poly for b in buildings for f in b.faces] + [b.footprint for b in buildings if b.footprint is not None]
    if geoms:
        bs = np.array([g.bounds for g in geoms if not g.is_empty])
        if len(bs):
            bb = (bs[:, 0].min() - 1, bs[:, 1].min() - 1, bs[:, 2].max() + 1, bs[:, 3].max() + 1)
    P = read_points(args.pc, p["classes"], bbox=bb, use_all=p["use_all"])
    summary, brows, frows = score(buildings, P, p)
    result = dict(
        tool="score_lod22", schema=SCHEMA_VERSION,
        mode="all" if p["use_all"] else "hold",
        params=p,
        inputs=dict(city=[dict(path=f, bytes=os.path.getsize(f)) for f in _collect_city_files(args.city)],
                    pc=[dict(path=f, bytes=os.path.getsize(f)) for f in args.pc]),
        summary=summary,
    )
    print_summary(summary)
    print(f"# mode={result['mode']} points_loaded={len(P.x)} time={time.time() - t0:.1f}s", file=sys.stderr)
    if args.out_json:
        with open(args.out_json, "w", encoding="utf-8") as fh:
            json.dump(_json_safe(result), fh, indent=1)
    if args.out_csv:
        _write_csv(args.out_csv, brows)
    if args.out_face_csv:
        _write_csv(args.out_face_csv, frows)
    return result, brows, frows


def run_compare(a_path, b_path, out=None):
    out = out or sys.stdout
    with open(a_path, encoding="utf-8") as fh:
        a = json.load(fh)
    with open(b_path, encoding="utf-8") as fh:
        b = json.load(fh)
    if a.get("params") != b.get("params"):
        diff = sorted(k for k in set(a.get("params", {})) | set(b.get("params", {}))
                      if a.get("params", {}).get(k) != b.get("params", {}).get(k))
        print(f"# WARNING: params differ: {', '.join(diff)}", file=out)
    sa, sb = a["summary"], b["summary"]
    keys = list(sa) + [k for k in sb if k not in sa]
    w = max(len(k) for k in keys)
    print(f"{'metric':<{w}}  {'A':>14}  {'B':>14}  {'B-A':>14}  {'rel':>8}", file=out)
    rows = {}
    for k in keys:
        va, vb = sa.get(k), sb.get(k)
        num = lambda v: isinstance(v, (int, float)) and v is not None
        if num(va) and num(vb):
            d = vb - va
            rel = f"{100 * d / abs(va):+.1f}%" if va else ""
            tag = ""
            if d and k in LOWER_BETTER:
                tag = " better" if d < 0 else " worse"
            elif d and k in HIGHER_BETTER:
                tag = " better" if d > 0 else " worse"
            rows[k] = d
            print(f"{k:<{w}}  {_fmt(va):>14}  {_fmt(vb):>14}  {_fmt(d) if isinstance(d, float) else f'{d:+d}':>14}  {rel:>8}{tag}", file=out)
        else:
            print(f"{k:<{w}}  {_fmt(va):>14}  {_fmt(vb):>14}", file=out)
    return rows


def build_parser():
    ap = argparse.ArgumentParser(prog="score_lod22.py", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    D = DEFAULTS
    ap.add_argument("--city", nargs="+", required=True, help="roofer .city.jsonl file(s) or directories")
    ap.add_argument("--pc", nargs="+", required=True, help="LAS/LAZ file(s) used for reconstruction")
    ap.add_argument("--classes", type=int, nargs="+", default=D["classes"], help="point classes to use")
    ap.add_argument("--all", dest="use_all", action="store_true", help="use all points for decisions and scores (no hold-out)")
    ap.add_argument("--id-attribute", default=None, help="building attribute to use as id (default: CityObject id)")
    ap.add_argument("--out-json", default=None)
    ap.add_argument("--out-csv", default=None, help="per-building CSV")
    ap.add_argument("--out-face-csv", default=None, help="per-face CSV")
    g = ap.add_argument_group("thresholds")
    for k, h in [("inlier_dist", "|r| inlier distance (m)"), ("rms_max", "face is bad if rms_in exceeds this (m)"),
                 ("cell", "raster cell size (m)"), ("missed_min", "min height above plane for missed (m)"),
                 ("missed_max", "max height above plane for missed (m)"), ("missed_margin", "negative buffer for missed cells (m)"),
                 ("blob_min", "min blob area (m2)"), ("blob_max", "max blob area (m2)"),
                 ("adj_dist", "adjacency distance (m)"), ("adj_len", "min shared boundary length (m)"),
                 ("split_factor", "split: merged rms <= max(own)*factor + add"), ("split_add", "split additive term (m)"),
                 ("split_keep_inliers", "split guard: merged inliers >= this * own inliers"),
                 ("off", "OFF residual threshold (m)"), ("usable_margin", "negative buffer for usable cells (m)"),
                 ("level_deg", "level face max tilt (deg)"), ("sliver_area", "sliver if area_3d below (m2)"),
                 ("sliver_thinness", "sliver if area_3d/perimeter_3d^2 below"), ("min_nz", "faces with n_z below are 'steep'")]:
        g.add_argument("--" + k.replace("_", "-"), dest=k, type=float, default=D[k], help=h)
    g.add_argument("--split-min-pts", dest="split_min_pts", type=int, default=D["split_min_pts"])
    return ap


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "compare":
        cp = argparse.ArgumentParser(prog="score_lod22.py compare")
        cp.add_argument("a")
        cp.add_argument("b")
        a = cp.parse_args(argv[1:])
        run_compare(a.a, a.b)
        return 0
    if argv and argv[0] == "score":
        argv = argv[1:]
    run_score(build_parser().parse_args(argv))
    return 0


if __name__ == "__main__":
    sys.exit(main())
