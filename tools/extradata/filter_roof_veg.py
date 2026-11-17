#!/usr/bin/env python3
"""Prototype: remove vegetation-like points from the building class before roofer.

Per class-`--bld-class` point we compute three cues:
  multi   : return_number < number_of_returns  (non-last return of a multi-return pulse)
  ndvi    : (NIR - Red) / (NIR + Red) from per-point NIR (LAS PF 8/10 or extra dim),
            or sampled from a DOP RGBI raster (--dop, see sample_dop_ndvi.py)
  planar  : local PCA on k nearest building points; surface variation
            l3/(l1+l2+l3) < --max-sv  AND  plane RMS < --max-rms  (green-roof / PV safety)

A point is reclassified (default to class 5, high vegetation; or dropped with --drop)
only if it is NOT on a locally planar patch AND (multi OR ndvi > --ndvi). Planar
high-NDVI points (extensive green roofs) and planar multi-return points (roof edges,
glass) are kept. A second, stricter rule removes planar points too when BOTH cues fire
and NDVI is very high (--ndvi-hard), aimed at dense flat canopy; disable with --ndvi-hard 2.

Writes stats to stdout. Footprint-agnostic: run per tile before roofer.
"""
import argparse
import sys

import laspy
import numpy as np
from scipy.spatial import cKDTree


def local_pca(xyz, k):
    """Return surface variation and plane RMS per point (vectorised, chunked)."""
    tree = cKDTree(xyz[:, :2])  # 2D neighbourhood: roof-like footprint in plan
    n = len(xyz)
    sv = np.empty(n, np.float32)
    rms = np.empty(n, np.float32)
    k = min(k, n)
    for s in range(0, n, 200_000):
        e = min(n, s + 200_000)
        _, idx = tree.query(xyz[s:e, :2], k=k)
        nb = xyz[idx]  # (m,k,3)
        c = nb - nb.mean(axis=1, keepdims=True)
        cov = np.einsum("mki,mkj->mij", c, c) / k
        w = np.linalg.eigvalsh(cov)  # ascending
        tot = w.sum(axis=1) + 1e-12
        sv[s:e] = w[:, 0] / tot
        rms[s:e] = np.sqrt(np.maximum(w[:, 0], 0))
    return sv, rms


def get_dim(las, names):
    have = {d.lower(): d for d in las.point_format.dimension_names}
    for n in names:
        if n.lower() in have:
            return np.asarray(las[have[n.lower()]], dtype=np.float64)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--bld-class", type=int, default=6)
    ap.add_argument("--to-class", type=int, default=5)
    ap.add_argument("--drop", action="store_true")
    ap.add_argument("--k", type=int, default=12)
    ap.add_argument("--max-sv", type=float, default=0.01)
    ap.add_argument("--max-rms", type=float, default=0.08, help="metres")
    ap.add_argument("--ndvi", type=float, default=0.25)
    ap.add_argument("--ndvi-hard", type=float, default=0.55)
    ap.add_argument("--no-multi", action="store_true", help="ignore return info")
    ap.add_argument("--dop", help="RGBI raster (COG/VRT/URL) to sample NDVI from")
    ap.add_argument("--dop-red-band", type=int, default=1)
    ap.add_argument("--dop-nir-band", type=int, default=4)
    a = ap.parse_args()

    las = laspy.read(a.inp)
    cls = np.asarray(las.classification)
    m = cls == a.bld_class
    xyz = np.column_stack([las.x, las.y, las.z])[m]
    nb = int(m.sum())
    if nb == 0:
        las.write(a.out)
        print("no building points")
        return

    # cue 1: multi-return non-last
    if a.no_multi:
        multi = np.zeros(nb, bool)
    else:
        rn = np.asarray(las.return_number)[m]
        nr = np.asarray(las.number_of_returns)[m]
        multi = (nr > 1) & (rn < nr)

    # cue 2: NDVI (per-point NIR first, else DOP raster)
    ndvi = None
    nir = get_dim(las, ["nir", "infrared", "near_infrared"])
    if nir is not None and "red" in las.point_format.dimension_names:
        red = np.asarray(las.red, np.float64)
        ndvi = ((nir - red) / np.maximum(nir + red, 1.0))[m]
        src = "per-point NIR"
    elif a.dop:
        from sample_dop_ndvi import sample_ndvi

        ndvi = sample_ndvi(a.dop, xyz[:, 0], xyz[:, 1], a.dop_red_band, a.dop_nir_band)
        src = "DOP raster"
    if ndvi is None:
        ndvi = np.full(nb, np.nan)
        src = "none"
    veg = np.nan_to_num(ndvi, nan=-1) > a.ndvi
    hard = np.nan_to_num(ndvi, nan=-1) > a.ndvi_hard

    # cue 3: local planarity (safety)
    sv, rms = local_pca(xyz, a.k)
    planar = (sv < a.max_sv) & (rms < a.max_rms)

    bad = (~planar & (multi | veg)) | (planar & multi & hard)

    print(f"building pts {nb}; ndvi source: {src}")
    print(f"  multi-nonlast {multi.mean():.3f}  ndvi>{a.ndvi} {veg.mean():.3f}  planar {planar.mean():.3f}")
    print(f"  green-roof-safe kept (planar & ndvi>{a.ndvi} & !bad): {(planar & veg & ~bad).sum()}")
    print(f"  reclassified/dropped: {bad.sum()} ({bad.mean():.3%})")

    idx = np.flatnonzero(m)[bad]
    if a.drop:
        keep = np.ones(len(cls), bool)
        keep[idx] = False
        las.points = las.points[keep]
    else:
        cls = cls.copy()
        cls[idx] = a.to_class
        las.classification = cls
    las.write(a.out)


if __name__ == "__main__":
    sys.exit(main())
