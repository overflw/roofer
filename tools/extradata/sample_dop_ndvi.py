#!/usr/bin/env python3
"""Sketch: sample NDVI from an RGBI DOP (COG / VRT / /vsicurl/ URL) at point XY.

Reads only the window covering the points (COG range requests), computes NDVI on a
3x3-median-smoothed window to absorb ~1 px true-ortho misregistration, then samples
nearest pixel. Points must be in the raster CRS (German DOPs: EPSG:25832/25833).

Examples of sources (not reachable from the sandbox, proxy 403):
  Niedersachsen: STAC https://dop.stac.lgln.niedersachsen.de/ -> COG asset href
  NRW:           https://www.opengeodata.nrw.de/produkte/geobasis/lusat/akt/dop/dop_jp2_f10/ (JP2, build a VRT)
Use '/vsicurl/<href>' with GDAL_DISABLE_READDIR_ON_OPEN=EMPTY_DIR for remote COGs.
"""
import numpy as np
import rasterio
from rasterio.windows import from_bounds


def _median3(a):
    p = np.pad(a, 1, mode="edge")
    st = np.stack([p[i:i + a.shape[0], j:j + a.shape[1]] for i in range(3) for j in range(3)])
    return np.median(st, axis=0)


def sample_ndvi(src_path, x, y, red_band=1, nir_band=4, smooth=True, pad=2.0):
    env = rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                       CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff,.jp2,.vrt")
    with env, rasterio.open(src_path) as ds:
        win = from_bounds(x.min() - pad, y.min() - pad, x.max() + pad, y.max() + pad,
                          ds.transform).round_offsets().round_lengths()
        r = ds.read(red_band, window=win, boundless=True, fill_value=0).astype(np.float32)
        n = ds.read(nir_band, window=win, boundless=True, fill_value=0).astype(np.float32)
        nd = (n - r) / np.maximum(n + r, 1.0)
        nd[(n + r) == 0] = np.nan
        if smooth:
            nd = _median3(nd)
        wt = ds.window_transform(win)
        rows, cols = rasterio.transform.rowcol(wt, x, y)
        rows = np.clip(np.asarray(rows), 0, nd.shape[0] - 1)
        cols = np.clip(np.asarray(cols), 0, nd.shape[1] - 1)
        return nd[rows, cols]


def shadow_mask(src_path, x, y, bands=(1, 2, 3), thresh=40):
    """Low-brightness pixels (shadow / dark roofs) -> NDVI unreliable; caller should
    treat NDVI as unknown there."""
    with rasterio.open(src_path) as ds:
        vals = np.stack([v for v in ds.sample(zip(x, y), indexes=list(bands))])
    return vals.mean(axis=1) < thresh


if __name__ == "__main__":
    import sys
    import laspy
    las = laspy.read(sys.argv[2])
    m = np.asarray(las.classification) == 6
    nd = sample_ndvi(sys.argv[1], np.asarray(las.x)[m], np.asarray(las.y)[m])
    print("NDVI p10/p50/p90 on class 6:", np.nanpercentile(nd, [10, 50, 90]))
