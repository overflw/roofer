"""Derive building outlines from class-6 points (stand-in for ALKIS footprints).
usage: make_footprints.py in.laz out.gpkg EPSG"""
import sys, numpy as np, laspy, scipy.ndimage as nd, rasterio.features as rf, shapely, geopandas as gpd
from affine import Affine
las = laspy.read(sys.argv[1]); m = las.classification == 6
x, y = np.asarray(las.x)[m], np.asarray(las.y)[m]
c = 0.5; x0, y1 = x.min(), y.max()
ix = ((x - x0) / c).astype(int); iy = ((y1 - y) / c).astype(int)
g = np.zeros((iy.max() + 1, ix.max() + 1), np.uint8); g[iy, ix] = 1
g = nd.binary_closing(g, iterations=2); g = nd.binary_fill_holes(g); g = nd.binary_opening(g, iterations=1)
lab, n = nd.label(g)
tr = Affine(c, 0, x0, 0, -c, y1)
polys = []
for geom, v in rf.shapes(lab.astype(np.int32), mask=lab > 0, transform=tr):
    p = shapely.geometry.shape(geom).simplify(0.4).buffer(0)
    if p.area >= 40: polys.append(p)
gdf = gpd.GeoDataFrame({"identificatie": [f"b{i}" for i in range(len(polys))]}, geometry=polys, crs=int(sys.argv[3]))
gdf.to_file(sys.argv[2], driver="GPKG")
print(len(polys), "footprints, area", round(gdf.area.sum()))
