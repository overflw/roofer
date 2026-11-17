"""Synthetic test tile: gable roof, flat green roof, tree crown over the gable roof.
Writes synth.laz (PF8 with NIR) and synth_dop.tif (4-band RGBI, 0.2 m)."""
import laspy
import numpy as np
import rasterio
from rasterio.transform import from_origin

rng = np.random.default_rng(0)
X0, Y0 = 500000.0, 5700000.0
P = []  # x,y,z,rn,nr,red,nir,cls,label
# gable roof 20x10 m, ridge along x, pitch 35 deg, 10 pts/m2
n = 2000
x = rng.uniform(0, 20, n); y = rng.uniform(0, 10, n)
z = 10 + np.tan(np.radians(35)) * (5 - np.abs(y - 5)) + rng.normal(0, 0.03, n)
P.append(np.column_stack([x, y, z, np.ones(n), np.ones(n), np.full(n, 30000), np.full(n, 33000), np.full(n, 6), np.zeros(n)]))
# flat green roof 15x10 m at x 30..45
n = 1500
x = rng.uniform(30, 45, n); y = rng.uniform(0, 10, n); z = 8 + rng.normal(0, 0.03, n)
P.append(np.column_stack([x, y, z, np.ones(n), np.ones(n), np.full(n, 12000), np.full(n, 36000), np.full(n, 6), np.ones(n)]))
# tree crown radius 3 m over gable at x=15,y=7, mostly first returns of 2-3 return pulses
n = 400
r = 3 * np.sqrt(rng.uniform(0, 1, n)); t = rng.uniform(0, 2 * np.pi, n)
x = 15 + r * np.cos(t); y = 7 + r * np.sin(t); z = 16 + rng.normal(0, 0.8, n) - 0.3 * r
nr = rng.choice([1, 2, 3], n, p=[0.1, 0.5, 0.4]); rn = np.where(nr > 1, rng.integers(1, np.maximum(nr, 2)), 1)
P.append(np.column_stack([x, y, z, rn, nr, np.full(n, 10000), np.full(n, 30000), np.full(n, 6), np.full(n, 2)]))
A = np.vstack(P)

h = laspy.LasHeader(point_format=8, version="1.4")
h.offsets = [X0, Y0, 0]; h.scales = [0.001] * 3
las = laspy.LasData(h)
las.x = A[:, 0] + X0; las.y = A[:, 1] + Y0; las.z = A[:, 2]
las.return_number = A[:, 3].astype(int); las.number_of_returns = A[:, 4].astype(int)
las.red = A[:, 5].astype(np.uint16); las.green = A[:, 5].astype(np.uint16); las.blue = A[:, 5].astype(np.uint16)
las.nir = A[:, 6].astype(np.uint16)
las.classification = A[:, 7].astype(np.uint8)
las.write("synth.laz")
np.save("synth_labels.npy", A[:, 8])

# DOP: 0.2 m RGBI, tree is drawn shifted 1 m (lean / misregistration test)
res = 0.2; W, H = int(50 / res), int(12 / res)
img = np.zeros((4, H, W), np.uint8)
img[0] = 120; img[3] = 130  # roofs: NDVI ~0.04
xs = (np.arange(W) + 0.5) * res; ys = 12 - (np.arange(H) + 0.5) * res
XX, YY = np.meshgrid(xs, ys)
green = (XX > 30) & (XX < 45) & (YY < 10)
img[0][green] = 50; img[3][green] = 150
tree = (XX - 16) ** 2 + (YY - 7) ** 2 < 9
img[0][tree] = 40; img[3][tree] = 140
with rasterio.open("synth_dop.tif", "w", driver="GTiff", width=W, height=H, count=4, dtype="uint8",
                   crs="EPSG:25832", transform=from_origin(X0, Y0 + 12, res, res)) as ds:
    ds.write(img)
