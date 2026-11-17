"""Fetch a bbox from a USGS EPT dataset into a LAZ in a target CRS.
usage: ept_fetch.py DATASET lon_min lat_min lon_max lat_max out.laz EPSG"""
import sys, json, io, urllib.request, concurrent.futures as cf
import numpy as np, laspy, pyproj
ds, x0, y0, x1, y1, out, epsg = sys.argv[1], *map(float, sys.argv[2:6]), sys.argv[6], int(sys.argv[7])
base = f"https://s3-us-west-2.amazonaws.com/usgs-lidar-public/{ds}/"
get = lambda p: urllib.request.urlopen(base + p, timeout=60).read()
ept = json.loads(get("ept.json"))
t = pyproj.Transformer.from_crs(4326, 3857, always_xy=True)
(mx0, my0), (mx1, my1) = t.transform(x0, y0), t.transform(x1, y1)
B = ept["bounds"]
def nb(key):
    d, x, y, z = map(int, key.split("-")); n = 2 ** d
    w = (B[3] - B[0]) / n; h = (B[4] - B[1]) / n
    return B[0] + x * w, B[1] + y * h, B[0] + (x + 1) * w, B[1] + (y + 1) * h
def hit(key):
    a, b, c, e = nb(key); return not (c < mx0 or a > mx1 or e < my0 or b > my1)
keys = []
def walk(hkey):
    h = json.loads(get(f"ept-hierarchy/{hkey}.json"))
    for k, v in h.items():
        if not hit(k): continue
        if v == -1:
            if k != hkey: walk(k)
        elif v > 0: keys.append(k)
walk("0-0-0-0")
print(len(keys), "nodes")
def fetch(k):
    las = laspy.read(io.BytesIO(get(f"ept-data/{k}.laz")))
    m = (las.x >= mx0) & (las.x <= mx1) & (las.y >= my0) & (las.y <= my1)
    return las.points[m]
with cf.ThreadPoolExecutor(16) as ex: parts = list(ex.map(fetch, keys))
hdr = laspy.LasHeader(point_format=parts[0].point_format.id, version="1.4")
X = np.concatenate([np.asarray(p.x) for p in parts]); Y = np.concatenate([np.asarray(p.y) for p in parts]); Z = np.concatenate([np.asarray(p.z) for p in parts])
tt = pyproj.Transformer.from_crs(3857, epsg, always_xy=True)
X2, Y2 = tt.transform(X, Y)
hdr.offsets = [np.floor(X2.min()), np.floor(Y2.min()), 0]; hdr.scales = [0.001] * 3
hdr.add_crs(pyproj.CRS.from_epsg(epsg))
o = laspy.LasData(hdr)
o.x, o.y, o.z = X2, Y2, Z
for f in ["classification", "return_number", "number_of_returns", "intensity", "red", "green", "blue", "gps_time"]:
    try: setattr(o, f, np.concatenate([np.asarray(p[f]) for p in parts]))
    except Exception as e: print("skip", f, e)
o.write(out)
area = (X2.max() - X2.min()) * (Y2.max() - Y2.min())
print(len(X2), "pts", len(X2) / area, "pts/m2", np.bincount(o.classification))
