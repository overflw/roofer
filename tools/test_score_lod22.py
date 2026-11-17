"""Self-test for score_lod22.py on a synthetic gable roof with a box dormer.

Run: pytest tools/test_score_lod22.py
"""

import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score_lod22 as S  # noqa: E402

X0, Y0 = 85000.0, 446000.0  # realistic projected offsets
EAVE, SLOPE = 6.0, 0.7  # z = 6 + 0.7 * distance from eave; ridge at y=5 -> 9.5
RIDGE = EAVE + 5 * SLOPE
DORMER = (4.0, 1.5, 8.0, 3.8)  # x0, y0, x1, y1 (local), on the south face
DORMER_Z = 9.4
DORMER_AREA = (DORMER[2] - DORMER[0]) * (DORMER[3] - DORMER[1])  # 9.2 m2


def z_south(y):
    return EAVE + SLOPE * y


def z_north(y):
    return RIDGE - SLOPE * (y - 5.0)


def make_points(path, density=30.0, noise=0.03, seed=1):
    import laspy

    rng = np.random.default_rng(seed)
    n = int(12 * 10 * density)
    x = rng.uniform(0, 12, n)
    y = rng.uniform(0, 10, n)
    z = np.where(y < 5, z_south(y), z_north(y))
    d = (x > DORMER[0]) & (x < DORMER[2]) & (y > DORMER[1]) & (y < DORMER[3])
    z[d] = DORMER_Z
    z += rng.normal(0, noise, n)
    cls = np.full(n, 6, np.uint8)
    # ground points around (class 2) must be ignored
    m = 2000
    gx = rng.uniform(-5, 17, m)
    gy = rng.uniform(-5, 15, m)
    gz = rng.normal(0, 0.02, m)
    # a few ground-class points *inside* the footprint at roof height would be
    # counted if the class filter were broken: make them obvious outliers
    gz[:50] = 20.0
    gx[:50], gy[:50] = rng.uniform(1, 11, 50), rng.uniform(1, 9, 50)
    hdr = laspy.LasHeader(point_format=3, version="1.2")
    hdr.scales = np.array([0.001, 0.001, 0.001])
    hdr.offsets = np.array([X0, Y0, 0.0])
    las = laspy.LasData(hdr)
    order = rng.permutation(n + m)
    las.x = (np.r_[x, gx] + X0)[order]
    las.y = (np.r_[y, gy] + Y0)[order]
    las.z = np.r_[z, gz][order]
    las.classification = np.r_[cls, np.full(m, 2, np.uint8)][order]
    las.write(path)


class CJ:
    """Tiny CityJSONFeature builder with roofer's layout (Building + BuildingPart)."""

    scale = (0.001, 0.001, 0.001)
    translate = (X0, Y0, 0.0)

    def __init__(self):
        self.verts = []
        self.vmap = {}

    def vi(self, p):
        q = tuple(int(round((c - 0) / s)) for c, s in zip(p, self.scale))
        if q not in self.vmap:
            self.vmap[q] = len(self.verts)
            self.verts.append(list(q))
        return self.vmap[q]

    def ring(self, pts):
        return [self.vi(p) for p in pts]


def south_face(cj, hole):
    ext = [(0, 0, z_south(0)), (12, 0, z_south(0)), (12, 5, RIDGE), (0, 5, RIDGE)]
    rings = [cj.ring(ext)]
    if hole:
        x0, y0, x1, y1 = DORMER
        rings.append(cj.ring([(x0, y0, z_south(y0)), (x0, y1, z_south(y1)), (x1, y1, z_south(y1)), (x1, y0, z_south(y0))]))
    return rings


def north_faces(cj, split):
    if not split:
        return [[cj.ring([(0, 5, RIDGE), (12, 5, RIDGE), (12, 10, z_north(10)), (0, 10, z_north(10))])]]
    return [[cj.ring([(0, 5, RIDGE), (6, 5, RIDGE), (6, 10, z_north(10)), (0, 10, z_north(10))])],
            [cj.ring([(6, 5, RIDGE), (12, 5, RIDGE), (12, 10, z_north(10)), (6, 10, z_north(10))])]]


def dormer_face(cj):
    x0, y0, x1, y1 = DORMER
    return [cj.ring([(x0, y0, DORMER_Z), (x1, y0, DORMER_Z), (x1, y1, DORMER_Z), (x0, y1, DORMER_Z)])]


def make_feature(variant, geom_type="Solid"):
    cj = CJ()
    surfaces, types = [], []
    surfaces.append([cj.ring([(0, 0, 0), (0, 10, 0), (12, 10, 0), (12, 0, 0)])])
    types.append(0)  # ground
    surfaces.append(south_face(cj, hole=(variant != "nodormer")))
    types.append(3)
    for f in north_faces(cj, split=(variant == "split")):
        surfaces.append(f)
        types.append(3)
    if variant != "nodormer":
        surfaces.append(dormer_face(cj))
        types.append(3)
        # a dormer wall (vertical; must be ignored)
        x0, y0, x1, _ = DORMER
        surfaces.append([cj.ring([(x0, y0, z_south(y0)), (x1, y0, z_south(y0)), (x1, y0, DORMER_Z), (x0, y0, DORMER_Z)])])
        types.append(1)
    sem_surfaces = [{"type": "GroundSurface"}, {"type": "WallSurface", "on_footprint_edge": True},
                    {"type": "WallSurface", "on_footprint_edge": False}]
    values = []
    for t in types:
        if t == 3:
            sem_surfaces.append({"type": "RoofSurface"})
            values.append(len(sem_surfaces) - 1)
        else:
            values.append(t)
    if geom_type == "Solid":
        geom = {"type": "Solid", "lod": "2.2", "boundaries": [surfaces], "semantics": {"surfaces": sem_surfaces, "values": [values]}}
    else:
        geom = {"type": "MultiSurface", "lod": "2.2", "boundaries": surfaces, "semantics": {"surfaces": sem_surfaces, "values": values}}
    fp = {"type": "MultiSurface", "lod": "0", "boundaries": [[cj.ring([(0, 0, 0), (12, 0, 0), (12, 10, 0), (0, 10, 0)])]]}
    feat = {
        "type": "CityJSONFeature", "id": "B1",
        "CityObjects": {
            "B1": {"type": "Building", "attributes": {"identificatie": "NL.B1", "rf_extrusion_mode": "standard",
                                                      "rf_roof_planes": 3}, "geometry": [fp], "children": ["B1-0"]},
            "B1-0": {"type": "BuildingPart", "parents": ["B1"], "geometry": [geom]},
        },
        "vertices": cj.verts,
    }
    # a fallback building with only a footprint (e.g. skipped / failed)
    cj2 = CJ()
    fb = {
        "type": "CityJSONFeature", "id": "B2",
        "CityObjects": {"B2": {"type": "Building", "attributes": {"identificatie": "NL.B2", "rf_extrusion_mode": "skip"},
                               "geometry": [{"type": "MultiSurface", "lod": "0", "boundaries": [[cj2.ring(
                                   [(30, 0, 0), (35, 0, 0), (35, 5, 0), (30, 5, 0)])]]}], "children": []}},
        "vertices": cj2.verts,
    }
    return [feat, fb]


def write_seq(path, feats):
    header = {"type": "CityJSON", "version": "2.0", "CityObjects": {}, "vertices": [],
              "transform": {"scale": list(CJ.scale), "translate": list(CJ.translate)}}
    with open(path, "w") as fh:
        fh.write(json.dumps(header) + "\n")
        for f in feats:
            fh.write(json.dumps(f) + "\n")


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    d = tmp_path_factory.mktemp("score")
    pc = str(d / "pts.laz")
    make_points(pc)
    out = {"pc": pc, "dir": d}
    for name, gt in [("good", "Solid"), ("split", "MultiSurface"), ("nodormer", "Solid")]:
        p = str(d / f"{name}.city.jsonl")
        write_seq(p, make_feature(name, gt))
        out[name] = p
    return out


def run(data, name, *extra):
    d = data["dir"]
    js = str(d / f"{name}{'_'.join(extra).replace('-', '')}.json")
    args = ["--city", data[name], "--pc", data["pc"], "--out-json", js,
            "--out-csv", str(d / f"{name}_b.csv"), "--out-face-csv", str(d / f"{name}_f.csv"), *extra]
    res, brows, frows = S.run_score(S.build_parser().parse_args(args))
    return res["summary"], brows, frows, js


def test_good(data):
    s, brows, frows, _ = run(data, "good")
    assert s["buildings"] == 2 and s["buildings_fallback"] == 1
    assert s["extrusion_mode_skip"] == 1
    assert s["faces"] == 3  # south (with hole), north, dormer roof; walls/ground ignored
    assert s["faces_level"] == 1  # flat dormer roof
    assert s["faces_nodata"] == 0
    assert s["splits"] == 0 and s["crease_m2"] == 0
    assert s["missed_m2"] == 0
    assert s["bad_m2"] == 0
    assert 0.02 < s["rms_in_face_p50"] < 0.045
    assert s["outlier_share"] < 0.01  # class-2 outliers at 20 m are filtered out
    assert abs(s["roof_area_2d_m2"] - 120.0) < 1e-6
    assert s["usable_frac_eligible"] > 0.9 and s["off_m2"] < 1.0
    b1 = [b for b in brows if b["building_id"] == "B1"][0]
    assert b1["n_parts"] == 1 and b1["n_faces"] == 3


def test_split(data):
    s, _, frows, _ = run(data, "split")
    assert s["faces"] == 4
    assert s["splits"] == 1
    assert abs(s["crease_m2"] - 30.0) < 1e-6
    assert sum(f["crease"] for f in frows) == 1
    assert s["missed_m2"] == 0


def test_nodormer(data):
    s, _, frows, _ = run(data, "nodormer")
    assert s["faces"] == 2
    assert s["splits"] == 0
    assert abs(s["missed_m2"] - DORMER_AREA) < 2.0, s["missed_m2"]
    assert s["missed_blobs"] == 1
    assert abs(s["bad_m2"] - s["missed_m2"]) < 1e-9  # no face is rms-bad
    assert s["off_m2"] > 0.7 * DORMER_AREA


def test_all_mode_and_compare(data, capsys):
    s_all, *_ = run(data, "nodormer", "--all")
    assert s_all["points_scored"] > 3000
    _, _, _, ja = run(data, "good")
    _, _, _, jb = run(data, "nodormer")
    rows = S.run_compare(ja, jb)
    out = capsys.readouterr().out
    assert "missed_m2" in out and rows["faces"] == -1
    assert S.main(["compare", ja, jb]) == 0


def test_id_attribute(data):
    args = ["--city", data["good"], "--pc", data["pc"], "--id-attribute", "identificatie"]
    _, brows, _ = S.run_score(S.build_parser().parse_args(args))
    assert sorted(b["building_id"] for b in brows) == ["NL.B1", "NL.B2"]
