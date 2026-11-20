#include <catch2/catch_test_macros.hpp>

#include <roofer/reconstruction/ArrangementOptimiser.hpp>

using namespace roofer;

namespace {
  // two 10 x 10 faces side by side, x in [0, 10) and [10, 20)
  Arrangement_2 two_faces() {
    Arrangement_2 arr;
    std::vector<Segment_2> segments = {{Point_2(0, 0), Point_2(20, 0)},
                                       {Point_2(20, 0), Point_2(20, 10)},
                                       {Point_2(20, 10), Point_2(0, 10)},
                                       {Point_2(0, 10), Point_2(0, 0)},
                                       {Point_2(10, 0), Point_2(10, 10)}};
    CGAL::insert(arr, segments.begin(), segments.end());
    for (auto face : arr.face_handles()) {
      if (!face->is_unbounded()) face->data().in_footprint = true;
    }
    return arr;
  }

  // ground level on the left face, z on the right one
  RasterTools::Raster heightfield(double z) {
    RasterTools::Raster raster(1, 0, 20, 0, 10);
    raster.prefill_arrays(RasterTools::MAX);
    for (double x = 0.5; x < 20; ++x)
      for (double y = 0.5; y < 10; ++y)
        raster.add_point(x, y, x < 10 ? 0 : z, RasterTools::MAX);
    return raster;
  }
}  // namespace

TEST_CASE("a capped data cost keeps a raised face off the ground",
          "[arrangementoptimiser]") {
  auto arr = two_faces();
  // the only roof plane is 2 m above the raised face, ground is at 0
  IndexedPlanesWithPoints roof = {{1, {Plane(0, 0, 1, -6), {}}}};
  IndexedPlanesWithPoints ground = {{1, {Plane(0, 0, 1, 0), {}}}};
  reconstruction::ArrangementOptimiserConfig cfg;
  cfg.max_data_distance = 1;
  reconstruction::createArrangementOptimiser()->compute(arr, heightfield(4),
                                                        roof, ground, cfg);

  size_t raised = 0, ground_faces = 0;
  for (auto face : arr.face_handles()) {
    if (face->is_unbounded()) continue;
    Point_2 c = CGAL::centroid(face->outer_ccb()->source()->point(),
                               face->outer_ccb()->target()->point(),
                               face->outer_ccb()->next()->target()->point());
    if (face->data().is_ground) {
      ++ground_faces;
      CHECK(c.x() < 10);
    } else {
      ++raised;
      CHECK(c.x() > 10);
    }
  }
  CHECK(raised == 1);
  CHECK(ground_faces == 1);
}
