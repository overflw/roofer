#include <catch2/catch_test_macros.hpp>

#include <roofer/reconstruction/ArrangementOptimiser.hpp>

namespace {
  using roofer::Point_2;
  using roofer::Segment_2;

  // A 10 x 10 footprint split into four 5 x 5 quadrants. The heightfield lies
  // on plane 1 (z = 10) except for the top right quadrant at z = 10.6, which
  // fits plane 2 (z = 11) slightly better.
  int top_right_segid(float complexity_factor) {
    roofer::Arrangement_2 arr;
    for (auto [a, b] : {std::pair{Point_2(0, 0), Point_2(10, 0)},
                        {Point_2(10, 0), Point_2(10, 10)},
                        {Point_2(10, 10), Point_2(0, 10)},
                        {Point_2(0, 10), Point_2(0, 0)},
                        {Point_2(5, 0), Point_2(5, 10)},
                        {Point_2(0, 5), Point_2(10, 5)}}) {
      CGAL::insert(arr, Segment_2(a, b));
    }
    for (auto face : arr.face_handles()) {
      if (!face->is_unbounded()) face->data().in_footprint = true;
    }

    roofer::RasterTools::Raster heightfield(1, 0, 10, 0, 10);
    heightfield.prefill_arrays(roofer::RasterTools::MAX);
    for (int x = 0; x < 10; ++x) {
      for (int y = 0; y < 10; ++y) {
        heightfield.add_point(x + 0.5, y + 0.5, x >= 5 && y >= 5 ? 10.6 : 10,
                              roofer::RasterTools::MAX);
      }
    }

    roofer::IndexedPlanesWithPoints roof_planes{
        {1, {roofer::Plane(0, 0, 1, -10), {}}},
        {2, {roofer::Plane(0, 0, 1, -11), {}}}};
    roofer::reconstruction::createArrangementOptimiser()->compute(
        arr, heightfield, roof_planes, {},
        {.complexity_factor = complexity_factor, .use_ground = false});

    for (auto face : arr.face_handles()) {
      if (face->is_unbounded()) continue;
      bool top_right = true;
      auto he = face->outer_ccb();
      do {
        top_right &= he->source()->point().x() >= 5 &&
                     he->source()->point().y() >= 5;
      } while (++he != face->outer_ccb());
      if (top_right) return face->data().segid;
    }
    return 0;
  }
}  // namespace

// The top right quadrant saves 0.2 data cost per cell on plane 2 and pays its
// 10 m of inner edges in smoothness, each metre counted twice (once per
// direction). That puts its switch to plane 2 between complexity factors 0.75
// and 0.85; counting each edge once would move it below 0.75.
TEST_CASE("graph cut weighs a shared edge twice its length") {
  CHECK(top_right_segid(0.75F) == 1);
  CHECK(top_right_segid(0.85F) == 2);
}
