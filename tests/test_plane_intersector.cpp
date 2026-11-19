#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <roofer/reconstruction/PlaneIntersector.hpp>

using namespace roofer;
using namespace roofer::reconstruction;

// a gable along y: plane 1 east of x = 0, plane 2 west of it, each with points
// along the ridge over [y0, y1]
static IndexedPlanesWithPoints gable(double y0_1, double y1_1, double y0_2,
                                     double y1_2) {
  IndexedPlanesWithPoints planes;
  auto side = [](double sx, double y0, double y1) {
    std::vector<Point> pts;
    for (double y = y0; y <= y1 + 1e-9; y += 0.25)
      for (double x = 0.1; x <= 0.5; x += 0.1)
        pts.emplace_back(sx * x, y, 10 - 0.5 * x);
    return pts;
  };
  planes[1] = {Plane(0.5, 0, 1, -10), side(1, y0_1, y1_1)};
  planes[2] = {Plane(-0.5, 0, 1, -10), side(-1, y0_2, y1_2)};
  return planes;
}

static const std::map<size_t, std::map<size_t, size_t>> adjacent = {
    {1, {{2, 10}}}};

TEST_CASE("an intersection line spans where both planes have points",
          "[plane_intersector]") {
  auto intersector = createPlaneIntersector();
  intersector->compute(gable(0, 5, 3, 12), adjacent);
  REQUIRE(intersector->segments.size() == 1);
  auto& s = intersector->segments[0];
  CHECK(std::min(s[0][1], s[1][1]) == Catch::Approx(3));
  CHECK(std::max(s[0][1], s[1][1]) == Catch::Approx(5));
}

TEST_CASE("planes with disjoint extents along their intersection give no line",
          "[plane_intersector]") {
  auto intersector = createPlaneIntersector();
  intersector->compute(gable(0, 5, 7, 12), adjacent);
  CHECK(intersector->segments.empty());
}
