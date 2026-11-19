#include <catch2/catch_test_macros.hpp>
#include <roofer/reconstruction/LineDetectorBase.hpp>

#include <set>

using namespace roofer::linedect;

static std::set<size_t> segments(const LineDetector& ld) {
  return {ld.point_segment_idx.begin(), ld.point_segment_idx.end()};
}

TEST_CASE("the line detector uses the given number of neighbours",
          "[line_detector]") {
  // two collinear runs of 5 points, 1 m apart: 5 neighbours (the point
  // itself included) don't reach across the gap, 10 do
  std::vector<Point> pts;
  for (double x : {0.0, 0.1, 0.2, 0.3, 0.4, 1.4, 1.5, 1.6, 1.7, 1.8})
    pts.emplace_back(x, 0, 0);
  for (auto [n, expected] : {std::pair{5, 2}, std::pair{10, 1}}) {
    auto copy = pts;
    LineDetector ld(copy, n);
    ld.detect();
    CHECK(segments(ld).size() == size_t(expected));
  }
}
