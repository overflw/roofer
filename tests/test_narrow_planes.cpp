#include <roofer/reconstruction/NarrowPlanes.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <random>

using namespace roofer::reconstruction;
namespace {
  constexpr double r = 0.5, epsilon = 0.2;
  // A flat roof 20 x 10 m at 10 m, 16 points/m², raised by `extra`; a
  // raised point is unsegmented (label 0), the others plane 1.
  template <typename F>
  void roof(F extra, std::vector<std::array<double, 3>>& pts,
            std::vector<int>& labels) {
    std::mt19937 rng(0);
    std::normal_distribution<double> noise(0, 0.03);
    for (double x = 0.125; x < 20; x += 0.25)
      for (double y = 0.125; y < 10; y += 0.25) {
        double h = extra(x, y);
        pts.push_back({x, y, 10 + h + noise(rng)});
        labels.push_back(h > epsilon ? 0 : 1);
      }
  }
  std::vector<std::array<double, 3>> unsegmented(
      const std::vector<std::array<double, 3>>& pts,
      const std::vector<int>& labels) {
    std::vector<std::array<double, 3>> out;
    for (size_t i = 0; i < pts.size(); ++i)
      if (labels[i] == 0) out.push_back(pts[i]);
    return out;
  }
}  // namespace

TEST_CASE("a wall on the roof stands, a facade under the eaves does not") {
  std::vector<std::array<double, 3>> pts;
  std::vector<int> labels;
  roof([](double x, double) { return std::abs(x - 8) < 0.2 ? 0.6 : 0.0; }, pts,
       labels);
  size_t n = pts.size();
  for (double y = 0.125; y < 10; y += 0.25) {  // a facade 2 m under the roof
    pts.push_back({20.1, y, 8});
    labels.push_back(0);
  }
  auto stands = standing(pts, labels, 2 * r, epsilon);
  size_t wall = 0, on_wall = 0, facade = 0;
  for (size_t i = 0; i < pts.size(); ++i) {
    if (i < n && labels[i] == 0) ++wall, on_wall += stands[i];
    if (i >= n) facade += stands[i];
  }
  CHECK(on_wall == wall);
  CHECK(facade == 0);
}

TEST_CASE("walls of one height split into runs, a broad block does not") {
  std::vector<std::array<double, 3>> pts;
  std::vector<int> labels;
  roof(
      [](double x, double y) {
        if (std::abs(x - 5) < 0.2 || std::abs(x - 10) < 0.2) return 0.6;
        return x > 14 && x < 18 && y > 2 && y < 6 ? 0.6 : 0.0;
      },
      pts, labels);
  auto xyz = unsegmented(pts, labels);
  size_t narrow = 0;
  for (auto& run : narrow_runs(xyz, r, 4 * r, epsilon)) {
    std::vector<Point2> xy;
    for (auto i : run) xy.push_back({xyz[i][0], xyz[i][1]});
    auto o = narrow_outline(xy, r);
    if (!o) continue;
    ++narrow;
    double cx = ((*o)[0][0] + (*o)[2][0]) / 2;
    CHECK((std::abs(cx - 5) < 0.3 || std::abs(cx - 10) < 0.3));
    CHECK(std::hypot((*o)[1][0] - (*o)[0][0], (*o)[1][1] - (*o)[0][1]) > 9);
  }
  CHECK(narrow == 2);
}

TEST_CASE("structures of two heights on one line make two runs") {
  std::vector<std::array<double, 3>> xyz;
  for (double x = 0; x < 10; x += 0.25) {
    xyz.push_back({x, 5, 10.6});
    xyz.push_back({x + 0.1, 5.1, 11.4});
  }
  auto runs = narrow_runs(xyz, r, 4 * r, epsilon);
  REQUIRE(runs.size() == 2);
  for (auto& run : runs) {
    for (auto i : run) CHECK(xyz[i][2] == xyz[run[0]][2]);
    CHECK(run.size() == 40);
  }
}

TEST_CASE("a structure runs on over its points and to the border") {
  Rectangle r0{Point2{4, 4.9}, Point2{8, 4.9}, Point2{8, 5.1}, Point2{4, 5.1}};
  std::vector<std::array<double, 3>> line;
  for (double x = 2; x < 12; x += 0.5) line.push_back({x, 5, 10.6});
  auto r1 = narrow_extend(
      r0, line, [](double, double) { return 10.6; }, epsilon, 2 * r);
  CHECK(r1[0][0] == Catch::Approx(2));
  CHECK(r1[1][0] == Catch::Approx(11.5));
  std::vector<std::vector<Point2>> rings{{{0, 0}, {13, 0}, {13, 10}, {0, 10}}};
  auto r2 = narrow_to_border(r1, rings, 4 * r);
  CHECK(r2[0][0] == Catch::Approx(0).margin(1e-9));
  CHECK(r2[1][0] == Catch::Approx(13));
  auto r3 = narrow_to_border(r0, rings, 4 * r);
  CHECK(r3[0][0] == Catch::Approx(4));
  CHECK(r3[1][0] == Catch::Approx(8));
}

TEST_CASE("a structure beside the border keeps its place") {
  Rectangle r0{Point2{4, 8.5}, Point2{8, 8.5}, Point2{8, 8.7}, Point2{4, 8.7}};
  std::vector<std::vector<Point2>> rings{{{0, 0}, {13, 0}, {13, 10}, {0, 10}}};
  auto r1 = narrow_to_border(r0, rings, 4 * r);
  for (int i = 0; i < 4; ++i) {
    CHECK(r1[i][0] == Catch::Approx(r0[i][0]));
    CHECK(r1[i][1] == Catch::Approx(r0[i][1]));
  }
}

TEST_CASE("a short structure's lines run on no farther than it is long") {
  Rectangle r0{Point2{4, 4.9}, Point2{6, 4.9}, Point2{6, 5.1}, Point2{4, 5.1}};
  auto edges = narrow_edges(r0, 10);
  CHECK(edges[0][0][0] == Catch::Approx(2));
  CHECK(edges[0][1][0] == Catch::Approx(8));
}

TEST_CASE("nothing stands on a steep roof plane") {
  std::vector<std::array<double, 3>> pts;
  std::vector<int> labels;
  roof([](double x, double) { return std::abs(x - 8) < 0.2 ? 0.6 : 0.0; }, pts,
       labels);
  auto stands = standing(pts, labels, 2 * r, epsilon, {1});
  for (size_t i = 0; i < pts.size(); ++i) CHECK_FALSE(stands[i]);
}

TEST_CASE("points on one line make a rectangle with a width") {
  std::vector<Point2> xy;
  for (double x = 0; x < 5; x += 0.25) xy.push_back({x, 2});
  auto o = narrow_outline(xy, r);
  REQUIRE(o);
  CHECK(std::hypot((*o)[3][0] - (*o)[0][0], (*o)[3][1] - (*o)[0][1]) ==
        Catch::Approx(r / 10));
}

TEST_CASE("a line on a ridge's flank is no wall, a wall on a flat roof is") {
  Rectangle line{Point2{2, 4.9}, Point2{18, 4.9}, Point2{18, 5.1},
                 Point2{2, 5.1}};
  std::vector<std::array<double, 3>> flank, wall;
  for (double x = 0.125; x < 20; x += 0.25)
    for (double y = 0.125; y < 10; y += 0.25) {
      flank.push_back({x, y, 10 + 0.8 * y});
      wall.push_back({x, y, std::abs(y - 5) < 0.2 ? 10.6 : 10.0});
    }
  CHECK(narrow_on_flank(narrow_profile(flank, line, r), epsilon));
  CHECK_FALSE(narrow_on_flank(narrow_profile(wall, line, r), epsilon));
  auto m = narrow_profile(wall, line, r);
  CHECK(m[1] - std::fmax(m[0], m[2]) == Catch::Approx(0.6));
}
