#include <roofer/reconstruction/WallTopDetector.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <random>

using namespace roofer::reconstruction;
namespace {
  // A flat roof 20 x 10 m at 10 m, 16 points/m², segmented as plane 1 but
  // where `extra` raises a point; a raised point stays unsegmented.
  template <typename F>
  std::vector<WallTop> detect(
      F extra, std::vector<std::vector<std::array<double, 2>>> rings = {}) {
    std::vector<std::array<float, 3>> pts;
    std::vector<int> labels;
    std::mt19937 rng(0);
    std::normal_distribution<float> noise(0, 0.03F);
    for (float x = 0.125F; x < 20; x += 0.25F)
      for (float y = 0.125F; y < 10; y += 0.25F) {
        float h = extra(x, y);
        pts.push_back({x, y, 10 + h + noise(rng)});
        labels.push_back(h > 0.2F ? 0 : 1);
      }
    return detect_wall_tops(pts, labels, {{1, {0, 0, 1, -10}}}, rings);
  }
}  // namespace

TEST_CASE("a wall across a roof is found, on top of the wall") {
  auto walls = detect([](float x, float y) {
    return std::abs(x - 8 - 0.2F * y) < 0.2F ? 0.6F : 0.0F;
  });
  REQUIRE(walls.size() == 1);
  const auto& w = walls[0];
  CHECK(std::hypot(w.x1 - w.x0, w.y1 - w.y0) > 8);
  CHECK(-(w.a * 10 + w.b * 5 + w.d) / w.c == Catch::Approx(10.6).margin(0.05));
}

TEST_CASE("PV rows and scattered points are no walls") {
  CHECK(detect([](float x, float y) {
          float m = std::fmod(x, 3.0F);
          return m < 1.6F ? 0.3F + 0.4F * m : 0.0F;
        }).empty());
  std::mt19937 rng(1);
  std::uniform_real_distribution<float> u(0, 1);
  CHECK(detect([&](float, float) {
          return u(rng) < 0.05F ? 0.5F : 0.0F;
        }).empty());
}

TEST_CASE("a broad raised band along the border is no wall") {
  CHECK(detect([](float, float y) { return y > 8.2F ? 0.8F : 0.0F; },
               {{{0, 0}, {20, 0}, {20, 10}, {0, 10}}})
            .empty());
}

TEST_CASE("a wall ending short of the border runs on to it") {
  WallTop w{{}, 0, 0, 1, -10, 1.2, 5, 18, 5, 0.5};
  snap_wall_top(w, {{{0, 0}, {20, 0}, {20, 10}, {0, 10}}}, 2.5);
  CHECK(w.x0 == Catch::Approx(-0.5));
  CHECK(w.x1 == Catch::Approx(20.5));
  WallTop free{{}, 0, 0, 1, -10, 5, 5, 12, 5, 0.5};
  snap_wall_top(free, {{{0, 0}, {20, 0}, {20, 10}, {0, 10}}}, 2.5);
  CHECK(free.x0 == Catch::Approx(5));
}
