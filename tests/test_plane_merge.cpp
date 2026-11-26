#include <catch2/catch_test_macros.hpp>
#include <roofer/reconstruction/PlaneDetectorBase.hpp>

#include <map>
#include <random>
#include <vector>

using namespace roofer;

namespace {
  struct Fixture {
    PointCollection points;
    vec3f normals;
    std::vector<planedect::PlaneRegion> regions{planedect::PlaneRegion(0)};
    std::map<size_t, std::map<size_t, size_t>> adjacencies;
    std::mt19937 rng{1};

    // a region of points z = z0 + slope * (x - x0) on [x0, x1) x [0, 10)
    void add_region(size_t id, float x0, float x1, float slope, float z0) {
      std::normal_distribution<float> noise(0.0F, 0.03F);
      planedect::PlaneRegion region(id);
      for (float x = x0; x < x1; x += 0.25F) {
        for (float y = 0; y < 10; y += 0.25F) {
          region.inliers.push_back(points.size());
          points.push_back({x, y, z0 + slope * (x - x0) + noise(rng)});
          normals.push_back({0, 0, 1});
        }
      }
      regions.push_back(region);
    }

    size_t merge() {
      planedect::PlaneDS pds(points, normals);
      for (auto& region : regions)
        if (region.get_region_id() != 0)
          pds.fit_plane(region.inliers, region.plane);
      return planedect::merge_coplanar_regions(pds, regions, adjacencies, 20,
                                               0.15, 0.02);
    }
  };
}  // namespace

TEST_CASE("a plane split in two is merged, not its steeper neighbour",
          "[planedetector]") {
  Fixture f;
  f.add_region(1, 0, 5, 0.6F, 0);
  f.add_region(2, 5, 10, 0.6F, 3);
  f.add_region(3, 10, 15, 0.4F, 6);
  f.adjacencies[2][1] = 10;
  f.adjacencies[3][2] = 10;
  CHECK(f.merge() == 1);
  REQUIRE(f.regions.size() == 3);
  CHECK(f.regions[1].get_region_id() == 1);
  CHECK(f.regions[1].inliers.size() == 2 * 20 * 40);
  CHECK(f.regions[2].get_region_id() == 3);
}

TEST_CASE("a slope with a slight kink stays two planes", "[planedetector]") {
  Fixture f;
  f.add_region(1, 0, 5, 0.6F, 0);
  f.add_region(2, 5, 10, 0.5F, 3);
  f.adjacencies[2][1] = 10;
  CHECK(f.merge() == 0);
  CHECK(f.regions.size() == 3);
}

TEST_CASE("parallel planes with a step stay apart", "[planedetector]") {
  Fixture f;
  f.add_region(1, 0, 5, 0.6F, 0);
  f.add_region(2, 5, 10, 0.6F, 3.4F);
  f.adjacencies[2][1] = 10;
  CHECK(f.merge() == 0);
}

TEST_CASE("coplanar planes that don't touch stay apart", "[planedetector]") {
  Fixture f;
  f.add_region(1, 0, 5, 0.6F, 0);
  f.add_region(2, 20, 25, 0.6F, 12);
  CHECK(f.merge() == 0);
}

TEST_CASE("no regions when region growing was skipped", "[planedetector]") {
  Fixture f;
  f.regions.clear();
  CHECK(f.merge() == 0);
  CHECK(f.regions.empty());
}
