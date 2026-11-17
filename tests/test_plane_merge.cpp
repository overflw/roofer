#include <catch2/catch_test_macros.hpp>
#include <roofer/reconstruction/PlaneDetectorBase.hpp>

#include <cmath>
#include <map>
#include <numbers>
#include <random>
#include <vector>

using namespace roofer;

namespace {
  struct Fixture {
    PointCollection points;
    vec3f normals;
    std::vector<size_t> region_ids;
    std::vector<planedect::PlaneRegion> regions;
    std::map<size_t, std::map<size_t, size_t>> adjacencies;

    // add a region of points z = slope * x + offset on [x0, x1) x [0, 10)
    void add_region(size_t id, float x0, float x1, float slope, float offset,
                    std::mt19937& rng) {
      std::normal_distribution<float> noise(0.0F, 0.03F);
      planedect::PlaneRegion region(id);
      for (float x = x0; x < x1; x += 0.25F) {
        for (float y = 0; y < 10; y += 0.25F) {
          region.inliers.push_back(points.size());
          region_ids.push_back(id);
          points.push_back({x, y, slope * x + offset + noise(rng)});
          normals.push_back({0, 0, 1});
        }
      }
      regions.push_back(region);
    }
  };

  // fit each region's plane, as the region grower does
  void fit_planes(planedect::PlaneDS& pds,
                  std::vector<planedect::PlaneRegion>& regions) {
    for (auto& region : regions) {
      if (region.get_region_id() != 0)
        pds.fit_plane(region.inliers, region.plane);
    }
  }

  planedect::RegionMergeParams default_params() {
    return {10.0 * std::numbers::pi / 180.0, 1.2, 0.02};
  }
}  // namespace

TEST_CASE("adjacent coplanar regions are merged", "[planedetector]") {
  std::mt19937 rng(1);
  Fixture f;
  f.regions.emplace_back(0);
  // one plane split in two regions, plus a steeper plane next to it
  f.add_region(1, 0, 5, 0.6F, 0, rng);
  f.add_region(2, 5, 10, 0.6F, 0, rng);
  f.add_region(3, 10, 15, -0.6F, 12, rng);
  f.adjacencies[2][1] = 10;
  f.adjacencies[3][2] = 10;
  planedect::PlaneDS pds(f.points, f.normals);
  fit_planes(pds, f.regions);

  auto merges = planedect::merge_coplanar_regions(
      pds, f.regions, f.adjacencies, f.region_ids, default_params());
  CHECK(merges == 1);
  REQUIRE(f.regions.size() == 3);
  CHECK(f.regions[1].get_region_id() == 1);
  CHECK(f.regions[1].inliers.size() == 2 * 20 * 40);
  CHECK(f.regions[2].get_region_id() == 3);
  CHECK(f.region_ids[20 * 40] == 1);
}

TEST_CASE("parallel regions with a step are not merged", "[planedetector]") {
  std::mt19937 rng(2);
  Fixture f;
  f.regions.emplace_back(0);
  f.add_region(1, 0, 5, 0.6F, 0, rng);
  // same slope, 0.4 m higher: a step in the roof
  f.add_region(2, 5, 10, 0.6F, 0.4F, rng);
  f.adjacencies[2][1] = 10;
  planedect::PlaneDS pds(f.points, f.normals);
  fit_planes(pds, f.regions);

  auto merges = planedect::merge_coplanar_regions(
      pds, f.regions, f.adjacencies, f.region_ids, default_params());
  CHECK(merges == 0);
  CHECK(f.regions.size() == 3);
}

TEST_CASE("coplanar regions that are not adjacent are not merged",
          "[planedetector]") {
  std::mt19937 rng(3);
  Fixture f;
  f.regions.emplace_back(0);
  f.add_region(1, 0, 5, 0.6F, 0, rng);
  f.add_region(2, 20, 25, 0.6F, 0, rng);
  planedect::PlaneDS pds(f.points, f.normals);
  fit_planes(pds, f.regions);

  auto merges = planedect::merge_coplanar_regions(
      pds, f.regions, f.adjacencies, f.region_ids, default_params());
  CHECK(merges == 0);
}
