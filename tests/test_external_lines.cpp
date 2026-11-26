#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <roofer/io/LineIndex.hpp>
#include <roofer/reconstruction/AlphaShaper.hpp>
#include <roofer/reconstruction/ArrangementBuilder.hpp>
#include <roofer/reconstruction/LineDetector.hpp>
#include <roofer/reconstruction/LineRegulariser.hpp>
#include <roofer/reconstruction/PlaneDetector.hpp>
#include <roofer/reconstruction/PlaneIntersector.hpp>

#include <filesystem>
#include <fstream>

using namespace roofer;

namespace {
  struct Source {
    std::filesystem::path path =
        std::filesystem::temp_directory_path() / "roofer-line-reader.geojson";
    explicit Source(std::string_view text) { std::ofstream(path) << text; }
    ~Source() {
      std::error_code error;
      std::filesystem::remove(path, error);
    }
  };

  LinearRing square(float size) {
    LinearRing ring;
    for (auto [x, y] :
         {std::pair{0.F, 0.F}, {size, 0.F}, {size, size}, {0.F, size}})
      ring.push_back({x, y, 0});
    return ring;
  }

  // Total length of the arrangement edges on the vertical line at x.
  double edge_length_at_x(const Arrangement_2& arrangement, double x) {
    double length = 0;
    for (auto e = arrangement.edges_begin(); e != arrangement.edges_end();
         ++e) {
      const auto& s = e->source()->point();
      const auto& t = e->target()->point();
      if (std::abs(CGAL::to_double(s.x()) - x) < 1e-3 &&
          std::abs(CGAL::to_double(t.x()) - x) < 1e-3)
        length += std::abs(CGAL::to_double(t.y() - s.y()));
    }
    return length;
  }
}  // namespace

TEST_CASE("line source is clipped to the footprint in local coordinates",
          "[external_lines]") {
  Source source(R"({"type":"FeatureCollection","features":[
    {"type":"Feature","properties":{},"geometry":{"type":"LineString",
      "coordinates":[[995,2005],[1015,2005]]}},
    {"type":"Feature","properties":{},"geometry":{"type":"MultiLineString",
      "coordinates":[[[1002,2001],[1002,2003]],[[1100,2100],[1101,2101]]]}},
    {"type":"Feature","properties":{},"geometry":{"type":"Point",
      "coordinates":[1005,2005]}},
    {"type":"Feature","properties":{},"geometry":{"type":"LineString",
      "coordinates":[[1200,2000],[1210,2000]]}}
  ]})");
  io::LineIndex index(source.path.string(), "");
  CHECK(index.size() == 4);

  auto footprint = square(10);
  footprint.interior_rings().push_back(
      {{4, 4, 0}, {4, 6, 0}, {6, 6, 0}, {6, 4, 0}});
  const auto lines = index.clip(footprint, {1000, 2000, 100});
  REQUIRE(lines.size() == 3);
  CHECK(lines[0] == std::array{arr3f{0, 5, 0}, arr3f{4, 5, 0}});
  CHECK(lines[1] == std::array{arr3f{6, 5, 0}, arr3f{10, 5, 0}});
  CHECK(lines[2] == std::array{arr3f{2, 1, 0}, arr3f{2, 3, 0}});

  CHECK_THROWS(io::LineIndex(source.path.string(), "missing"));
  CHECK_THROWS(io::LineIndex(source.path.string(), "", "EPSG:28992"));
}

TEST_CASE("an external line becomes an arrangement edge", "[external_lines]") {
  using namespace roofer::reconstruction;
  // a flat 10 x 10 m roof: one plane, no plane intersections, and nothing in
  // the point cloud along x = 5
  PointCollection points;
  for (float x = 0.25F; x < 10; x += 0.5F)
    for (float y = 0.25F; y < 10; y += 0.5F) points.push_back({x, y, 3});
  auto footprint = square(10);

  auto planes = createPlaneDetector();
  planes->detect(points);
  REQUIRE(planes->pts_per_roofplane.size() == 1);
  auto alpha_shaper = createAlphaShaper();
  alpha_shaper->compute(planes->pts_per_roofplane);
  auto line_detector = createLineDetector();
  line_detector->detect(alpha_shaper->alpha_rings, alpha_shaper->roofplane_ids,
                        planes->pts_per_roofplane);
  auto plane_intersector = createPlaneIntersector();
  plane_intersector->compute(planes->pts_per_roofplane,
                             planes->plane_adjacencies);

  auto edge_length = [&](const SegmentCollection& external, int priority) {
    LineRegulariserConfig config;
    config.external_line_priority = priority;
    auto regulariser = createLineRegulariser();
    regulariser->compute(line_detector->edge_segments,
                         plane_intersector->segments, config, external);
    Arrangement_2 arrangement;
    createArrangementBuilder()->compute(arrangement, footprint,
                                        regulariser->exact_regularised_edges);
    return edge_length_at_x(arrangement, 5);
  };

  SegmentCollection external;
  CHECK(edge_length(external, 0) == 0);
  external.push_back({arr3f{5, 0, 0}, arr3f{5, 10, 0}});
  for (int priority : {0, 1, 2}) {
    CAPTURE(priority);
    CHECK(edge_length(external, priority) == Catch::Approx(10).margin(1e-3));
  }
}

namespace {
  LinearRing ring(std::initializer_list<std::pair<float, float>> points) {
    LinearRing result;
    for (auto [x, y] : points) result.push_back({x, y, 0});
    return result;
  }

  bool has_arc(const SegmentCollection& arcs, arr3f a, arr3f b) {
    auto near = [](const arr3f& p, const arr3f& q) {
      return std::abs(p[0] - q[0]) < 1e-3 && std::abs(p[1] - q[1]) < 1e-3;
    };
    for (const auto& s : arcs)
      if ((near(s[0], a) && near(s[1], b)) || (near(s[0], b) && near(s[1], a)))
        return true;
    return false;
  }
}  // namespace

TEST_CASE("footprint skeleton of a rectangle is a hipped roof",
          "[external_lines]") {
  auto arcs = reconstruction::footprint_skeleton(
      ring({{0, 0}, {10, 0}, {10, 4}, {0, 4}}));
  CHECK(arcs.size() == 5);
  CHECK(has_arc(arcs, {2, 2, 0}, {8, 2, 0}));
  CHECK(has_arc(arcs, {0, 0, 0}, {2, 2, 0}));
  CHECK(has_arc(arcs, {10, 4, 0}, {8, 2, 0}));

  // clockwise and closed gives the same arcs
  auto cw = reconstruction::footprint_skeleton(
      ring({{0, 0}, {0, 4}, {10, 4}, {10, 0}, {0, 0}}));
  CHECK(cw.size() == 5);
  CHECK(has_arc(cw, {2, 2, 0}, {8, 2, 0}));

  // the hips are 2.83 m, so only the ridge is at least 3 m
  auto long_arcs = reconstruction::footprint_skeleton(
      ring({{0, 0}, {10, 0}, {10, 4}, {0, 4}}), 3.F);
  CHECK(long_arcs.size() == 1);
  CHECK(has_arc(long_arcs, {2, 2, 0}, {8, 2, 0}));
}

TEST_CASE("footprint skeleton of an L-shape joins both ridges",
          "[external_lines]") {
  auto arcs = reconstruction::footprint_skeleton(
      ring({{0, 0}, {10, 0}, {10, 4}, {4, 4}, {4, 10}, {0, 10}}));
  CHECK(arcs.size() == 8);
  CHECK(has_arc(arcs, {2, 2, 0}, {8, 2, 0}));
  CHECK(has_arc(arcs, {2, 2, 0}, {2, 8, 0}));
  CHECK(has_arc(arcs, {4, 4, 0}, {2, 2, 0}));
}

TEST_CASE("footprint skeleton respects holes and rejects invalid rings",
          "[external_lines]") {
  auto frame = ring({{0, 0}, {10, 0}, {10, 10}, {0, 10}});
  frame.interior_rings().push_back(ring({{3, 3}, {7, 3}, {7, 7}, {3, 7}}));
  auto arcs = reconstruction::footprint_skeleton(frame);
  // a ring of corner arcs and the frame's centre line
  CHECK(arcs.size() == 12);
  CHECK(has_arc(arcs, {0, 0, 0}, {1.5, 1.5, 0}));
  CHECK(has_arc(arcs, {3, 3, 0}, {1.5, 1.5, 0}));
  CHECK(has_arc(arcs, {1.5, 1.5, 0}, {8.5, 1.5, 0}));

  auto bow_tie = ring({{0, 0}, {10, 10}, {10, 0}, {0, 10}});
  CHECK(reconstruction::footprint_skeleton(bow_tie).empty());
}
