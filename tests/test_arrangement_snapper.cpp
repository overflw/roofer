#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <variant>

#include <catch2/catch_test_macros.hpp>

#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/centroid.h>

#include <roofer/reconstruction/ArrangementBase.hpp>
#include <roofer/reconstruction/ArrangementSnapper.hpp>
#include <roofer/reconstruction/ArrangementExtruder.hpp>

namespace {
  using roofer::Arrangement_2;
  using roofer::Point_2;
  using roofer::Segment_2;
  using Face_const_handle = Arrangement_2::Face_const_handle;

  Arrangement_2 cross_arrangement(const std::array<double, 4>& heights) {
    Arrangement_2 arrangement;
    const std::array<Point_2, 4> corners = {Point_2(0, 0), Point_2(10, 0),
                                            Point_2(10, 10), Point_2(0, 10)};
    for (std::size_t i = 0; i < corners.size(); ++i) {
      CGAL::insert(arrangement,
                   Segment_2(corners[i], corners[(i + 1) % corners.size()]));
    }
    CGAL::insert(arrangement, Segment_2(Point_2(0, 5), Point_2(10, 5)));
    CGAL::insert(arrangement, Segment_2(Point_2(5, 0), Point_2(5, 10)));

    using PointLocation =
        CGAL::Arr_walk_along_line_point_location<Arrangement_2>;
    PointLocation point_location(arrangement);
    const std::array<Point_2, 4> samples = {
        Point_2(7.5, 7.5), Point_2(2.5, 7.5), Point_2(2.5, 2.5),
        Point_2(7.5, 2.5)};
    for (std::size_t i = 0; i < samples.size(); ++i) {
      auto object = point_location.locate(samples[i]);
      auto face = std::get_if<Face_const_handle>(&object);
      REQUIRE(face != nullptr);
      auto mutable_face = arrangement.non_const_handle(*face);
      mutable_face->data().in_footprint = true;
      mutable_face->data().segid = static_cast<int>(i + 1);
      mutable_face->data().plane = roofer::Plane(0, 0, 1, -heights[i]);
    }
    return arrangement;
  }

  Arrangement_2::Face_handle face_at(Arrangement_2& arrangement,
                                     const Point_2& point) {
    using PointLocation =
        CGAL::Arr_walk_along_line_point_location<Arrangement_2>;
    PointLocation point_location(arrangement);
    auto object = point_location.locate(point);
    auto face = std::get_if<Face_const_handle>(&object);
    REQUIRE(face != nullptr);
    return arrangement.non_const_handle(*face);
  }

  Arrangement_2 repeated_face_arrangement() {
    auto arrangement = cross_arrangement({1, 2, 1, 3});
    auto repeated_face = face_at(arrangement, Point_2(7.5, 7.5));
    auto other_face = face_at(arrangement, Point_2(2.5, 2.5));
    other_face->data().segid = repeated_face->data().segid;
    other_face->data().plane = repeated_face->data().plane;
    return arrangement;
  }

  Arrangement_2 boundary_junction_arrangement() {
    Arrangement_2 arrangement;
    const std::array<Point_2, 4> corners = {Point_2(0, 0), Point_2(10, 0),
                                            Point_2(10, 10), Point_2(0, 10)};
    for (std::size_t i = 0; i < corners.size(); ++i) {
      CGAL::insert(arrangement,
                   Segment_2(corners[i], corners[(i + 1) % corners.size()]));
    }
    CGAL::insert(arrangement, Segment_2(Point_2(5, 0), Point_2(2, 10)));
    CGAL::insert(arrangement, Segment_2(Point_2(5, 0), Point_2(8, 10)));

    auto left = face_at(arrangement, Point_2(1, 5));
    left->data().in_footprint = true;
    left->data().segid = 1;
    left->data().plane = roofer::Plane(0, 0, 1, -15);

    auto middle = face_at(arrangement, Point_2(5, 5));
    middle->data().in_footprint = true;
    middle->data().segid = 2;
    middle->data().plane = roofer::Plane(0, 0, 1, -5);

    auto right = face_at(arrangement, Point_2(9, 5));
    right->data().in_footprint = true;
    right->data().segid = 3;
    right->data().plane = roofer::Plane(0, 0, 1, -15);

    arrangement.unbounded_face()->data().in_footprint = false;
    arrangement.unbounded_face()->data().segid = 0;
    // Exterior plane data is not used by extrusion and may retain an
    // unrelated roof plane. The snapper must use the ground provider instead.
    arrangement.unbounded_face()->data().plane = roofer::Plane(0, 0, 1, -100);

    return arrangement;
  }

  Arrangement_2 finite_exterior_junction_arrangement() {
    auto arrangement = cross_arrangement({15, 5, 15, 0});

    auto hole = face_at(arrangement, Point_2(7.5, 2.5));
    hole->data().in_footprint = false;
    hole->data().is_footprint_hole = true;
    hole->data().segid = 0;
    hole->data().plane = roofer::Plane(0, 0, 1, -100);
    return arrangement;
  }

  Arrangement_2 non_footprint_bridge_arrangement() {
    auto arrangement = cross_arrangement({0, 5, 0, 15});

    for (const auto& sample : {Point_2(7.5, 7.5), Point_2(2.5, 2.5)}) {
      auto hole = face_at(arrangement, sample);
      hole->data().in_footprint = false;
      hole->data().is_footprint_hole = true;
      hole->data().segid = 0;
      // Exterior plane data is intentionally stale; non-footprint sectors use
      // the exterior/base elevation provider instead.
      hole->data().plane = roofer::Plane(0, 0, 1, -100);
    }

    return arrangement;
  }

  Arrangement_2 narrow_angle_junction_arrangement() {
    Arrangement_2 arrangement;
    const std::array<Point_2, 4> corners = {Point_2(0, 0), Point_2(10, 0),
                                            Point_2(10, 10), Point_2(0, 10)};
    for (std::size_t i = 0; i < corners.size(); ++i) {
      CGAL::insert(arrangement,
                   Segment_2(corners[i], corners[(i + 1) % corners.size()]));
    }
    CGAL::insert(arrangement, Segment_2(Point_2(5, 5), Point_2(10, 5)));
    CGAL::insert(arrangement, Segment_2(Point_2(5, 5), Point_2(10, 5.1)));
    CGAL::insert(arrangement, Segment_2(Point_2(5, 5), Point_2(5, 10)));
    CGAL::insert(arrangement, Segment_2(Point_2(5, 5), Point_2(0, 5)));
    CGAL::insert(arrangement, Segment_2(Point_2(5, 5), Point_2(5, 0)));

    auto set_roof = [&](const Point_2& sample, int segid, double height) {
      auto face = face_at(arrangement, sample);
      face->data().in_footprint = true;
      face->data().segid = segid;
      face->data().plane = roofer::Plane(0, 0, 1, -height);
    };

    set_roof(Point_2(7.5, 2.5), 1, 6);
    set_roof(Point_2(7.5, 5.025), 2, 15);
    set_roof(Point_2(7.5, 7.5), 3, 15);
    set_roof(Point_2(2.5, 7.5), 4, 5);
    set_roof(Point_2(2.5, 2.5), 5, 15);

    arrangement.unbounded_face()->data().in_footprint = false;
    arrangement.unbounded_face()->data().segid = 0;
    arrangement.unbounded_face()->data().plane = roofer::Plane(0, 0, 1, 0);

    return arrangement;
  }

  Arrangement_2 clearance_regression_arrangement() {
    auto arrangement = cross_arrangement({15, 5, 15, 6});
    CGAL::insert(arrangement, Segment_2(Point_2(5.05, 5.55), Point_2(5.05, 9)));
    CGAL::insert(arrangement, Segment_2(Point_2(5.55, 5.05), Point_2(9, 5.05)));
    return arrangement;
  }

  Arrangement_2 gapped_split_arrangement() {
    Arrangement_2 arrangement;
    const std::array<Point_2, 4> corners = {Point_2(0, 0), Point_2(10, 0),
                                            Point_2(10, 10), Point_2(0, 10)};
    for (std::size_t i = 0; i < corners.size(); ++i) {
      CGAL::insert(arrangement,
                   Segment_2(corners[i], corners[(i + 1) % corners.size()]));
    }
    CGAL::insert(arrangement, Segment_2(Point_2(0, 5), Point_2(5, 5)));
    CGAL::insert(arrangement, Segment_2(Point_2(5.01, 5), Point_2(10, 5)));

    auto roof = face_at(arrangement, Point_2(5, 7.5));
    roof->data().in_footprint = true;
    roof->data().segid = 1;
    roof->data().plane = roofer::Plane(0, 0, 1, -1);

    arrangement.unbounded_face()->data().in_footprint = false;
    arrangement.unbounded_face()->data().segid = 0;
    arrangement.unbounded_face()->data().plane = roofer::Plane(0, 0, 1, 0);

    return arrangement;
  }

  Arrangement_2 sliver_hole_arrangement(
      const std::array<Point_2, 3>& hole_corners, const Point_2& inside_hole,
      bool roof_first) {
    Arrangement_2 arrangement;
    const std::array<Point_2, 4> corners = {Point_2(0, 0), Point_2(10, 0),
                                            Point_2(10, 10), Point_2(0, 10)};
    auto insert_ring = [&](const auto& points) {
      for (std::size_t i = 0; i < points.size(); ++i) {
        CGAL::insert(arrangement,
                     Segment_2(points[i], points[(i + 1) % points.size()]));
      }
    };
    if (roof_first) {
      insert_ring(corners);
      insert_ring(hole_corners);
    } else {
      insert_ring(hole_corners);
      insert_ring(corners);
    }

    auto roof = face_at(arrangement, Point_2(1, 1));
    roof->data().in_footprint = true;
    roof->data().segid = 1;
    roof->data().plane = roofer::Plane(0, 0, 1, -1);

    auto hole = face_at(arrangement, inside_hole);
    hole->data().in_footprint = true;
    hole->data().segid = 2;
    hole->data().plane = roofer::Plane(0, 0, 1, -2);

    arrangement.unbounded_face()->data().in_footprint = false;
    arrangement.unbounded_face()->data().segid = 0;
    arrangement.unbounded_face()->data().plane = roofer::Plane(0, 0, 1, 0);
    return arrangement;
  }

  Arrangement_2::Vertex_handle vertex_at(Arrangement_2& arrangement,
                                         const Point_2& point) {
    for (auto vertex : arrangement.vertex_handles()) {
      if (vertex->point() == point) return vertex;
    }
    return {};
  }

  bool has_vertex_near(Arrangement_2& arrangement, const Point_2& point,
                       double tolerance) {
    const double tolerance_sq = tolerance * tolerance;
    for (auto vertex : arrangement.vertex_handles()) {
      if (CGAL::to_double(CGAL::squared_distance(vertex->point(), point)) <=
          tolerance_sq) {
        return true;
      }
    }
    return false;
  }

  std::size_t count_vertices_near(Arrangement_2& arrangement,
                                  const Point_2& point, double tolerance) {
    std::size_t count = 0;
    const double tolerance_sq = tolerance * tolerance;
    for (auto vertex : arrangement.vertex_handles()) {
      if (CGAL::to_double(CGAL::squared_distance(vertex->point(), point)) <=
          tolerance_sq) {
        ++count;
      }
    }
    return count;
  }

  void check_two_manifold_edges(const roofer::Mesh& mesh) {
    using Edge = std::pair<roofer::arr3f, roofer::arr3f>;
    std::map<Edge, std::size_t> edge_incidence;

    auto add_ring = [&](const roofer::vec3f& ring) {
      for (std::size_t i = 0; i < ring.size(); ++i) {
        auto first = ring[i];
        auto second = ring[(i + 1) % ring.size()];
        if (second < first) std::swap(first, second);
        if (first != second) ++edge_incidence[{first, second}];
      }
    };
    for (const auto& polygon : mesh.get_polygons()) {
      add_ring(polygon);
      for (const auto& hole : polygon.interior_rings()) add_ring(hole);
    }
    for (const auto& [edge, incidence] : edge_incidence) {
      CAPTURE(edge.first, edge.second);
      CHECK(incidence == 2);
    }
  }

  // A roofer building's arrangement; the snapper reroutes the edge between
  // face 3 and an exterior sliver past two close vertices.
  Arrangement_2 rerouted_sliver_arrangement() {
    constexpr std::array<std::array<double, 4>, 51> edges = {{
        {171.09100341796875, 211.20899963378906, 170.80099487304688,
         211.13299560546875},
        {171.11500549316406, 211.12399291992188, 171.09100341796875,
         211.20899963378906},
        {181.31399536132812, 215.76199340820312, 181.41099548339844,
         215.78900146484375},
        {169.69466947839328, 215.05152958469449, 169.85499572753906,
         214.48399353027344},
        {169.40926725042706, 214.09575427834744, 169.92399597167969,
         214.24099731445312},
        {166.3489990234375, 216.41900634765625, 168.58296450523352,
         217.04901212241106},
        {169.01028363580144, 217.47804531128992, 168.66299438476562,
         218.70899963378906},
        {176.25276267357802, 212.57307435569547, 174.13307281887404,
         220.29245660000817},
        {179.99887065630543, 220.56739859118119, 179.63800048828125,
         221.88600158691406},
        {170.80099487304688, 211.13299560546875, 170.30499381585525,
         212.89077675092551},
        {181.85800170898438, 214.15400695800781, 181.47303585331241,
         215.56207920245589},
        {180.59991330443484, 218.37121639987919, 180.34261801127951,
         218.29382780905883},
        {169.69466947839328, 215.05152958469449, 169.16267161032346,
         216.93473943095142},
        {180.39007320101359, 218.06166796833139, 180.34261801127951,
         218.29382780905883},
        {180.66484720334745, 218.13395091988798, 180.59991330443484,
         218.37121639987919},
        {180.66484720334745, 218.13395091988798, 180.39007320101359,
         218.06166796833139},
        {179.99887065630543, 220.56739859118119, 177.61015673791695,
         220.0768976694535},
        {181.47303585331241, 215.56207920245589, 181.41099548339844,
         215.78900146484375},
        {167.18099975585938, 213.46699523925781, 166.3489990234375,
         216.41900634765625},
        {169.16267161032346, 216.93473943095142, 169.08999633789062,
         217.19200134277344},
        {169.16267161032346, 216.93473943095142, 168.64210154487114,
         216.83765263316775},
        {169.68455252850327, 215.08812334992064, 169.01028363580144,
         217.47804531128992},
        {176.25276267357802, 212.57307435569547, 171.11500549316406,
         211.12399291992188},
        {170.89508822646258, 213.03192240926046, 170.30499381585525,
         212.89077675092551},
        {167.18099975585938, 213.46699523925781, 169.40926725042706,
         214.09575427834744},
        {170.30499381585525, 212.89077675092551, 169.92399597167969,
         214.24099731445312},
        {171.99234008789062, 212.59428405761719, 172.72092768067728,
         214.00649673247761},
        {171.99234008789062, 212.59428405761719, 170.89508822646258,
         213.03192240926046},
        {168.64210154487114, 216.83765263316775, 168.58296450523352,
         217.04901212241106},
        {169.08999633789062, 217.19200134277344, 168.58296450523352,
         217.04901212241106},
        {169.42393745905582, 218.0208195807383, 169.01028363580144,
         217.47804531128992},
        {171.79591438585612, 217.66236794330047, 171.58288399147875,
         217.79008561031338},
        {169.69466947839328, 215.05152958469449, 169.20652948246789,
         214.82035177847212},
        {169.85499572753906, 214.48399353027344, 169.68455252850327,
         215.08812334992064},
        {172.83095931086609, 216.03449117518716, 169.68455252850327,
         215.08812334992064},
        {172.83095931086609, 216.03449117518716, 172.54816854361667,
         216.65977955280778},
        {180.59991330443484, 218.37121639987919, 179.99887065630543,
         220.56739859118119},
        {169.40926725042706, 214.09575427834744, 169.20652948246789,
         214.82035177847212},
        {171.58288399147875, 217.79008561031338, 169.42393745905582,
         218.0208195807383},
        {181.31399536132812, 215.76199340820312, 180.66484720334745,
         218.13395091988798},
        {174.13307281887404, 220.29245660000817, 179.63800048828125,
         221.88600158691406},
        {177.29395297456506, 217.24719241227177, 180.39007320101359,
         218.06166796833139},
        {177.29395297456506, 217.24719241227177, 177.29395008670139,
         217.24719625968709},
        {179.00237153677926, 214.97111691707698, 177.29395297456506,
         217.24719241227177},
        {179.00237153677926, 214.97111691707698, 181.47303585331241,
         215.56207920245589},
        {172.54816854361667, 216.65977955280778, 171.79591438585612,
         217.66236794330047},
        {181.85800170898438, 214.15400695800781, 176.25276267357802,
         212.57307435569547},
        {177.61015673791695, 220.0768976694535, 177.29395008670139,
         217.24719625968709},
        {172.94528535007856, 215.78170094347524, 172.72092768067728,
         214.00649673247761},
        {172.94528535007856, 215.78170094347524, 172.83095931086609,
         216.03449117518716},
        {174.13307281887404, 220.29245660000817, 168.66299438476562,
         218.70899963378906},
    }};
    Arrangement_2 arrangement;
    for (const auto& e : edges)
      CGAL::insert(arrangement,
                   Segment_2(Point_2(e[0], e[1]), Point_2(e[2], e[3])));
    struct Face {
      double x, y;
      int segid;
      double a, b, c, d;
    };
    constexpr std::array<Face, 9> faces = {{
        {168.022, 215.735, 4, 0.060728655003513272, 0.013524297732812789,
         -0.99806268532206854, 24.718584634409524},
        {174.224, 214.894, 7, 0.68071572979596806, 0.18799908947352215,
         -0.70801302075983685, -128.68210469934277},
        {171.009, 216.347, 3, -0.10983477899695213, 0.3428764198590552,
         0.93293734089022695, -94.215307904259106},
        {180.501, 218.214, 5, -0.68133048076688507, -0.18600874905613654,
         -0.70794740006058099, 190.8833935412124},
        {176.564, 216.405, 5, -0.68133048076688507, -0.18600874905613654,
         -0.70794740006058099, 190.8833935412124},
        {168.873, 216.992, 3, -0.10983477899695213, 0.3428764198590552,
         0.93293734089022695, -94.215307904259106},
        {178.941, 219.224, 2, 0.093576370149018248, -0.35490137743154654,
         -0.93020883421247125, 99.505924606880626},
        {171.050, 214.362, 6, 0.10478944625439066, -0.37045503942239122,
         0.9229204926320862, 23.215258976046044},
        {179.474, 216.518, 1, 0.093668005252051004, -0.35689069581616334,
         0.92943818300732428, 21.963916399275568},
    }};
    for (const auto& f : faces) {
      auto face = face_at(arrangement, Point_2(f.x, f.y));
      face->data().in_footprint = true;
      face->data().segid = f.segid;
      face->data().plane = roofer::Plane(f.a, f.b, f.c, f.d);
    }
    arrangement.unbounded_face()->data().in_footprint = false;
    arrangement.unbounded_face()->data().segid = 0;
    return arrangement;
  }
}  // namespace

TEST_CASE("snapper repairs an alternating four-way roof junction") {
  // Heights in NE, NW, SW, SE order. The uniquely lowest NW face must absorb
  // the repair area.
  auto arrangement = cross_arrangement({15, 5, 15, 6});
  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(arrangement, {.distance_threshold = 0.001F,
                                 .repair_non_manifold_vertices = true,
                                 .manifold_repair_radius = 0.5F,
                                 .manifold_height_tolerance = 1e-4F});

  CHECK(vertex_at(arrangement, Point_2(5, 5)) ==
        Arrangement_2::Vertex_handle());
  CHECK(face_at(arrangement, Point_2(5, 5))->data().segid == 2);

  std::size_t roof_faces = 0;
  for (auto face : arrangement.face_handles()) {
    if (face->data().in_footprint) ++roof_faces;
  }
  CHECK(roof_faces == 4);
  for (auto vertex : arrangement.vertex_handles()) {
    CHECK(vertex->degree() <= 3);
  }

  auto extruder = roofer::reconstruction::createArrangementExtruder();
  extruder->compute(arrangement, 0.0F);
  REQUIRE(extruder->meshes.size() == 1);
  check_two_manifold_edges(extruder->meshes.front());
}

TEST_CASE("snapper removes dangling constraint chains before repair") {
  auto expected = cross_arrangement({15, 5, 15, 6});
  auto arrangement = cross_arrangement({15, 5, 15, 6});
  CGAL::insert(arrangement, Segment_2(Point_2(5, 5), Point_2(6, 5.5)));
  CGAL::insert(arrangement, Segment_2(Point_2(6, 5.5), Point_2(7, 6)));

  const roofer::reconstruction::ArrangementSnapperConfig config{
      .distance_threshold = 0.001F,
      .repair_non_manifold_vertices = true,
      .manifold_repair_radius = 0.5F,
      .manifold_height_tolerance = 1e-4F};
  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(expected, config);
  snapper->compute(arrangement, config);

  CHECK(arrangement.number_of_vertices() == expected.number_of_vertices());
  CHECK(arrangement.number_of_edges() == expected.number_of_edges());
  CHECK(arrangement.number_of_faces() == expected.number_of_faces());
  for (auto vertex : expected.vertex_handles()) {
    CHECK(vertex_at(arrangement, vertex->point()) !=
          Arrangement_2::Vertex_handle());
  }
  for (auto vertex : arrangement.vertex_handles()) {
    CAPTURE(vertex->point());
    CHECK(vertex_at(expected, vertex->point()) !=
          Arrangement_2::Vertex_handle());
  }
}

TEST_CASE("snapper preserves a courtyard when its narrow entrance collapses") {
  // The 1 cm entrance collapses to a bridge with a ring at each end.
  Arrangement_2 arrangement;
  const std::array boundary = {
      Point_2(0, 0),    Point_2(10, 0),   Point_2(10, 10), Point_2(0, 10),
      Point_2(0, 5.01), Point_2(3, 5.01), Point_2(3, 6),   Point_2(5, 6),
      Point_2(5, 4),    Point_2(3, 4),    Point_2(3, 5),   Point_2(0, 5)};
  for (std::size_t i = 0; i < boundary.size(); ++i) {
    CGAL::insert(arrangement,
                 Segment_2(boundary[i], boundary[(i + 1) % boundary.size()]));
  }
  auto roof = face_at(arrangement, Point_2(8, 5));
  roof->data().in_footprint = true;
  roof->data().segid = 1;
  roof->data().plane = roofer::Plane(0, 0, 1, -5);

  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(arrangement, {.distance_threshold = 0.02F});

  auto snapped_roof = face_at(arrangement, Point_2(8, 5));
  CHECK(snapped_roof->data().in_footprint);
  CHECK_FALSE(face_at(arrangement, Point_2(4, 5))->data().in_footprint);
  const auto polygon = roofer::reconstruction::arr_face_polygon(snapped_roof);
  CHECK(polygon.outer_boundary().is_simple());
  CHECK(polygon.outer_boundary().area() == 100);
  REQUIRE(polygon.number_of_holes() == 1);
  CHECK(polygon.holes_begin()->is_simple());
  CHECK(polygon.holes_begin()->area() == -4);
}

TEST_CASE("snapper ignores non-constrained triangulation edges for clearance") {
  auto arrangement = clearance_regression_arrangement();
  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(arrangement, {.distance_threshold = 0.001F,
                                 .repair_non_manifold_vertices = true,
                                 .manifold_repair_radius = 0.5F,
                                 .manifold_height_tolerance = 1e-4F});

  CHECK(vertex_at(arrangement, Point_2(5, 5)) ==
        Arrangement_2::Vertex_handle());
  CHECK(has_vertex_near(arrangement, Point_2(5.5, 5), 1e-9));
  CHECK(has_vertex_near(arrangement, Point_2(5, 5.5), 1e-9));
}

TEST_CASE("snapper staggers close repair-cell vertices") {
  auto arrangement = narrow_angle_junction_arrangement();
  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(arrangement, {.distance_threshold = 0.001F,
                                 .repair_non_manifold_vertices = true,
                                 .manifold_repair_radius = 0.5F,
                                 .manifold_height_tolerance = 1e-4F});

  CHECK(vertex_at(arrangement, Point_2(5, 5)) ==
        Arrangement_2::Vertex_handle());

  const double near_ray_length = std::sqrt(5.0 * 5.0 + 0.1 * 0.1);
  const Point_2 east_outer_split(5.5, 5);
  const Point_2 east_inner_split(5 + 0.5 * 0.45, 5);
  const Point_2 near_outer_split(5 + 0.5 * 5 / near_ray_length,
                                 5 + 0.5 * 0.1 / near_ray_length);
  const Point_2 near_inner_split(5 + 0.5 * 0.45 * 5 / near_ray_length,
                                 5 + 0.5 * 0.45 * 0.1 / near_ray_length);
  const bool east_outer_near_inner =
      has_vertex_near(arrangement, east_outer_split, 1e-9) &&
      has_vertex_near(arrangement, near_inner_split, 1e-9);
  const bool east_inner_near_outer =
      has_vertex_near(arrangement, east_inner_split, 1e-9) &&
      has_vertex_near(arrangement, near_outer_split, 1e-9);

  CHECK(count_vertices_near(arrangement, Point_2(5, 5), 0.55) == 5);
  CHECK((east_outer_near_inner || east_inner_near_outer));
  CHECK_FALSE((has_vertex_near(arrangement, east_outer_split, 1e-9) &&
               has_vertex_near(arrangement, near_outer_split, 1e-9)));
  if (east_outer_near_inner) {
    CHECK(CGAL::to_double(CGAL::squared_distance(
              east_outer_split, near_inner_split)) >= 0.25 * 0.25);
  }
  if (east_inner_near_outer) {
    CHECK(CGAL::to_double(CGAL::squared_distance(
              east_inner_split, near_outer_split)) >= 0.25 * 0.25);
  }

  auto extruder = roofer::reconstruction::createArrangementExtruder();
  extruder->compute(arrangement, 0.0F);
  REQUIRE(extruder->meshes.size() == 1);
  check_two_manifold_edges(extruder->meshes.front());
}

TEST_CASE("snapper collapses unconstrained short triangulation edges") {
  auto arrangement = gapped_split_arrangement();
  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(arrangement, {.distance_threshold = 0.02F,
                                 .repair_non_manifold_vertices = false});

  CHECK_FALSE(has_vertex_near(arrangement, Point_2(5, 5), 1e-9));
  CHECK_FALSE(has_vertex_near(arrangement, Point_2(5.01, 5), 1e-9));
  CHECK(has_vertex_near(arrangement, Point_2(5.005, 5), 1e-9));
  CHECK(face_at(arrangement, Point_2(5, 7.5))->data().segid == 1);
}

TEST_CASE("snapper keeps the roof label when a hole in it collapses") {
  const std::array central = {Point_2(3, 5), Point_2(7, 5), Point_2(7, 5.01)};
  const std::array near_edge = {Point_2(0.2, 9.75), Point_2(7, 9.75),
                                Point_2(7, 9.76)};
  for (const auto& corners : {central, near_edge}) {
    const auto inside_hole = CGAL::centroid(corners[0], corners[1], corners[2]);
    for (bool roof_first : {false, true}) {
      CAPTURE(inside_hole, roof_first);
      auto arrangement =
          sliver_hole_arrangement(corners, inside_hole, roof_first);
      auto snapper = roofer::reconstruction::createArrangementSnapper();
      snapper->compute(arrangement, {.distance_threshold = 0.1F,
                                     .repair_non_manifold_vertices = false});

      CHECK(arrangement.number_of_faces() == 2);
      CHECK(face_at(arrangement, Point_2(1, 1))->data().segid == 1);
      CHECK(face_at(arrangement, inside_hole)->data().segid == 1);
    }
  }
}

TEST_CASE("snapper preserves a manifold four-way roof junction") {
  auto arrangement = cross_arrangement({0, 1, 2, 3});
  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(arrangement, {.distance_threshold = 0.001F,
                                 .repair_non_manifold_vertices = true,
                                 .manifold_repair_radius = 0.5F,
                                 .manifold_height_tolerance = 1e-4F});

  auto center = vertex_at(arrangement, Point_2(5, 5));
  REQUIRE(center != Arrangement_2::Vertex_handle());
  CHECK(center->degree() == 4);
}

TEST_CASE("snapper repairs a repeated-face junction") {
  auto arrangement = repeated_face_arrangement();
  auto repeated_segid = face_at(arrangement, Point_2(7.5, 7.5))->data().segid;

  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(arrangement, {.distance_threshold = 0.001F,
                                 .repair_non_manifold_vertices = true,
                                 .manifold_repair_radius = 0.5F,
                                 .manifold_height_tolerance = 1e-4F});

  CHECK(vertex_at(arrangement, Point_2(5, 5)) ==
        Arrangement_2::Vertex_handle());
  CHECK(face_at(arrangement, Point_2(5, 5))->data().segid == repeated_segid);
}

TEST_CASE("snapper keeps a boundary repair cell inside the footprint") {
  auto arrangement = boundary_junction_arrangement();

  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(arrangement, 0.0F,
                   {.distance_threshold = 0.001F,
                    .repair_non_manifold_vertices = true,
                    .manifold_repair_radius = 0.5F,
                    .manifold_height_tolerance = 1e-4F});

  CHECK(vertex_at(arrangement, Point_2(5, 0)) ==
        Arrangement_2::Vertex_handle());
  CHECK(face_at(arrangement, Point_2(5, 0.05))->data().in_footprint);
}

TEST_CASE("snapper keeps a finite exterior repair cell inside the footprint") {
  auto arrangement = finite_exterior_junction_arrangement();

  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(arrangement, 0.0F,
                   {.distance_threshold = 0.001F,
                    .repair_non_manifold_vertices = true,
                    .manifold_repair_radius = 0.5F,
                    .manifold_height_tolerance = 1e-4F});

  CHECK(vertex_at(arrangement, Point_2(5, 5)) ==
        Arrangement_2::Vertex_handle());
  CHECK(face_at(arrangement, Point_2(5, 5))->data().segid == 2);
}

TEST_CASE("snapper assigns a bridge repair cell to non-footprint") {
  auto arrangement = non_footprint_bridge_arrangement();

  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(arrangement, 0.0F,
                   {.distance_threshold = 0.001F,
                    .repair_non_manifold_vertices = true,
                    .manifold_repair_radius = 0.5F,
                    .manifold_height_tolerance = 1e-4F});

  CHECK(vertex_at(arrangement, Point_2(5, 5)) ==
        Arrangement_2::Vertex_handle());
  auto repair_face = face_at(arrangement, Point_2(5, 5));
  CHECK_FALSE(repair_face->data().in_footprint);
  CHECK(repair_face->data().segid == 0);

  auto extruder = roofer::reconstruction::createArrangementExtruder();
  extruder->compute(arrangement, 0.0F);
  REQUIRE(!extruder->meshes.empty());
  for (const auto& mesh : extruder->meshes) check_two_manifold_edges(mesh);
}

TEST_CASE("snapper rejects a non-finite exterior height") {
  auto arrangement = boundary_junction_arrangement();
  auto snapper = roofer::reconstruction::createArrangementSnapper();

  CHECK_THROWS(snapper->compute(arrangement,
                                std::numeric_limits<float>::quiet_NaN(),
                                {.distance_threshold = 0.001F,
                                 .repair_non_manifold_vertices = true,
                                 .manifold_repair_radius = 0.5F,
                                 .manifold_height_tolerance = 1e-4F}));
}

TEST_CASE("snapper keeps the outer labels of the edges a collapse merges") {
  // Collapsing v1-v2 (or v1-v3-v2) merges the edges to p, which removes the
  // slivers between them. With q, it also rebuilds the triangle v1-q-p from
  // new faces only.
  for (bool mirrored : {false, true}) {
    for (bool with_q : {true, false}) {
      for (bool with_v3 : {false, true}) {
        CAPTURE(mirrored, with_q, with_v3);
        auto at = [&](double x, double y) {
          return Point_2(mirrored ? 10 - x : x, y);
        };
        const Point_2 v1 = at(5, 5), v2 = at(5.01, 5), v3 = at(5.005, 5.008),
                      p = at(5, 10), q = at(3, 8);
        Arrangement_2 arrangement;
        const std::array<Point_2, 4> corners = {
            Point_2(0, 0), Point_2(10, 0), Point_2(10, 10), Point_2(0, 10)};
        for (std::size_t i = 0; i < corners.size(); ++i) {
          CGAL::insert(
              arrangement,
              Segment_2(corners[i], corners[(i + 1) % corners.size()]));
        }
        for (const auto& segment :
             {Segment_2(v1, v2), Segment_2(v1, p), Segment_2(v2, p)}) {
          CGAL::insert(arrangement, segment);
        }
        if (with_v3) {
          CGAL::insert(arrangement, Segment_2(v1, v3));
          CGAL::insert(arrangement, Segment_2(v3, v2));
          CGAL::insert(arrangement, Segment_2(v3, p));
        }
        if (with_q) {
          CGAL::insert(arrangement, Segment_2(v1, q));
          CGAL::insert(arrangement, Segment_2(q, p));
        }
        const Point_2 inside_small = CGAL::centroid(v1, q, p);
        auto set_roof = [&](const Point_2& sample, int segid) {
          auto face = face_at(arrangement, sample);
          face->data().in_footprint = true;
          face->data().segid = segid;
          face->data().plane = roofer::Plane(0, 0, 1, -segid);
        };
        set_roof(Point_2(1, 1), 1);
        if (with_q) set_roof(inside_small, 2);
        set_roof(CGAL::centroid(v3, v2, p), 3);
        if (with_v3) {
          set_roof(CGAL::centroid(v1, v3, p), 4);
          set_roof(CGAL::centroid(v1, v2, v3), 5);
        }
        arrangement.unbounded_face()->data().in_footprint = false;
        arrangement.unbounded_face()->data().segid = 0;
        arrangement.unbounded_face()->data().plane = roofer::Plane(0, 0, 1, 0);

        auto snapper = roofer::reconstruction::createArrangementSnapper();
        snapper->compute(arrangement, {.distance_threshold = 0.02F,
                                       .repair_non_manifold_vertices = false});

        CHECK_FALSE(has_vertex_near(arrangement, v1, 1e-9));
        CHECK(arrangement.number_of_faces() == (with_q ? 3 : 2));
        CHECK(face_at(arrangement, inside_small)->data().segid ==
              (with_q ? 2 : 1));
      }
    }
  }
}

TEST_CASE("snapper keeps a roof face whose edge it reroutes") {
  auto arrangement = rerouted_sliver_arrangement();
  auto snapper = roofer::reconstruction::createArrangementSnapper();
  snapper->compute(arrangement);

  struct Sample {
    double x, y;
    int segid;
  };
  for (auto [x, y, segid] :
       {Sample{171.009, 216.347, 3}, Sample{168.022, 215.735, 4},
        Sample{174.224, 214.894, 7}, Sample{171.05, 214.362, 6},
        Sample{176.564, 216.405, 5}, Sample{178.941, 219.224, 2},
        Sample{179.474, 216.518, 1}}) {
    CAPTURE(x, y);
    auto face = face_at(arrangement, Point_2(x, y));
    CHECK(face->data().in_footprint);
    CHECK(face->data().segid == segid);
  }
}
