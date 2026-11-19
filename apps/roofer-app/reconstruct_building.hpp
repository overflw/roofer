// Copyright (c) 2018-2025 TU Delft 3D geoinformation group, Ravi Peters (3DGI),
// and Balazs Dukai (3DGI)

// This file is part of roofer (https://github.com/3DBAG/roofer)

// geoflow-roofer was created as part of the 3DBAG project by the TU Delft 3D
// geoinformation group (3d.bk.tudelf.nl) and 3DGI (3dgi.nl)

// geoflow-roofer is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version. geoflow-roofer is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details. You should have received a copy of the GNU
// General Public License along with geoflow-roofer. If not, see
// <https://www.gnu.org/licenses/>.

// Author(s):
// Ravi Peters
#pragma once
#ifdef RF_USE_RERUN
// Adapters so we can log eigen vectors as rerun positions:
template <>
struct rerun::CollectionAdapter<rerun::Position3D, roofer::PointCollection> {
  /// Borrow for non-temporary.
  Collection<rerun::Position3D> operator()(
      const roofer::PointCollection& container) {
    return Collection<rerun::Position3D>::borrow(container.data(),
                                                 container.size());
  }

  // Do a full copy for temporaries (otherwise the data might be deleted when
  // the temporary is destroyed).
  Collection<rerun::Position3D> operator()(
      roofer::PointCollection&& container) {
    std::vector<rerun::Position3D> positions(container.size());
    memcpy(positions.data(), container.data(),
           container.size() * sizeof(roofer::arr3f));
    return Collection<rerun::Position3D>::take_ownership(std::move(positions));
  }
};
template <>
struct rerun::CollectionAdapter<rerun::Position3D, roofer::TriangleCollection> {
  /// Borrow for non-temporary.
  Collection<rerun::Position3D> operator()(
      const roofer::TriangleCollection& container) {
    return Collection<rerun::Position3D>::borrow(container[0].data(),
                                                 container.vertex_count());
  }

  // Do a full copy for temporaries (otherwise the data might be deleted when
  // the temporary is destroyed).
  Collection<rerun::Position3D> operator()(
      roofer::TriangleCollection&& container) {
    std::vector<rerun::Position3D> positions(container.size());
    memcpy(positions.data(), container[0].data(),
           container.vertex_count() * sizeof(roofer::arr3f));
    return Collection<rerun::Position3D>::take_ownership(std::move(positions));
  }
};
#endif

enum LOD { LOD11 = 11, LOD12 = 12, LOD13 = 13, LOD22 = 22 };

void add_ms_to_bbox(roofer::Box& box,
                    std::unordered_map<int, roofer::Mesh>& multisolid) {
  for (auto& [i, mesh] : multisolid) {
    auto& faces = mesh.get_polygons();
    auto& labels = mesh.get_labels();
    assert(faces.size() == labels.size());
    size_t n_faces = faces.size();

    auto& attributes = mesh.get_attributes();
    assert(attributes.size() == n_faces);
    for (size_t i = 0; i < n_faces; ++i) {
      // 1 is the roof label
      if (labels[i] == 1) {
        box.add(faces[i].box());
      }
    }
  }
}

void compute_mesh_properties(
    std::unordered_map<int, roofer::Mesh>& multisolid_lod12,
    std::unordered_map<int, roofer::Mesh>& multisolid_lod13,
    std::unordered_map<int, roofer::Mesh>& multisolid_lod22, float z_offset,
    RooferConfig* cfg) {
  auto MeshPropertyCalculator = roofer::misc::createMeshPropertyCalculator();

  // make sure bbox is initialized properly
  roofer::Box bbox;
  if (multisolid_lod12.size()) {
    add_ms_to_bbox(bbox, multisolid_lod12);
  }
  if (multisolid_lod13.size()) {
    add_ms_to_bbox(bbox, multisolid_lod13);
  }
  if (multisolid_lod22.size()) {
    add_ms_to_bbox(bbox, multisolid_lod22);
  }

  auto heightmap = MeshPropertyCalculator->get_heightmap(
      multisolid_lod22, bbox, cfg->metres_to_input_units(cfg->cellsize));

  for (auto& [i, mesh22] : multisolid_lod22) {
    mesh22.get_attributes().resize(mesh22.get_polygons().size());
    MeshPropertyCalculator->calculate_h_attr(mesh22, heightmap,
                                             {.z_offset = z_offset,
                                              .h_50p = cfg->a_h_roof_50p,
                                              .h_70p = cfg->a_h_roof_70p,
                                              .h_min = cfg->a_h_roof_min,
                                              .h_max = cfg->a_h_roof_max});

    MeshPropertyCalculator->compute_roof_orientation(
        mesh22, {.slope = cfg->a_slope, .azimuth = cfg->a_azimuth});
  }
  for (auto& [i, mesh12] : multisolid_lod12) {
    mesh12.get_attributes().resize(mesh12.get_polygons().size());
    MeshPropertyCalculator->calculate_h_attr(mesh12, heightmap,
                                             {.z_offset = z_offset,
                                              .h_50p = cfg->a_h_roof_50p,
                                              .h_70p = cfg->a_h_roof_70p,
                                              .h_min = cfg->a_h_roof_min,
                                              .h_max = cfg->a_h_roof_max});
  }
  for (auto& [i, mesh13] : multisolid_lod13) {
    mesh13.get_attributes().resize(mesh13.get_polygons().size());
    MeshPropertyCalculator->calculate_h_attr(mesh13, heightmap,
                                             {.z_offset = z_offset,
                                              .h_50p = cfg->a_h_roof_50p,
                                              .h_70p = cfg->a_h_roof_70p,
                                              .h_min = cfg->a_h_roof_min,
                                              .h_max = cfg->a_h_roof_max});
  }
}

void multisolid_post_process(BuildingObject& building, RooferConfig* cfg,
                             LOD lod,
                             std::unordered_map<int, roofer::Mesh>& multisolid,
                             std::optional<float>& rmse,
                             std::optional<float>& volume,
                             std::optional<std::string>& attr_val3dity) {
  auto MeshTriangulator =
      roofer::reconstruction::createMeshTriangulatorLegacy();
  MeshTriangulator->compute(multisolid);
  volume = MeshTriangulator->volumes.at(0);
  // logger.debug("Completed MeshTriangulator");
#ifdef RF_USE_RERUN
  if (cfg->use_rerun) {
    const auto& rec = rerun::RecordingStream::current();
    std::string worldname = fmt::format("world/lod{}/", (int)lod);
    rec.log(worldname + "MeshTriangulator",
            rerun::Mesh3D(MeshTriangulator->triangles)
                .with_vertex_normals(MeshTriangulator->normals)
                .with_class_ids(MeshTriangulator->ring_ids));
  }
#endif

  auto PC2MeshDistCalculator = roofer::misc::createPC2MeshDistCalculator();
  PC2MeshDistCalculator->compute(building.pointcloud_building,
                                 MeshTriangulator->multitrianglecol,
                                 MeshTriangulator->ring_ids);
  rmse = PC2MeshDistCalculator->rms_error;
  // logger.debug("Completed PC2MeshDistCalculator. RMSE={}",
  //  PC2MeshDistCalculator->rms_error);
#ifdef RF_USE_RERUN
// rec.log(worldname+"PC2MeshDistCalculator",
// rerun::Mesh3D(PC2MeshDistCalculator->triangles).with_vertex_normals(MeshTriangulator->normals).with_class_ids(MeshTriangulator->ring_ids));
#endif

#ifdef RF_USE_VAL3DITY
  if (multisolid.size() > 0) {
    auto Val3dator = roofer::misc::createVal3dator();
    Val3dator->compute(multisolid);
    attr_val3dity = Val3dator->errors.front();
  }
  // logger.debug("Completed Val3dator. Errors={}", Val3dator->errors.front());
#endif
}

std::unordered_map<int, roofer::Mesh> extrude_lod22(
    roofer::Arrangement_2 arrangement, BuildingObject& building,
    RooferConfig* cfg,
    roofer::reconstruction::SegmentRasteriserInterface* SegmentRasteriser,
    LOD lod, std::optional<float>& rmse, std::optional<float>& volume,
    std::optional<std::string>& attr_val3dity) {
  bool dissolve_step_edges = false;
  bool dissolve_all_interior = false;
  bool extrude_LoD2 = true;
  if (lod == LOD12) {
    dissolve_all_interior = true, extrude_LoD2 = false;
  } else if (lod == LOD13) {
    dissolve_step_edges = true, extrude_LoD2 = false;
  }
#ifdef RF_USE_RERUN
  const auto& rec = rerun::RecordingStream::current();
  std::string worldname = fmt::format("world/lod{}/", (int)lod);
#endif

  auto& logger = roofer::logger::Logger::get_logger();
  const auto reconstruction = cfg->reconstruction_in_input_units();

  auto elevation_provider =
      roofer::reconstruction::createElevationProvider(*building.h_ground);

  auto ArrangementDissolver =
      roofer::reconstruction::createArrangementDissolver();
  auto dissolver_config = reconstruction.arrangement_dissolver;
  dissolver_config.dissolve_step_edges = dissolve_step_edges;
  dissolver_config.dissolve_all_interior = dissolve_all_interior;
  dissolver_config.step_height_threshold = reconstruction.lod13_step_height;
  ArrangementDissolver->compute(arrangement, SegmentRasteriser->heightfield,
                                *elevation_provider, dissolver_config);
  // logger.debug("Completed ArrangementDissolver");
  // logger.debug("Roof partition has {} faces", arrangement.number_of_faces());
#ifdef RF_USE_RERUN
  if (cfg->use_rerun) {
    rec.log(
        worldname + "ArrangementDissolver",
        rerun::LineStrips3D(roofer::reconstruction::arr2polygons(arrangement)));
  }
#endif
  auto ArrangementSnapper = roofer::reconstruction::createArrangementSnapper();
  ArrangementSnapper->compute(arrangement, *elevation_provider,
                              reconstruction.arrangement_snapper);
  // logger.debug("Completed ArrangementSnapper");
#ifdef RF_USE_RERUN
// rec.log(worldname+"ArrangementSnapper", rerun::LineStrips3D(
// roofer::reconstruction::arr2polygons(arrangement) ));
#endif

  auto ArrangementExtruder =
      roofer::reconstruction::createArrangementExtruder();
  auto extruder_config = reconstruction.arrangement_extruder;
  extruder_config.lod2 = extrude_LoD2;
  ArrangementExtruder->compute(
      arrangement, *elevation_provider,
      lod == LOD22 ? std::span<const roofer::Underpass>(building.underpasses)
                   : std::span<const roofer::Underpass>{},
      extruder_config);
  for (const auto& diagnostic : ArrangementExtruder->diagnostics) {
    const auto& input = building.underpasses.at(diagnostic.input_index);
    const auto id =
        input.source_id.empty()
            ? fmt::format("{}/input {}", building.id, diagnostic.input_index)
            : input.source_id;
    if (diagnostic.code ==
        roofer::reconstruction::UnderpassDiagnosticCode::ceiling_lowered)
      logger.info("Lowering underpass {} ceiling from {} to {}", id,
                  input.ceiling_elevation,
                  *diagnostic.adjusted_ceiling_elevation);
    else
      logger.warning("Skipping underpass {}: {}", id, diagnostic.reason);
  }

#ifdef RF_USE_RERUN
  if (cfg->use_rerun) {
    rec.log(worldname + "ArrangementExtruder",
            rerun::LineStrips3D(ArrangementExtruder->faces)
                .with_class_ids(ArrangementExtruder->labels));
  }
#endif

  multisolid_post_process(building, cfg, lod, ArrangementExtruder->multisolid,
                          rmse, volume, attr_val3dity);

  return ArrangementExtruder->multisolid;
}

void extrude_lod11(BuildingObject& building, float extrusion_h,
                   RooferConfig* cfg) {
  auto SimplePolygonExtruder =
      roofer::reconstruction::createSimplePolygonExtruder();
  SimplePolygonExtruder->compute(building.footprint, *building.h_ground,
                                 extrusion_h);
  // std::vector<std::unordered_map<int, roofer::Mesh>> multisolidvec;
  building.multisolids_lod12 = SimplePolygonExtruder->multisolid;
  building.multisolids_lod13 = SimplePolygonExtruder->multisolid;
  building.multisolids_lod22 = SimplePolygonExtruder->multisolid;

  compute_mesh_properties(building.multisolids_lod12,
                          building.multisolids_lod13,
                          building.multisolids_lod22, building.z_offset, cfg);

  multisolid_post_process(
      building, cfg, LOD11, SimplePolygonExtruder->multisolid,
      building.rmse_lod12, building.volume_lod12, building.val3dity_lod12);
  building.rmse_lod13 = building.rmse_lod12;
  building.rmse_lod22 = building.rmse_lod12;
  building.volume_lod13 = building.volume_lod12;
  building.volume_lod22 = building.volume_lod12;
#ifdef RF_USE_VAL3DITY
  building.val3dity_lod13 = building.val3dity_lod12;
  building.val3dity_lod22 = building.val3dity_lod12;
#endif

  building.extrusion_mode = LOD11_FALLBACK;
  building.roof_elevation_70p = building.h_pc_roof_70p + building.z_offset;
}

void reconstruct_building(BuildingObject& building, RooferConfig* cfg) {
  auto& logger = roofer::logger::Logger::get_logger();
  const auto reconstruction = cfg->reconstruction_in_input_units();

#ifdef RF_USE_RERUN
  // const auto& rec = rerun::RecordingStream::current();
  const auto rec = rerun::RecordingStream(building.output_label());
  // rec.save(building.jsonl_path.string() + ".rrd");
  if (cfg->use_rerun) {
    rec.spawn().exit_on_failure();
    rec.set_thread_local();
  }
  // roofer::vec3f testpts;
  // testpts.push_back({1.0f, 2.0f, 3.0f});
  // rec.log("test", rerun::Points3D(testpts));

#endif

#ifdef RF_USE_RERUN
  // rec.log("world/raw_points",
  //         rerun::AnnotationContext({
  //             rerun::AnnotationInfo(6, "BUILDING", rerun::Rgba32(255, 0, 0)),
  //             rerun::AnnotationInfo(2, "GROUND"),
  //             rerun::AnnotationInfo(1, "UNCLASSIFIED"),
  //         }));
  if (cfg->use_rerun) {
    rec.log("world/building_points",
            rerun::Points3D(building.pointcloud_building));
    rec.log("world/ground_points", rerun::Points3D(building.pointcloud_ground));
  }
#endif

  std::unordered_map<std::string, std::chrono::duration<double>> timings;

  // pointcloud_insufficient is set by StreamCropper when a footprint's
  // building-class point density is below the absolute min_building_density
  // floor.
  if ((building.pointcloud_insufficient && cfg->clear_if_insufficient)) {
    building.extrusion_mode = SKIP;
  }

  if (building.extrusion_mode == SKIP) {
    // If possible still emit a LoD 1.1 block model rather than
    // nothing. This needs a ground elevation for the
    // floor and an explicit height fallback. Without both no geometry
    // can be built, so the building is left empty.
    if (building.h_ground.has_value() && building.roof_h_fallback.has_value()) {
      extrude_lod11(building, *building.roof_h_fallback, cfg);
    }
    return;
  } else if (building.extrusion_mode == LOD11_FALLBACK) {
    extrude_lod11(building, building.h_pc_roof_70p, cfg);
    return;
  } else if (building.extrusion_mode == STANDARD) {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto PlaneDetector = roofer::reconstruction::createPlaneDetector();
    auto PlaneDetector_ground = roofer::reconstruction::createPlaneDetector();
    try {
      auto plane_detector_cfg = reconstruction.plane_detector;
      PlaneDetector->detect(building.pointcloud_building, plane_detector_cfg);
      timings["PlaneDetector"] = std::chrono::high_resolution_clock::now() - t0;
      t0 = std::chrono::high_resolution_clock::now();
      PlaneDetector_ground->detect(building.pointcloud_ground,
                                   plane_detector_cfg);
      timings["PlaneDetector_ground"] =
          std::chrono::high_resolution_clock::now() - t0;

      building.roof_type = PlaneDetector->roof_type;
      building.roof_elevation_50p =
          PlaneDetector->roof_elevation_50p + building.z_offset;
      building.roof_elevation_70p =
          PlaneDetector->roof_elevation_70p + building.z_offset;
      building.roof_elevation_min =
          PlaneDetector->roof_elevation_min + building.z_offset;
      building.roof_elevation_max =
          PlaneDetector->roof_elevation_max + building.z_offset;
      building.roof_n_planes = PlaneDetector->pts_per_roofplane.size();

      bool pointcloud_insufficient = PlaneDetector->roof_type == "no points" ||
                                     PlaneDetector->roof_type == "no planes";
      if (pointcloud_insufficient) {
        building.extrusion_mode = SKIP;
        building.pointcloud_insufficient = true;
        if (building.roof_h_fallback.has_value()) {
          extrude_lod11(building, *building.roof_h_fallback, cfg);
        }
        return;
      }
    } catch (const std::exception& e) {
      // Any failure during plane detection / region growing (the deterministic
      // region-count limit, but also CGAL preconditions, allocation failures,
      // etc.) degrades gracefully to a LoD 1.1 block model rather than failing
      // the building outright.
      extrude_lod11(building, building.h_pc_roof_70p, cfg);
      logger.warning("[reconstructor] {}, LoD1.1 fallback: {}",
                     building.output_label(), e.what());
      return;
    }
    // #ifdef RF_USE_RERUN
    //     rec.log("world/segmented_points",
    //             rerun::Collection{rerun::components::AnnotationContext{
    //                 rerun::datatypes::AnnotationInfo(
    //                     0, "no plane", rerun::datatypes::Rgba32(30, 30,
    //                     30))}});
    //     rec.log(
    //         "world/segmented_points",
    //         rerun::Points3D(points_roof).with_class_ids(PlaneDetector->plane_id));
    // #endif
    // Narrow structures (narrow-planes): a wall's top, a parapet is too
    // narrow for the region grower, its points stay unsegmented or form
    // narrow planes of their own. Those standing above the roof beside them
    // are split into straight runs, each a plane outlined by its rectangle,
    // too narrow for an alpha shape: as triangles in the heightfield and as
    // lines in the arrangement.
    std::map<int, roofer::reconstruction::Rectangle> narrow;
    const double alpha_radius = std::sqrt(reconstruction.alpha_shaper.alpha);
    if (reconstruction.plane_detector.narrow_planes) {
      t0 = std::chrono::high_resolution_clock::now();
      auto& labels = PlaneDetector->plane_id;
      auto& planes = PlaneDetector->pts_per_roofplane;
      const auto& pd = reconstruction.plane_detector;
      const auto& pts = building.pointcloud_building;
      auto outline = [&](int k) {
        std::vector<roofer::reconstruction::Point2> xy;
        for (auto& p : planes[k].second) xy.push_back({p.x(), p.y()});
        return roofer::reconstruction::narrow_outline(xy, alpha_radius);
      };
      auto drop = [&](int k) {
        narrow.erase(k);
        planes.erase(k);
        PlaneDetector->plane_adjacencies.erase(k);
        for (auto& [l, adjacent] : PlaneDetector->plane_adjacencies)
          adjacent.erase(k);
        for (auto& l : labels)
          if (l == k) l = 0;
      };
      std::set<int> steep;
      for (auto& [k, v] : planes)
        if (std::abs(v.first.c()) <
            std::cos(pd.narrow_max_slope * M_PI / 180) *
                std::hypot(v.first.a(), v.first.b(), v.first.c()))
          steep.insert(k);
      auto stands = roofer::reconstruction::standing(
          pts, labels, 2 * alpha_radius, pd.plane_epsilon, steep);
      std::map<int, std::array<size_t, 2>> standing;
      for (size_t i = 0; i < labels.size(); ++i)
        if (labels[i] > 0) {
          auto& s = standing[labels[i]];
          s[0] += stands[i], s[1] += 1;
        }
      auto runs = [&](const std::vector<std::array<double, 3>>& xyz) {
        return roofer::reconstruction::narrow_runs(
            xyz, alpha_radius, 4 * alpha_radius, pd.plane_epsilon);
      };
      auto is_narrow = [&](const std::vector<std::array<double, 3>>& xyz,
                           const std::vector<size_t>& run) {
        std::vector<roofer::reconstruction::Point2> xy;
        for (auto i : run) xy.push_back({xyz[i][0], xyz[i][1]});
        auto o = roofer::reconstruction::narrow_outline(xy, alpha_radius);
        if (!o) return false;
        auto m = roofer::reconstruction::narrow_profile(pts, *o, alpha_radius);
        return m[1] - std::fmax(m[0], m[2]) >= pd.narrow_min_height &&
               !roofer::reconstruction::narrow_on_flank(m, pd.plane_epsilon);
      };
      // the narrow runs of a plane standing on the roof, the tops of walls of
      // one height, go to the loose points; the plane keeps the rest, if enough
      std::vector<bool> released(labels.size(), false);
      const auto labels0 = labels;
      const auto planes0 = planes;
      const auto adjacencies0 = PlaneDetector->plane_adjacencies;
      for (auto& [k, s] : standing) {
        if (2 * s[0] <= s[1]) continue;
        std::vector<size_t> index;
        std::vector<std::array<double, 3>> xyz;
        for (size_t i = 0; i < labels.size(); ++i)
          if (labels[i] == k) {
            index.push_back(i);
            xyz.push_back({pts[i][0], pts[i][1], pts[i][2]});
          }
        for (auto& run : runs(xyz))
          if (is_narrow(xyz, run))
            for (auto j : run) labels[index[j]] = 0, released[index[j]] = true;
        auto& kept = planes[k].second;
        kept.clear();
        for (auto i : index)
          if (labels[i] == k)
            kept.emplace_back(pts[i][0], pts[i][1], pts[i][2]);
        if (kept.size() < size_t(pd.min_plane_points)) drop(k);
      }
      std::vector<size_t> loose;
      std::vector<std::array<double, 3>> xyz;
      for (size_t i = 0; i < labels.size(); ++i)
        if (labels[i] == 0 && (stands[i] || released[i])) {
          loose.push_back(i);
          xyz.push_back({pts[i][0], pts[i][1], pts[i][2]});
        }
      int id = 0;
      for (auto& [k, v] : planes) id = std::max(id, k);
      for (auto& run : runs(xyz)) {
        if (!is_narrow(xyz, run)) continue;
        ++id;
        for (auto i : run) {
          planes[id].second.emplace_back(xyz[i][0], xyz[i][1], xyz[i][2]);
          labels[loose[i]] = id;
        }
        narrow[id] = *outline(id);
      }
      auto level_across = [&](int k) {
        std::vector<std::array<double, 3>> xyz;
        for (auto& p : planes[k].second) xyz.push_back({p.x(), p.y(), p.z()});
        auto c = roofer::reconstruction::narrow_plane(xyz, narrow[k]);
        planes[k].first = roofer::Plane(c[0], c[1], c[2], c[3]);
      };
      for (auto& [k, r] : narrow) level_across(k);
      // pieces of one narrow structure at one height become one plane
      auto height = [&](int k, double x, double y) {
        const auto& pl = planes[k].first;
        return -(pl.a() * x + pl.b() * y + pl.d()) / pl.c();
      };
      for (bool merged = true; merged;) {
        merged = false;
        for (auto a = narrow.begin(); a != narrow.end() && !merged;) {
          auto next = std::next(a);
          for (auto b = std::next(a); b != narrow.end(); ++b) {
            double x = (a->second[0][0] + b->second[2][0]) / 2,
                   y = (a->second[0][1] + b->second[2][1]) / 2;
            if (!roofer::reconstruction::narrow_in_line(a->second, b->second,
                                                        alpha_radius) ||
                std::abs(height(a->first, x, y) - height(b->first, x, y)) >
                    pd.plane_epsilon)
              continue;
            std::vector<roofer::reconstruction::Point2> xy;
            for (auto k : {a->first, b->first})
              for (auto& p : planes[k].second) xy.push_back({p.x(), p.y()});
            auto r = roofer::reconstruction::narrow_outline(xy, alpha_radius);
            if (!r) continue;
            auto& pb = planes[b->first].second;
            planes[a->first].second.insert(planes[a->first].second.end(),
                                           pb.begin(), pb.end());
            for (auto& l : labels)
              if (l == b->first) l = a->first;
            drop(b->first);
            a->second = *r;
            level_across(a->first);
            merged = true;
            break;
          }
          if (!merged) a = next;
        }
      }
      std::vector<std::vector<roofer::reconstruction::Point2>> rings(1);
      for (auto& p : building.footprint) rings[0].push_back({p[0], p[1]});
      for (auto& hole : building.footprint.interior_rings()) {
        rings.emplace_back();
        for (auto& p : hole) rings.back().push_back({p[0], p[1]});
      }
      std::vector<roofer::arr3f> unsegmented;
      for (size_t i = 0; i < labels.size(); ++i)
        if (labels[i] == 0) unsegmented.push_back(pts[i]);
      std::vector<int> short_or_on_border;
      for (auto& [k, r] : narrow) {
        r = roofer::reconstruction::narrow_extend(
            r, unsegmented, [&](double x, double y) { return height(k, x, y); },
            pd.plane_epsilon, 2 * alpha_radius);
        r = roofer::reconstruction::narrow_to_border(r, rings,
                                                     4 * alpha_radius);
        auto mid = [&](int i, int j) {
          return roofer::reconstruction::Point2{(r[i][0] + r[j][0]) / 2,
                                                (r[i][1] + r[j][1]) / 2};
        };
        auto inside = [&](roofer::reconstruction::Point2 p) {
          return roofer::reconstruction::to_border(rings, p[0], p[1]) >
                 2 * alpha_radius;
        };
        if (std::hypot(r[1][0] - r[0][0], r[1][1] - r[0][1]) <
                4 * alpha_radius ||
            !(inside(mid(0, 2)) || inside(mid(0, 3)) || inside(mid(1, 2))))
          short_or_on_border.push_back(k);
      }
      for (auto k : short_or_on_border) drop(k);
      if (narrow.empty()) {
        labels = labels0;
        planes = planes0;
        PlaneDetector->plane_adjacencies = adjacencies0;
      }
      timings["NarrowPlanes"] = std::chrono::high_resolution_clock::now() - t0;
    }

    t0 = std::chrono::high_resolution_clock::now();
    auto AlphaShaper = roofer::reconstruction::createAlphaShaper();
    AlphaShaper->compute(PlaneDetector->pts_per_roofplane,
                         reconstruction.alpha_shaper);
    timings["AlphaShaper"] = std::chrono::high_resolution_clock::now() - t0;
    if (!narrow.empty()) {
      std::vector<roofer::LinearRing> alpha_rings;
      roofer::vec1i ids;
      for (size_t i = 0; i < AlphaShaper->alpha_rings.size(); ++i) {
        if (narrow.count(AlphaShaper->roofplane_ids[i])) continue;
        alpha_rings.push_back(AlphaShaper->alpha_rings[i]);
        ids.push_back(AlphaShaper->roofplane_ids[i]);
      }
      AlphaShaper->alpha_rings = std::move(alpha_rings);
      AlphaShaper->roofplane_ids = std::move(ids);
      for (const auto& [id, r] : narrow) {
        const auto& pl = PlaneDetector->pts_per_roofplane.at(id).first;
        auto z = [&](const roofer::reconstruction::Point2& p) {
          return roofer::arr3f{
              float(p[0]), float(p[1]),
              float(-(pl.a() * p[0] + pl.b() * p[1] + pl.d()) / pl.c())};
        };
        AlphaShaper->alpha_triangles.push_back({z(r[0]), z(r[1]), z(r[2])});
        AlphaShaper->alpha_triangles.push_back({z(r[0]), z(r[2]), z(r[3])});
      }
    }
    // logger.debug("Completed AlphaShaper (roof), found {} rings, {} labels",
    //  AlphaShaper->alpha_rings.size(),
    //  AlphaShaper->roofplane_ids.size());
#ifdef RF_USE_RERUN
    if (cfg->use_rerun) {
      rec.log("world/alpha_rings_roof",
              rerun::LineStrips3D(AlphaShaper->alpha_rings)
                  .with_class_ids(AlphaShaper->roofplane_ids));
    }
#endif
    t0 = std::chrono::high_resolution_clock::now();
    auto AlphaShaper_ground = roofer::reconstruction::createAlphaShaper();
    AlphaShaper_ground->compute(PlaneDetector_ground->pts_per_roofplane,
                                reconstruction.alpha_shaper);
    timings["AlphaShaper_ground"] =
        std::chrono::high_resolution_clock::now() - t0;
    // logger.debug("Completed AlphaShaper (ground), found {} rings, {} labels",
    //  AlphaShaper_ground->alpha_rings.size(),
    //  AlphaShaper_ground->roofplane_ids.size());
#ifdef RF_USE_RERUN
    if (cfg->use_rerun) {
      rec.log("world/alpha_rings_ground",
              rerun::LineStrips3D(AlphaShaper_ground->alpha_rings)
                  .with_class_ids(AlphaShaper_ground->roofplane_ids));
    }
#endif
    t0 = std::chrono::high_resolution_clock::now();
    auto LineDetector = roofer::reconstruction::createLineDetector();
    LineDetector->detect(AlphaShaper->alpha_rings, AlphaShaper->roofplane_ids,
                         PlaneDetector->pts_per_roofplane,
                         reconstruction.line_detector);
    timings["LineDetector"] = std::chrono::high_resolution_clock::now() - t0;
    // logger.debug("Completed LineDetector");
#ifdef RF_USE_RERUN
    if (cfg->use_rerun) {
      rec.log("world/boundary_lines",
              rerun::LineStrips3D(LineDetector->edge_segments));
    }
#endif

    t0 = std::chrono::high_resolution_clock::now();
    auto PlaneIntersector = roofer::reconstruction::createPlaneIntersector();
    PlaneIntersector->compute(PlaneDetector->pts_per_roofplane,
                              PlaneDetector->plane_adjacencies,
                              reconstruction.plane_intersector);
    timings["PlaneIntersector"] =
        std::chrono::high_resolution_clock::now() - t0;

    size_t hr_i;
    float hr_z;
    building.roof_n_ridgelines =
        PlaneIntersector->find_highest_ridgeline(hr_z, hr_i);
    if (building.roof_n_ridgelines > 0) {
      building.roof_elevation_ridge = hr_z;
    }
    // logger.debug("Completed PlaneIntersector");
#ifdef RF_USE_RERUN
    if (cfg->use_rerun) {
      rec.log("world/intersection_lines",
              rerun::LineStrips3D(PlaneIntersector->segments));
    }
#endif

    t0 = std::chrono::high_resolution_clock::now();
    auto LineRegulariser = roofer::reconstruction::createLineRegulariser();
    LineRegulariser->compute(
        LineDetector->edge_segments, PlaneIntersector->segments,
        reconstruction.line_regulariser, building.external_lines);
    timings["LineRegulariser"] = std::chrono::high_resolution_clock::now() - t0;
    // past the regulariser, which would merge a narrow outline's sides
    for (const auto& [id, r] : narrow) {
      for (auto& [p, q] : roofer::reconstruction::narrow_edges(
               r, reconstruction.line_regulariser.extension)) {
        LineRegulariser->exact_regularised_edges.emplace_back(
            roofer::EPECK::Point_2(p[0], p[1]),
            roofer::EPECK::Point_2(q[0], q[1]));
      }
    }
    // logger.debug("Completed LineRegulariser");
#ifdef RF_USE_RERUN
    if (cfg->use_rerun) {
      rec.log("world/regularised_lines",
              rerun::LineStrips3D(LineRegulariser->regularised_edges));
    }
#endif

    t0 = std::chrono::high_resolution_clock::now();
    auto SegmentRasteriser = roofer::reconstruction::createSegmentRasteriser();
    auto rasteriser_config = reconstruction.segment_rasteriser;
    rasteriser_config.use_ground =
        !building.pointcloud_ground.empty() && reconstruction.clip_terrain;
    SegmentRasteriser->compute(AlphaShaper->alpha_triangles,
                               AlphaShaper_ground->alpha_triangles,
                               rasteriser_config);
    timings["SegmentRasteriser"] =
        std::chrono::high_resolution_clock::now() - t0;
    // logger.debug("Completed SegmentRasteriser");

#ifdef RF_USE_RERUN
    if (cfg->use_rerun) {
      auto heightfield_copy = SegmentRasteriser->heightfield;
      heightfield_copy.set_nodata(0);
      rec.log(
          "world/heightfield",
          rerun::DepthImage(heightfield_copy.vals_->data(),
                            {static_cast<uint32_t>(heightfield_copy.dimx_),
                             static_cast<uint32_t>(heightfield_copy.dimy_)}));
    }
#endif

    t0 = std::chrono::high_resolution_clock::now();
    roofer::Arrangement_2 arrangement;
    auto ArrangementBuilder =
        roofer::reconstruction::createArrangementBuilder();
    ArrangementBuilder->compute(arrangement, building.footprint,
                                LineRegulariser->exact_regularised_edges,
                                reconstruction.arrangement_builder);
    timings["ArrangementBuilder"] =
        std::chrono::high_resolution_clock::now() - t0;
    // logger.debug("Completed ArrangementBuilder");
    // logger.debug("Roof partition has {} faces",
    // arrangement.number_of_faces());
#ifdef RF_USE_RERUN
    if (cfg->use_rerun) {
      rec.log("world/initial_partition",
              rerun::LineStrips3D(
                  roofer::reconstruction::arr2polygons(arrangement)));
    }
#endif

    t0 = std::chrono::high_resolution_clock::now();
    auto ArrangementOptimiser =
        roofer::reconstruction::createArrangementOptimiser();
    auto optimiser_config = reconstruction.arrangement_optimiser;
    optimiser_config.use_ground =
        !building.pointcloud_ground.empty() && reconstruction.clip_terrain;
    ArrangementOptimiser->compute(arrangement, SegmentRasteriser->heightfield,
                                  PlaneDetector->pts_per_roofplane,
                                  PlaneDetector_ground->pts_per_roofplane,
                                  optimiser_config);
    timings["ArrangementOptimiser"] =
        std::chrono::high_resolution_clock::now() - t0;
    // A narrow plane can fit a cell off its outline as well as the roof
    // does; there the cell takes the roof plane that fits it best.
    if (!narrow.empty()) {
      std::vector<std::pair<int, roofer::Plane>>
          roof_labels;  // the optimiser's
      for (auto& [k, v] : PlaneDetector->pts_per_roofplane)
        if (k >= 1) roof_labels.emplace_back(k, v.first);
      for (auto f : arrangement.face_handles()) {
        auto n = narrow.find(f->data().segid);
        if (!f->data().in_footprint || n == narrow.end()) continue;
        double cx = 0, cy = 0;
        int count = 0;
        auto he = f->outer_ccb(), first = he;
        do {
          cx += CGAL::to_double(he->source()->point().x());
          cy += CGAL::to_double(he->source()->point().y());
          ++count;
        } while (++he != first);
        if (roofer::reconstruction::in_rectangle(n->second, cx / count,
                                                 cy / count, 0))
          continue;
        const auto& cost = f->data().vertex_label_cost;
        size_t best = cost.size();
        for (size_t i = 0; i < roof_labels.size() && i < cost.size(); ++i)
          if (!narrow.count(roof_labels[i].first) &&
              (best == cost.size() || cost[i] < cost[best]))
            best = i;
        if (best == cost.size()) continue;
        f->data().label = best;
        f->data().segid = roof_labels[best].first;
        f->data().plane = roof_labels[best].second;
      }
    }
    // logger.debug("Completed ArrangementOptimiser");
    // rec.log("world/optimised_partition", rerun::LineStrips3D(
    // roofer::reconstruction::arr2polygons(arrangement) ));

    // LoDs
    // attributes to be filled during reconstruction
    // logger.debug("LoD={}", cfg->lod);
    t0 = std::chrono::high_resolution_clock::now();
    if (cfg->reconstruction.lod12) {
      building.multisolids_lod12 = extrude_lod22(
          arrangement, building, cfg, SegmentRasteriser.get(), LOD12,
          building.rmse_lod12, building.volume_lod12, building.val3dity_lod12);
    }

    if (cfg->reconstruction.lod13) {
      building.multisolids_lod13 = extrude_lod22(
          arrangement, building, cfg, SegmentRasteriser.get(), LOD13,
          building.rmse_lod13, building.volume_lod13, building.val3dity_lod13);
    }

    if (cfg->reconstruction.lod22) {
      building.multisolids_lod22 = extrude_lod22(
          arrangement, building, cfg, SegmentRasteriser.get(), LOD22,
          building.rmse_lod22, building.volume_lod22, building.val3dity_lod22);
      compute_mesh_properties(
          building.multisolids_lod12, building.multisolids_lod13,
          building.multisolids_lod22, building.z_offset, cfg);
      // a roof polygon in a narrow outline and on its plane is that structure
      if (!cfg->a_narrow.empty() && !narrow.empty()) {
        for (auto& [i, mesh] : building.multisolids_lod22) {
          auto& polygons = mesh.get_polygons();
          auto& labels = mesh.get_labels();
          auto& attributes = mesh.get_attributes();
          for (size_t j = 0; j < polygons.size(); ++j) {
            if (labels[j] != 1 || polygons[j].empty()) continue;
            std::array<double, 3> c{0, 0, 0};
            for (auto& p : polygons[j])
              for (int k = 0; k < 3; ++k) c[k] += p[k] / polygons[j].size();
            for (const auto& [id, r] : narrow) {
              const auto& pl = PlaneDetector->pts_per_roofplane.at(id).first;
              if (roofer::reconstruction::in_rectangle(r, c[0], c[1], 0.05) &&
                  std::abs(c[2] + (pl.a() * c[0] + pl.b() * c[1] + pl.d()) /
                                      pl.c()) < 0.1) {
                attributes[j].insert(cfg->a_narrow, true);
                break;
              }
            }
          }
        }
      }
    }
    timings["extrude"] = std::chrono::high_resolution_clock::now() - t0;

    std::string timings_str =
        fmt::format("[reconstructor t] {} (", building.output_label());
    for (const auto& [key, value] : timings) {
      auto ms = static_cast<int>(
          std::chrono::duration_cast<std::chrono::milliseconds>(value).count());
      timings_str += fmt::format("({}, {}),", key, ms);
    }
    logger.debug("{})", timings_str);
  }
}
