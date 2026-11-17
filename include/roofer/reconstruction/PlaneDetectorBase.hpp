// Copyright (c) 2018-2026 TU Delft 3D geoinformation group, Ravi Peters (3DGI),
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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Plane_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <deque>
#include <queue>
#include <algorithm>
#include <cmath>
#include <map>
#include <roofer/common/datastructures.hpp>
#include <set>
#include <utility>
#include <vector>
#include <roofer/reconstruction/RegionGrower.hpp>
#include <roofer/reconstruction/RegionGrower_DS_CGAL.hpp>
#include <roofer/reconstruction/cgal_shared_definitions.hpp>

namespace roofer {

  namespace planedect {

    class PlaneDS : public regiongrower::CGAL_RegionGrowerDS {
     public:
      roofer::vec3f& normals;
      std::vector<Plane> seed_planes;

      PlaneDS(roofer::PointCollection& points, roofer::vec3f& normals,
              size_t N = 15)
          : CGAL_RegionGrowerDS(points, N), normals(normals){};

      // Note this crashes when idx.size()==1;
      inline double fit_plane(std::vector<size_t>& idx, Plane& plane) {
        std::vector<Point> neighbor_points;
        for (auto i : idx)
          neighbor_points.push_back(
              Point(points[i][0], points[i][1], points[i][2]));
        double quality = linear_least_squares_fitting_3(
            neighbor_points.begin(), neighbor_points.end(), plane,
            CGAL::Dimension_tag<0>());
        return quality;
      }

      virtual std::deque<size_t> get_seeds() override {
        // seed generation
        typedef std::pair<size_t, double> index_dist_pair;
        auto cmp = [](index_dist_pair left, index_dist_pair right) {
          return left.second < right.second;
        };
        std::priority_queue<index_dist_pair, std::vector<index_dist_pair>,
                            decltype(cmp)>
            pq(cmp);

        size_t i = 0;
        Plane plane;
        for (size_t pi = 0; pi < size; ++pi) {
          auto neighbours = get_neighbours(pi);
          auto quality = fit_plane(neighbours, plane);
          pq.push(index_dist_pair(i++, quality));
        }

        std::deque<size_t> seed_order;
        while (pq.size() > 0) {
          seed_order.push_back(pq.top().first);
          pq.pop();
        }
        return seed_order;
      }
    };

    class PlaneRegion : public regiongrower::Region {
     public:
      using regiongrower::Region::Region;
      Plane plane;
      std::vector<size_t> inliers;
    };

    class DistAndNormalTester {
     public:
      float dist_thres;
      float normal_thres;
      size_t n_refit, refit_counter = 0;

      DistAndNormalTester(float dist_thres = 0.04, float normal_thres = 0.9,
                          size_t n_refit = 5)
          : dist_thres(dist_thres),
            normal_thres(normal_thres),
            n_refit(n_refit){};

      bool is_valid(PlaneDS& cds, size_t candidate, size_t neighbour,
                    PlaneRegion& shape) {
        Point p_c = Point(cds.points[candidate][0], cds.points[candidate][1],
                          cds.points[candidate][2]);
        Vector n_c =
            Vector(cds.normals[candidate][0], cds.normals[candidate][1],
                   cds.normals[candidate][2]);
        Point p = Point(cds.points[neighbour][0], cds.points[neighbour][1],
                        cds.points[neighbour][2]);
        Vector n = Vector(cds.normals[neighbour][0], cds.normals[neighbour][1],
                          cds.normals[neighbour][2]);

        if (shape.inliers.size() == 0) {
          shape.plane = Plane(p_c, n_c);
          shape.inliers.push_back(candidate);
        }

        bool valid =
            (CGAL::squared_distance(shape.plane, p) < dist_thres) &&
            (std::abs(shape.plane.orthogonal_vector() * n) > normal_thres);
        if (valid) {
          shape.inliers.push_back(neighbour);
          if (shape.inliers.size() % n_refit == 0) {
            cds.fit_plane(shape.inliers, shape.plane);
          }
        }
        return valid;
      }
    };

    struct RegionMergeParams {
      // maximum angle between the normals of two regions, in radians
      double max_angle;
      // a merge is accepted when, for both regions, the RMS distance of their
      // points to the merged plane is at most
      //   rms_factor * max(rms_a, rms_b) + rms_offset
      double rms_factor;
      double rms_offset;
    };

    inline double plane_rms(PlaneDS& cds, const std::vector<size_t>& idx,
                            const Plane& plane) {
      if (idx.empty()) return 0;
      double sum = 0;
      for (auto i : idx) {
        sum += CGAL::squared_distance(
            plane, Point(cds.points[i][0], cds.points[i][1], cds.points[i][2]));
      }
      return std::sqrt(sum / idx.size());
    }

    /**
     * Merge adjacent regions that are one plane.
     *
     * Region growing can split one roof plane into several regions, for
     * example when the plane drifts while growing or when noisy normals stop
     * the growth. This refits every region on all its inliers and then
     * greedily merges the adjacent pair whose merged plane fits best, as long
     * as the merged plane fits the points of both regions about as well as
     * their own planes do. Merged regions keep the lowest region id.
     *
     * @param regions regions from the region grower; index 0 is the
     * unsegmented placeholder and is left untouched
     * @param adjacencies adjacency counts keyed by region id (any order)
     * @param region_ids region id per point, updated for merged points
     * @return number of merges
     */
    template <typename RegionVec>
    size_t merge_coplanar_regions(
        PlaneDS& cds, RegionVec& regions,
        const std::map<size_t, std::map<size_t, size_t>>& adjacencies,
        std::vector<size_t>& region_ids, const RegionMergeParams& params) {
      std::map<size_t, size_t> id_to_index;
      std::map<size_t, double> rms;
      for (size_t i = 0; i < regions.size(); ++i) {
        auto& region = regions[i];
        if (region.get_region_id() == 0 || region.inliers.size() < 3) continue;
        id_to_index[region.get_region_id()] = i;
        cds.fit_plane(region.inliers, region.plane);
        rms[region.get_region_id()] =
            plane_rms(cds, region.inliers, region.plane);
      }

      std::map<size_t, std::set<size_t>> neighbours;
      for (const auto& [a, adj] : adjacencies) {
        if (!id_to_index.count(a)) continue;
        for (const auto& [b, count] : adj) {
          if (a == b || !id_to_index.count(b)) continue;
          neighbours[a].insert(b);
          neighbours[b].insert(a);
        }
      }

      const double min_cos = std::cos(params.max_angle);
      struct Candidate {
        bool valid = false;
        double cost = 0;
        Plane plane;
      };
      std::map<std::pair<size_t, size_t>, Candidate> cache;
      auto evaluate = [&](size_t a, size_t b) -> Candidate {
        auto key = std::make_pair(std::min(a, b), std::max(a, b));
        auto it = cache.find(key);
        if (it != cache.end()) return it->second;
        Candidate c;
        auto& ra = regions[id_to_index[a]];
        auto& rb = regions[id_to_index[b]];
        auto na = ra.plane.orthogonal_vector();
        auto nb = rb.plane.orthogonal_vector();
        double cos_angle = std::abs(CGAL::to_double(na * nb)) /
                           std::sqrt(CGAL::to_double(na.squared_length() *
                                                     nb.squared_length()));
        if (cos_angle >= min_cos) {
          std::vector<size_t> both(ra.inliers);
          both.insert(both.end(), rb.inliers.begin(), rb.inliers.end());
          cds.fit_plane(both, c.plane);
          double limit =
              params.rms_factor * std::max(rms[a], rms[b]) + params.rms_offset;
          double rms_a = plane_rms(cds, ra.inliers, c.plane);
          double rms_b = plane_rms(cds, rb.inliers, c.plane);
          if (rms_a <= limit && rms_b <= limit) {
            c.valid = true;
            c.cost = std::max(rms_a - rms[a], rms_b - rms[b]);
          }
        }
        cache[key] = c;
        return c;
      };

      size_t n_merges = 0;
      while (true) {
        bool found = false;
        size_t best_a = 0, best_b = 0;
        Candidate best;
        for (const auto& [a, nbs] : neighbours) {
          for (auto b : nbs) {
            if (b <= a) continue;
            auto c = evaluate(a, b);
            if (c.valid && (!found || c.cost < best.cost)) {
              found = true;
              best = c;
              best_a = a;
              best_b = b;
            }
          }
        }
        if (!found) break;

        // merge best_b into best_a
        auto& ra = regions[id_to_index[best_a]];
        auto& rb = regions[id_to_index[best_b]];
        for (auto i : rb.inliers) region_ids[i] = best_a;
        ra.inliers.insert(ra.inliers.end(), rb.inliers.begin(),
                          rb.inliers.end());
        rb.inliers.clear();
        ra.plane = best.plane;
        rms[best_a] = plane_rms(cds, ra.inliers, ra.plane);
        ++n_merges;

        for (auto c : neighbours[best_b]) {
          neighbours[c].erase(best_b);
          if (c != best_a) {
            neighbours[c].insert(best_a);
            neighbours[best_a].insert(c);
          }
        }
        neighbours.erase(best_b);
        neighbours[best_a].erase(best_a);
        for (auto it = cache.begin(); it != cache.end();) {
          if (it->first.first == best_a || it->first.second == best_a ||
              it->first.first == best_b || it->first.second == best_b) {
            it = cache.erase(it);
          } else {
            ++it;
          }
        }
      }

      // drop the emptied regions
      regions.erase(std::remove_if(regions.begin() + 1, regions.end(),
                                   [](auto& r) { return r.inliers.empty(); }),
                    regions.end());
      return n_merges;
    }

  }  // namespace planedect

}  // namespace roofer
