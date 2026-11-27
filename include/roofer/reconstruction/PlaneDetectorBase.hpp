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

#include <algorithm>
#include <cmath>
#include <deque>
#include <map>
#include <numeric>
#include <queue>
#include <tuple>
#include <vector>
#include <roofer/common/datastructures.hpp>
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

    struct PlaneFit {
      double off_share, rms;
    };

    inline PlaneFit plane_fit(PlaneDS& cds, const std::vector<size_t>& idx,
                              const Plane& plane, double off) {
      size_t n_off = 0;
      double sum = 0;
      for (auto i : idx) {
        double d2 = CGAL::squared_distance(
            plane, Point(cds.points[i][0], cds.points[i][1], cds.points[i][2]));
        n_off += d2 > off * off;
        sum += d2;
      }
      return {double(n_off) / idx.size(), std::sqrt(sum / idx.size())};
    }

    // Merge adjacent regions that one plane fits about as well as their own;
    // max_angle in degrees.
    inline size_t merge_coplanar_regions(
        PlaneDS& cds, std::vector<PlaneRegion>& regions,
        const std::map<size_t, std::map<size_t, size_t>>& adjacencies,
        double max_angle, double off, double rms_slack) {
      std::map<size_t, size_t> index;
      std::vector<PlaneFit> own(regions.size());
      for (size_t i = 0; i < regions.size(); ++i) {
        if (regions[i].get_region_id() == 0) continue;
        index[regions[i].get_region_id()] = i;
        Plane plane;
        cds.fit_plane(regions[i].inliers, plane);
        own[i] = plane_fit(cds, regions[i].inliers, plane, off);
      }
      std::vector<std::tuple<double, size_t, size_t>> pairs;
      for (const auto& [a, adjacent] : adjacencies) {
        for (const auto& [b, count] : adjacent) {
          auto ia = index.find(a), ib = index.find(b);
          if (ia == index.end() || ib == index.end()) continue;
          size_t i = ia->second, j = ib->second;
          double angle =
              CGAL::approximate_angle(regions[i].plane.orthogonal_vector(),
                                      regions[j].plane.orthogonal_vector());
          angle = std::min(angle, 180.0 - angle);
          if (angle < max_angle) pairs.emplace_back(angle, i, j);
        }
      }
      std::sort(pairs.begin(), pairs.end());

      std::vector<size_t> group_of(regions.size());
      std::iota(group_of.begin(), group_of.end(), 0);
      std::vector<std::vector<size_t>> members(regions.size());
      for (size_t i = 0; i < regions.size(); ++i) members[i] = {i};
      size_t n_merges = 0;
      for (auto [angle, i, j] : pairs) {
        auto [ri, rj] = std::minmax({group_of[i], group_of[j]});
        if (ri == rj) continue;
        std::vector<size_t> points;
        for (auto m : {ri, rj})
          for (auto q : members[m])
            points.insert(points.end(), regions[q].inliers.begin(),
                          regions[q].inliers.end());
        Plane plane;
        cds.fit_plane(points, plane);
        auto fits = [&](size_t q) {
          auto fit = plane_fit(cds, regions[q].inliers, plane, off);
          return fit.off_share <= own[q].off_share + 0.01 &&
                 fit.rms <= own[q].rms + rms_slack;
        };
        if (!std::ranges::all_of(members[ri], fits) ||
            !std::ranges::all_of(members[rj], fits))
          continue;
        for (auto q : members[rj]) group_of[q] = ri;
        members[ri].insert(members[ri].end(), members[rj].begin(),
                           members[rj].end());
        regions[ri].plane = plane;
        ++n_merges;
      }

      for (size_t r = 0; r < regions.size(); ++r) {
        if (group_of[r] != r) continue;
        for (auto q : members[r]) {
          if (q == r) continue;
          regions[r].inliers.insert(regions[r].inliers.end(),
                                    regions[q].inliers.begin(),
                                    regions[q].inliers.end());
          regions[q].inliers.clear();
        }
      }
      regions.erase(std::remove_if(regions.begin(), regions.end(),
                                   [](auto& r) {
                                     return r.get_region_id() != 0 &&
                                            r.inliers.empty();
                                   }),
                    regions.end());
      return n_merges;
    }

  }  // namespace planedect

}  // namespace roofer
