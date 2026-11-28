// Copyright (c) 2018-2026 TU Delft 3D geoinformation group, Ravi Peters (3DGI),
// and Balazs Dukai (3DGI)

// This file is part of roofer (https://github.com/3DBAG/roofer)

// geoflow-roofer is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version. geoflow-roofer is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details. You should have received a copy of the GNU
// General Public License along with geoflow-roofer. If not, see
// <https://www.gnu.org/licenses/>.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <map>
#include <optional>
#include <set>
#include <vector>

namespace roofer::reconstruction {

  using Point2 = std::array<double, 2>;
  using Rectangle = std::array<Point2, 4>;

  // The unit vector along `r` (corners 0 to 1) and its length.
  inline std::array<double, 3> rectangle_axis(const Rectangle& r) {
    double ux = r[1][0] - r[0][0], uy = r[1][1] - r[0][1],
           l = std::hypot(ux, uy);
    return {ux / l, uy / l, l};
  }

  // The midpoint of corners `i` and `j` of `r`.
  inline Point2 mid(const Rectangle& r, int i, int j) {
    return {(r[i][0] + r[j][0]) / 2, (r[i][1] + r[j][1]) / 2};
  }

  // The height of `plane` (a, b, c, d) at (x, y).
  template <typename Plane>
  double plane_z(const Plane& plane, double x, double y) {
    return -(plane.a() * x + plane.b() * y + plane.d()) / plane.c();
  }

  /**
   * The outline of a narrow plane: a wall's top, a parapet, a pipe. An alpha
   * shape of `alpha_radius` breaks a region narrower than about twice that
   * radius into fragments, so a plane whose points spread across it by less
   * than half the radius (standard deviation) and run on at least four times
   * as long takes the rectangle around them along their main direction; the
   * few that lie off their line, as at a corner, give it no width.
   *
   * @param xy the plane's points in plan
   * @return the rectangle, corners in order, or nothing for a plane that is
   * not narrow
   */
  inline std::optional<Rectangle> narrow_outline(const std::vector<Point2>& xy,
                                                 double alpha_radius) {
    // points off its line, as at a corner with another wall, give no width
    auto line = xy;
    for (bool trimmed = false;; trimmed = true) {
      if (line.size() < 3) return std::nullopt;
      double mx = 0, my = 0;
      for (auto& p : line) mx += p[0], my += p[1];
      mx /= line.size(), my /= line.size();
      double sxx = 0, syy = 0, sxy = 0;
      for (auto& p : line) {
        double dx = p[0] - mx, dy = p[1] - my;
        sxx += dx * dx, syy += dy * dy, sxy += dx * dy;
      }
      sxx /= line.size(), syy /= line.size(), sxy /= line.size();
      double t = std::atan2(2 * sxy, sxx - syy) / 2, ux = std::cos(t),
             uy = std::sin(t);
      auto v_of = [&](const Point2& p) {
        return -(p[0] - mx) * uy + (p[1] - my) * ux;
      };
      double vmin = INFINITY, vmax = -INFINITY, vv = 0;
      for (auto& p : line) {
        double v = v_of(p);
        vmin = std::min(vmin, v), vmax = std::max(vmax, v), vv += v * v;
      }
      if (!trimmed) {
        double spread = std::sqrt(vv / line.size());
        if (spread >= alpha_radius / 2) return std::nullopt;
        std::erase_if(line, [&](const Point2& p) {
          return std::abs(v_of(p)) > 2 * spread;
        });
        continue;
      }
      double umin = INFINITY, umax = -INFINITY;
      for (auto& p : xy) {
        double u = (p[0] - mx) * ux + (p[1] - my) * uy;
        umin = std::min(umin, u), umax = std::max(umax, u);
      }
      if (umax - umin < 2 * alpha_radius || umax - umin < 4 * (vmax - vmin))
        return std::nullopt;
      // points on one line (a scan line) still make a rectangle, not a segment
      double pad = std::max(0.0, alpha_radius / 10 - (vmax - vmin)) / 2;
      vmin -= pad, vmax += pad;
      auto at = [&](double u, double v) {
        return Point2{mx + u * ux - v * uy, my + u * uy + v * ux};
      };
      return Rectangle{at(umin, vmin), at(umax, vmin), at(umax, vmax),
                       at(umin, vmax)};
    }
  }

  /**
   * The lines a narrow plane's outline puts in the arrangement: its long
   * sides, carried on like the regularised lines by `extension`, but no
   * farther than the outline is long, and its ends.
   */
  inline std::vector<std::array<Point2, 2>> narrow_edges(const Rectangle& r,
                                                         double extension) {
    auto [ux, uy, l] = rectangle_axis(r);
    extension = std::min(extension, l);
    double ex = ux * extension, ey = uy * extension;
    return {{Point2{r[0][0] - ex, r[0][1] - ey},
             Point2{r[1][0] + ex, r[1][1] + ey}},
            {Point2{r[3][0] - ex, r[3][1] - ey},
             Point2{r[2][0] + ex, r[2][1] + ey}},
            {r[0], r[3]},
            {r[1], r[2]}};
  }

  /**
   * The plane (a, b, c, d) of a narrow structure's points `xyz`: its width
   * leaves the slope across it undetermined, so it is fitted along the
   * outline `r` and level across.
   */
  template <typename Points>
  std::array<double, 4> narrow_plane(const Points& xyz, const Rectangle& r) {
    auto [ux, uy, l] = rectangle_axis(r);
    double n = 0, su = 0, sz = 0, suu = 0, suz = 0;
    for (auto& p : xyz) {
      double u = p[0] * ux + p[1] * uy;
      n += 1, su += u, sz += p[2], suu += u * u, suz += u * p[2];
    }
    double den = n * suu - su * su,
           g = std::abs(den) > 1e-9 ? (n * suz - su * sz) / den : 0,
           z0 = (sz - g * su) / n;
    return {-g * ux, -g * uy, 1, -z0};  // z = z0 + g * u
  }

  /**
   * Whether two narrow outlines are pieces of one structure, broken where its
   * points thin out: within 10 degrees of each other, side by side within
   * `alpha_radius`, and end to end within four times that.
   */
  inline bool narrow_in_line(const Rectangle& a, const Rectangle& b,
                             double alpha_radius) {
    auto ua = rectangle_axis(a), ub = rectangle_axis(b);
    if (std::abs(ua[0] * ub[0] + ua[1] * ub[1]) <
        std::cos(10 * std::numbers::pi / 180))
      return false;
    auto ca = mid(a, 0, 2), cb = mid(b, 0, 2);
    double dx = cb[0] - ca[0], dy = cb[1] - ca[1];
    if (std::abs(-dx * ua[1] + dy * ua[0]) > alpha_radius) return false;
    double half_a = ua[2] / 2, half_b = ub[2] / 2;
    return std::abs(dx * ua[0] + dy * ua[1]) - half_a - half_b <=
           4 * alpha_radius;
  }

  // The point of `rings` nearest to `p`.
  inline Point2 nearest_on_border(const std::vector<std::vector<Point2>>& rings,
                                  Point2 p) {
    double near = INFINITY;
    Point2 q = p;
    for (const auto& ring : rings)
      for (size_t i = 0; i < ring.size(); ++i) {
        const auto &a = ring[i], &b = ring[(i + 1) % ring.size()];
        double ex = b[0] - a[0], ey = b[1] - a[1], l2 = ex * ex + ey * ey;
        double s =
            l2 > 0 ? std::clamp(((p[0] - a[0]) * ex + (p[1] - a[1]) * ey) / l2,
                                0.0, 1.0)
                   : 0.0;
        Point2 c{a[0] + s * ex, a[1] + s * ey};
        double d = std::hypot(c[0] - p[0], c[1] - p[1]);
        if (d < near) near = d, q = c;
      }
    return q;
  }

  // `r` with its ends carried along its axis by `carry(end, direction)`.
  template <typename Carry>
  Rectangle carry_ends(const Rectangle& r, Carry carry) {
    auto [ux, uy, l] = rectangle_axis(r);
    double hx = (r[3][0] - r[0][0]) / 2, hy = (r[3][1] - r[0][1]) / 2;
    auto e0 = carry(mid(r, 0, 3), -ux, -uy), e1 = carry(mid(r, 1, 2), ux, uy);
    return Rectangle{
        Point2{e0[0] - hx, e0[1] - hy}, Point2{e1[0] - hx, e1[1] - hy},
        Point2{e1[0] + hx, e1[1] + hy}, Point2{e0[0] + hx, e0[1] + hy}};
  }

  /**
   * `r` carried on along its axis to the footprint's border at either end
   * within `reach` of it: straight on to the border ahead, or as far as the
   * nearest border point lies ahead (a corner). The faces either side of a
   * structure stay one unless it reaches the border.
   *
   * @param rings the footprint's rings, each a closed sequence of (x, y)
   */
  inline Rectangle narrow_to_border(
      const Rectangle& r, const std::vector<std::vector<Point2>>& rings,
      double reach) {
    return carry_ends(r, [&](Point2 p, double dx, double dy) {
      double best_t = INFINITY;
      for (const auto& ring : rings)
        for (size_t i = 0; i < ring.size(); ++i) {
          const auto &a = ring[i], &b = ring[(i + 1) % ring.size()];
          double ex = b[0] - a[0], ey = b[1] - a[1], den = dx * ey - dy * ex;
          if (std::abs(den) > 1e-12) {
            double t = ((a[0] - p[0]) * ey - (a[1] - p[1]) * ex) / den,
                   s = ((a[0] - p[0]) * dy - (a[1] - p[1]) * dx) / den;
            if (t >= 0 && s >= 0 && s <= 1) best_t = std::min(best_t, t);
          }
        }
      Point2 q = nearest_on_border(rings, p);
      if (best_t > reach && std::hypot(q[0] - p[0], q[1] - p[1]) <= reach)
        best_t = (q[0] - p[0]) * dx + (q[1] - p[1]) * dy;
      if (best_t <= 0 || best_t > reach) return p;
      return Point2{p[0] + best_t * dx, p[1] + best_t * dy};
    });
  }

  /**
   * `r` carried on along its axis at either end across another structure
   * whose line, run on by `reach`, it meets, to that structure's far side
   * within `reach`, if that structure's end would reach as far across `r`:
   * walls meeting at a T or a corner join, pieces of one line (within 10
   * degrees) do not.
   */
  inline Rectangle narrow_to_crossing(const Rectangle& r,
                                      const std::vector<Rectangle>& others,
                                      double reach) {
    auto half_width = [](const Rectangle& q) {
      return std::hypot(q[3][0] - q[0][0], q[3][1] - q[0][1]) / 2;
    };
    return carry_ends(r, [&](Point2 p, double dx, double dy) {
      double best = INFINITY;
      for (const auto& o : others) {
        auto [ox, oy, l] = rectangle_axis(o);
        auto a = mid(o, 0, 3);
        double den = dx * oy - dy * ox, ax = a[0] - p[0], ay = a[1] - p[1];
        if (std::abs(den) < std::sin(10 * std::numbers::pi / 180)) continue;
        // as far as `o` would run on to the far side of `r`
        double s = (ax * dy - ay * dx) / den,
               o_reach = reach - half_width(r) / std::abs(den),
               far = (ax * oy - ay * ox) / den + half_width(o) / std::abs(den);
        if (s >= -o_reach && s <= l + o_reach && far > 0 && far <= reach)
          best = std::min(best, far);
      }
      if (best == INFINITY) return p;
      return Point2{p[0] + best * dx, p[1] + best * dy};
    });
  }

  /**
   * Points `xyz` split into straight runs, each in a band as wide as
   * `narrow_outline` allows (sqrt(3) r) with few of all points beside it: a
   * layer of the band's points (heights `epsilon` apart at most), each
   * following the last within `gap`. The line holding most points comes first.
   *
   * @return per run the indices of its points in `xyz`
   */
  inline std::vector<std::vector<size_t>> narrow_runs(
      const std::vector<std::array<double, 3>>& xyz, double r, double gap,
      double epsilon) {
    const double half = std::sqrt(3.0) / 2 * r, bin = r / 2;
    std::vector<std::vector<size_t>> out;
    std::vector<size_t> left(xyz.size());
    for (size_t i = 0; i < xyz.size(); ++i) left[i] = i;
    std::vector<size_t> bins;
    while (left.size() >= 3) {
      double best_t = 0;
      size_t best_n = 0, best_b = 0;
      std::vector<size_t> bin_of(left.size()), best_bin_of(left.size());
      for (int deg = 0; deg < 180; ++deg) {
        double t = deg * std::numbers::pi / 180, nx = -std::sin(t),
               ny = std::cos(t);
        double lo = INFINITY, hi = -INFINITY;
        for (auto i : left) {
          double s = nx * xyz[i][0] + ny * xyz[i][1];
          lo = std::min(lo, s), hi = std::max(hi, s);
        }
        bins.assign(size_t((hi - lo) / bin) + 1, 0);
        bool better = false;
        for (size_t j = 0; j < left.size(); ++j) {
          auto i = left[j];
          size_t b = bin_of[j] = (nx * xyz[i][0] + ny * xyz[i][1] - lo) / bin;
          if (++bins[b] > best_n)
            best_n = bins[b], best_t = t, best_b = b, better = true;
        }
        if (better) std::swap(best_bin_of, bin_of);
      }
      if (best_n < 3) break;
      double ux = std::cos(best_t), uy = std::sin(best_t), nx = -uy, ny = ux;
      std::vector<double> in_bin;
      for (size_t j = 0; j < left.size(); ++j)
        if (best_bin_of[j] == best_b)
          in_bin.push_back(nx * xyz[left[j]][0] + ny * xyz[left[j]][1]);
      std::nth_element(in_bin.begin(), in_bin.begin() + in_bin.size() / 2,
                       in_bin.end());
      double centre = in_bin[in_bin.size() / 2];
      auto off = [&](size_t i) {
        return std::abs(nx * xyz[i][0] + ny * xyz[i][1] - centre);
      };
      auto u = [&](size_t i) { return ux * xyz[i][0] + uy * xyz[i][1]; };
      std::vector<size_t> on, rest, near;
      for (auto i : left) (off(i) < half ? on : rest).push_back(i);
      for (size_t i = 0; i < xyz.size(); ++i)
        if (off(i) >= half && off(i) < 2 * r) near.push_back(i);
      std::sort(on.begin(), on.end(),
                [&](size_t a, size_t b) { return xyz[a][2] < xyz[b][2]; });
      std::vector<size_t> unused;
      for (size_t layer = 0, end = 1; end <= on.size(); ++end) {
        if (end < on.size() && xyz[on[end]][2] - xyz[on[end - 1]][2] <= epsilon)
          continue;
        std::vector<size_t> level(on.begin() + layer, on.begin() + end);
        std::sort(level.begin(), level.end(),
                  [&](size_t a, size_t b) { return u(a) < u(b); });
        for (size_t from = 0, k = 1; k <= level.size(); ++k) {
          if (k < level.size() && u(level[k]) - u(level[k - 1]) <= gap)
            continue;
          double lo = u(level[from]), hi = u(level[k - 1]);
          size_t beside = 0;
          for (auto i : near) beside += u(i) >= lo && u(i) <= hi;
          if (k - from >= 3 && 2 * beside < k - from)
            out.emplace_back(level.begin() + from, level.begin() + k);
          else
            unused.insert(unused.end(), level.begin() + from,
                          level.begin() + k);
          from = k;
        }
        layer = end;
      }
      if (unused.size() < on.size())
        rest.insert(rest.end(), unused.begin(), unused.end());
      left = rest;
    }
    return out;
  }

  // The distance from (x, y) to the nearest of `rings`.
  inline double to_border(const std::vector<std::vector<Point2>>& rings,
                          double x, double y) {
    auto q = nearest_on_border(rings, {x, y});
    return std::hypot(q[0] - x, q[1] - y);
  }

  /**
   * Per point, whether it stands above the roof beside it: a point of
   * another plane (`labels` > 0) within `beside` lies more than `epsilon`
   * lower, and none more than `epsilon` higher. A wall's top, a parapet
   * stand; a facade under the eaves does not. Planes in `steep` don't count as
   * the roof below: their step edges and upper faces are no walls.
   */
  template <typename Points>
  std::vector<bool> standing(const Points& pts, const std::vector<int>& labels,
                             double beside, double epsilon,
                             const std::set<int>& steep = {}) {
    auto cell = [&](double x, double y) {
      return std::array<long long, 2>{(long long)std::floor(x / beside),
                                      (long long)std::floor(y / beside)};
    };
    std::map<std::array<long long, 2>, std::vector<size_t>> grid;
    for (size_t i = 0; i < pts.size(); ++i)
      if (labels[i] > 0) grid[cell(pts[i][0], pts[i][1])].push_back(i);
    std::vector<bool> out(pts.size(), false);
    for (size_t i = 0; i < pts.size(); ++i) {
      auto [cx, cy] = cell(pts[i][0], pts[i][1]);
      bool lower = false, higher = false;
      for (long long dx = -1; dx <= 1 && !higher; ++dx)
        for (long long dy = -1; dy <= 1 && !higher; ++dy) {
          auto it = grid.find({cx + dx, cy + dy});
          if (it == grid.end()) continue;
          for (auto j : it->second) {
            if (labels[j] == labels[i] ||
                std::hypot(pts[j][0] - pts[i][0], pts[j][1] - pts[i][1]) >
                    beside)
              continue;
            double h = pts[i][2] - pts[j][2];
            if (!lower && h > epsilon && !steep.count(labels[j])) lower = true;
            if (h < -epsilon) higher = true;
          }
        }
      out[i] = lower && !higher;
    }
    return out;
  }

  /**
   * `r` run on along its axis over the points `xyz` on its line, within
   * `epsilon` of its plane `z(x, y)`, as long as they follow each other
   * within `gap`: a structure's runs break where its points thin out.
   */
  template <typename Points, typename Height>
  Rectangle narrow_extend(Rectangle r, const Points& xyz, Height z,
                          double epsilon, double gap) {
    auto [ux, uy, l] = rectangle_axis(r);
    double half =
        std::max(gap / 4, std::hypot(r[3][0] - r[0][0], r[3][1] - r[0][1]) / 2);
    Point2 c = mid(r, 0, 2);
    std::vector<double> ahead, behind;
    for (auto& p : xyz) {
      double dx = p[0] - c[0], dy = p[1] - c[1], u = dx * ux + dy * uy;
      if (std::abs(-dx * uy + dy * ux) > half ||
          std::abs(p[2] - z(p[0], p[1])) > epsilon)
        continue;
      if (u > l / 2) ahead.push_back(u - l / 2);
      if (u < -l / 2) behind.push_back(-l / 2 - u);
    }
    auto reach = [&](std::vector<double>& d) {
      std::sort(d.begin(), d.end());
      double at = 0;
      for (auto x : d) {
        if (x - at > gap) break;
        at = x;
      }
      return at;
    };
    double a = reach(ahead), b = reach(behind);
    for (int i : {1, 2}) r[i] = {r[i][0] + a * ux, r[i][1] + a * uy};
    for (int i : {0, 3}) r[i] = {r[i][0] - b * ux, r[i][1] - b * uy};
    return r;
  }

  /**
   * The median height of all points `pts` across `r` in three bands r/2
   * apart: beside it on one side, on it, beside it on the other side; NaN
   * for a band without points.
   */
  template <typename Points>
  std::array<double, 3> narrow_profile(const Points& pts, const Rectangle& r,
                                       double alpha_radius) {
    auto [ux, uy, l] = rectangle_axis(r);
    double bin = alpha_radius / 2;
    Point2 c = mid(r, 0, 3);
    std::array<std::vector<double>, 3> band;
    for (const auto& p : pts) {
      double dx = p[0] - c[0], dy = p[1] - c[1], u = dx * ux + dy * uy,
             v = -dx * uy + dy * ux;
      if (u >= 0 && u <= l && std::abs(v) < 3 * bin)
        band[v < -bin ? 0 : v < bin ? 1 : 2].push_back(p[2]);
    }
    std::array<double, 3> m;
    for (int i = 0; i < 3; ++i) {
      auto& z = band[i];
      if (z.empty()) {
        m[i] = NAN;
        continue;
      }
      std::nth_element(z.begin(), z.begin() + z.size() / 2, z.end());
      m[i] = z[z.size() / 2];
    }
    return m;
  }

  // Whether (x, y) lies in `r`, within `margin`.
  inline bool in_rectangle(const Rectangle& r, double x, double y,
                           double margin) {
    auto [ux, uy, lu] = rectangle_axis(r);
    double vx = r[3][0] - r[0][0], vy = r[3][1] - r[0][1],
           lv = std::hypot(vx, vy);
    double u = (x - r[0][0]) * ux + (y - r[0][1]) * uy,
           v = ((x - r[0][0]) * vx + (y - r[0][1]) * vy) / lv;
    return u >= -margin && u <= lu + margin && v >= -margin && v <= lv + margin;
  }

  /**
   * Narrow structures in a building's segmentation (narrow-planes): wall tops
   * and parapets, each a plane of its own outlined by a rectangle, too narrow
   * for an alpha shape. `labels`, `planes` and `adjacencies` change only if
   * one is found.
   *
   * @return per narrow plane its outline
   */
  template <typename Points, typename Labels, typename Planes,
            typename Adjacencies, typename Footprint, typename Config>
  std::map<int, Rectangle> narrow_structures(const Points& pts, Labels& labels0,
                                             Planes& planes0,
                                             Adjacencies& adjacencies0,
                                             const Footprint& footprint,
                                             const Config& config,
                                             double alpha_radius) {
    using Plane = typename Planes::mapped_type::first_type;
    auto labels = labels0;
    auto planes = planes0;
    auto adjacencies = adjacencies0;
    std::map<int, Rectangle> narrow;
    auto drop = [&](int k) {
      narrow.erase(k);
      planes.erase(k);
      adjacencies.erase(k);
      for (auto& [l, adjacent] : adjacencies) adjacent.erase(k);
      for (auto& l : labels)
        if (l == k) l = 0;
    };
    auto refill = [&](int k) {
      auto& kept = planes[k].second;
      kept.clear();
      for (size_t i = 0; i < labels.size(); ++i)
        if (labels[i] == k) kept.emplace_back(pts[i][0], pts[i][1], pts[i][2]);
    };
    std::set<int> steep;
    for (auto& [k, v] : planes)
      if (std::abs(v.first.c()) <
          std::cos(config.narrow_max_slope * std::numbers::pi / 180) *
              std::hypot(v.first.a(), v.first.b(), v.first.c()))
        steep.insert(k);
    auto stands =
        standing(pts, labels, 2 * alpha_radius, config.plane_epsilon, steep);
    std::map<int, std::array<size_t, 2>> stand_count;
    for (size_t i = 0; i < labels.size(); ++i)
      if (labels[i] > 0 && planes.count(labels[i])) {
        auto& s = stand_count[labels[i]];
        s[0] += stands[i], s[1] += 1;
      }
    auto runs = [&](const std::vector<std::array<double, 3>>& xyz) {
      return narrow_runs(xyz, alpha_radius, 4 * alpha_radius,
                         config.plane_epsilon);
    };
    // the highest roof: the median height of the highest plane
    double roof_top = NAN;
    for (auto& [k, v] : planes) {
      if (v.second.empty()) continue;
      std::vector<double> z;
      for (auto& p : v.second) z.push_back(p.z());
      std::nth_element(z.begin(), z.begin() + z.size() / 2, z.end());
      roof_top = std::fmax(roof_top, z[z.size() / 2]);
    }
    // a run's outline, if it is narrow, stands over the roof on both sides and
    // not far over the highest roof: a wire or a ridge whose roof went
    // undetected is no wall
    auto outline = [&](const std::vector<std::array<double, 3>>& xyz,
                       const std::vector<size_t>& run) {
      std::vector<Point2> xy;
      for (auto i : run) xy.push_back({xyz[i][0], xyz[i][1]});
      auto rect = narrow_outline(xy, alpha_radius);
      if (rect) {
        auto profile = narrow_profile(pts, *rect, alpha_radius);
        if (!(profile[1] - std::fmax(profile[0], profile[2]) >=
              config.narrow_min_height) ||
            profile[1] - roof_top > config.narrow_max_height)
          rect.reset();
      }
      return rect;
    };
    // a standing plane's narrow runs go loose, it keeps the rest if enough
    std::vector<bool> released(labels.size(), false);
    for (auto& [k, s] : stand_count) {
      if (2 * s[0] <= s[1]) continue;
      std::vector<size_t> index;
      std::vector<std::array<double, 3>> xyz;
      for (size_t i = 0; i < labels.size(); ++i)
        if (labels[i] == k) {
          index.push_back(i);
          xyz.push_back({pts[i][0], pts[i][1], pts[i][2]});
        }
      for (auto& run : runs(xyz))
        if (outline(xyz, run))
          for (auto j : run) labels[index[j]] = 0, released[index[j]] = true;
      refill(k);
      if (planes[k].second.size() < size_t(config.min_plane_points)) drop(k);
    }
    std::vector<size_t> loose;
    std::vector<std::array<double, 3>> xyz;
    for (size_t i = 0; i < labels.size(); ++i)
      if (labels[i] == 0 && (stands[i] || released[i])) {
        loose.push_back(i);
        xyz.push_back({pts[i][0], pts[i][1], pts[i][2]});
      }
    int id = 0;  // past every plane, also those dropped above
    for (auto& [k, v] : planes0) id = std::max(id, k);
    for (auto& run : runs(xyz)) {
      auto rect = outline(xyz, run);
      if (!rect) continue;
      ++id;
      for (auto i : run) {
        planes[id].second.emplace_back(xyz[i][0], xyz[i][1], xyz[i][2]);
        labels[loose[i]] = id;
      }
      narrow[id] = *rect;
    }
    auto level_across = [&](int k) {
      std::vector<std::array<double, 3>> xyz;
      for (auto& p : planes[k].second) xyz.push_back({p.x(), p.y(), p.z()});
      auto c = narrow_plane(xyz, narrow[k]);
      planes[k].first = Plane(c[0], c[1], c[2], c[3]);
    };
    for (auto& [k, r] : narrow) level_across(k);
    auto height = [&](int k, double x, double y) {
      return plane_z(planes[k].first, x, y);
    };
    // pieces of one narrow structure at one height become one plane
    for (bool merged = true; merged;) {
      merged = false;
      for (auto a = narrow.begin(); a != narrow.end() && !merged; ++a) {
        for (auto b = std::next(a); b != narrow.end(); ++b) {
          double x = (a->second[0][0] + b->second[2][0]) / 2,
                 y = (a->second[0][1] + b->second[2][1]) / 2;
          if (!narrow_in_line(a->second, b->second, alpha_radius) ||
              std::abs(height(a->first, x, y) - height(b->first, x, y)) >
                  config.plane_epsilon)
            continue;
          std::vector<Point2> xy;
          for (auto k : {a->first, b->first})
            for (auto& p : planes[k].second) xy.push_back({p.x(), p.y()});
          auto r = narrow_outline(xy, alpha_radius);
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
      }
    }
    std::vector<std::vector<Point2>> rings(1);
    for (auto& p : footprint) rings[0].push_back({p[0], p[1]});
    for (auto& hole : footprint.interior_rings()) {
      rings.emplace_back();
      for (auto& p : hole) rings.back().push_back({p[0], p[1]});
    }
    std::vector<typename Points::value_type> unsegmented;
    for (size_t i = 0; i < labels.size(); ++i)
      if (labels[i] == 0) unsegmented.push_back(pts[i]);
    std::vector<int> short_or_on_border;
    for (auto& [k, r] : narrow) {
      r = narrow_extend(
          r, unsegmented, [&](double x, double y) { return height(k, x, y); },
          config.plane_epsilon, 2 * alpha_radius);
      r = narrow_to_border(r, rings, 4 * alpha_radius);
      auto inside = [&](Point2 p) {
        return to_border(rings, p[0], p[1]) > 2 * alpha_radius;
      };
      double length = rectangle_axis(r)[2];
      auto profile = narrow_profile(pts, r, alpha_radius);
      bool isolated =
          inside(mid(r, 0, 3)) && inside(mid(r, 1, 2)) &&
          std::abs(profile[0] - profile[2]) <= config.narrow_min_height;
      if (length < 4 * alpha_radius ||
          (isolated && length < config.narrow_min_length) ||
          !(inside(mid(r, 0, 2)) || inside(mid(r, 0, 3)) ||
            inside(mid(r, 1, 2))))
        short_or_on_border.push_back(k);
    }
    for (auto k : short_or_on_border) drop(k);
    if (narrow.empty()) return narrow;
    // walls meeting at a T or a corner join
    for (auto& [k, r] : narrow) {
      std::vector<Rectangle> others;
      for (auto& [l, o] : narrow)
        if (l != k) others.push_back(o);
      r = narrow_to_crossing(r, others, 4 * alpha_radius);
    }
    // leftover points go back to their plane, if enough for a dropped one
    std::map<int, size_t> back;
    for (size_t i = 0; i < labels.size(); ++i)
      if (labels[i] == 0 && labels0[i] > 0 && planes0.count(labels0[i]))
        ++back[labels0[i]];
    for (auto [k, n] : back) {
      bool dropped = !planes.count(k);
      if (dropped && n < size_t(config.min_plane_points)) continue;
      for (size_t i = 0; i < labels.size(); ++i)
        if (labels[i] == 0 && labels0[i] == k) labels[i] = k;
      planes[k].first = planes0.at(k).first;
      refill(k);
      if (!dropped) continue;
      if (adjacencies0.count(k))
        for (auto& [l, n] : adjacencies0.at(k))
          if (planes.count(l)) adjacencies[k][l] = n;
      for (auto& [l, adjacent] : adjacencies0)
        if (adjacent.count(k) && planes.count(l))
          adjacencies[l][k] = adjacent.at(k);
    }
    labels0 = std::move(labels);
    planes0 = std::move(planes);
    adjacencies0 = std::move(adjacencies);
    return narrow;
  }

}  // namespace roofer::reconstruction
