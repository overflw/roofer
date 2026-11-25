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
#include <map>
#include <optional>
#include <set>
#include <vector>

namespace roofer::reconstruction {

  using Point2 = std::array<double, 2>;
  using Rectangle = std::array<Point2, 4>;

  /**
   * The outline of a narrow plane: a wall's top, a parapet, a pipe. An alpha
   * shape of `alpha_radius` breaks a region narrower than about twice that
   * radius into fragments, so a plane whose points spread across it by less
   * than half the radius (standard deviation) and run on at least four times
   * as long takes the rectangle around them along their main direction.
   *
   * @param xy the plane's points in plan
   * @return the rectangle, corners in order, or nothing for a plane that is
   * not narrow
   */
  inline std::optional<Rectangle> narrow_outline(const std::vector<Point2>& xy,
                                                 double alpha_radius) {
    if (xy.size() < 3) return std::nullopt;
    double mx = 0, my = 0;
    for (auto& p : xy) mx += p[0], my += p[1];
    mx /= xy.size(), my /= xy.size();
    double sxx = 0, syy = 0, sxy = 0;
    for (auto& p : xy) {
      double dx = p[0] - mx, dy = p[1] - my;
      sxx += dx * dx, syy += dy * dy, sxy += dx * dy;
    }
    sxx /= xy.size(), syy /= xy.size(), sxy /= xy.size();
    double t = std::atan2(2 * sxy, sxx - syy) / 2, ux = std::cos(t),
           uy = std::sin(t);
    double umin = INFINITY, umax = -INFINITY, vmin = INFINITY, vmax = -INFINITY,
           vv = 0;
    for (auto& p : xy) {
      double u = (p[0] - mx) * ux + (p[1] - my) * uy,
             v = -(p[0] - mx) * uy + (p[1] - my) * ux;
      umin = std::min(umin, u), umax = std::max(umax, u);
      vmin = std::min(vmin, v), vmax = std::max(vmax, v), vv += v * v;
    }
    double spread = std::sqrt(vv / xy.size());
    if (spread >= alpha_radius / 2 || umax - umin < 2 * alpha_radius ||
        umax - umin < 4 * (vmax - vmin))
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

  /**
   * The lines a narrow plane's outline puts in the arrangement: its long
   * sides, carried on like the regularised lines by `extension`, but no
   * farther than the outline is long, and its ends.
   */
  inline std::vector<std::array<Point2, 2>> narrow_edges(const Rectangle& r,
                                                         double extension) {
    double ux = r[1][0] - r[0][0], uy = r[1][1] - r[0][1],
           l = std::hypot(ux, uy);
    extension = std::min(extension, l);
    double ex = ux / l * extension, ey = uy / l * extension;
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
    double ux = r[1][0] - r[0][0], uy = r[1][1] - r[0][1],
           l = std::hypot(ux, uy);
    ux /= l, uy /= l;
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
    auto axis = [](const Rectangle& r) {
      double ux = r[1][0] - r[0][0], uy = r[1][1] - r[0][1],
             l = std::hypot(ux, uy);
      return std::array<double, 2>{ux / l, uy / l};
    };
    auto centre = [](const Rectangle& r) {
      return Point2{(r[0][0] + r[2][0]) / 2, (r[0][1] + r[2][1]) / 2};
    };
    auto ua = axis(a), ub = axis(b);
    if (std::abs(ua[0] * ub[0] + ua[1] * ub[1]) < std::cos(10 * M_PI / 180))
      return false;
    auto ca = centre(a), cb = centre(b);
    double dx = cb[0] - ca[0], dy = cb[1] - ca[1];
    if (std::abs(-dx * ua[1] + dy * ua[0]) > alpha_radius) return false;
    double half_a = std::hypot(a[1][0] - a[0][0], a[1][1] - a[0][1]) / 2,
           half_b = std::hypot(b[1][0] - b[0][0], b[1][1] - b[0][1]) / 2;
    return std::abs(dx * ua[0] + dy * ua[1]) - half_a - half_b <=
           4 * alpha_radius;
  }

  /**
   * `r` carried on along its axis to the footprint's border, at either end
   * that lies within `reach` of it: straight on where the border lies ahead,
   * else as far on as its nearest point lies ahead, a corner the structure
   * runs into. Points next
   * to a facade are few, so a structure stops a little short of the border
   * in them; unless it reaches it, the faces on either side of it stay one.
   *
   * @param rings the footprint's rings, each a closed sequence of (x, y)
   */
  inline Rectangle narrow_to_border(
      Rectangle r, const std::vector<std::vector<Point2>>& rings,
      double reach) {
    auto carry = [&](Point2 p, double dx, double dy) {
      double best_t = INFINITY, near = INFINITY;
      Point2 q = p;
      for (const auto& ring : rings)
        for (size_t i = 0; i < ring.size(); ++i) {
          const auto &a = ring[i], &b = ring[(i + 1) % ring.size()];
          double ex = b[0] - a[0], ey = b[1] - a[1], den = dx * ey - dy * ex;
          if (std::abs(den) > 1e-12) {
            double t = ((a[0] - p[0]) * ey - (a[1] - p[1]) * ex) / den,
                   s = ((a[0] - p[0]) * dy - (a[1] - p[1]) * dx) / den;
            if (t >= 0 && s >= 0 && s <= 1) best_t = std::min(best_t, t);
          }
          double l2 = ex * ex + ey * ey,
                 s = l2 > 0
                         ? std::clamp(
                               ((p[0] - a[0]) * ex + (p[1] - a[1]) * ey) / l2,
                               0.0, 1.0)
                         : 0.0;
          Point2 c{a[0] + s * ex, a[1] + s * ey};
          double d = std::hypot(c[0] - p[0], c[1] - p[1]);
          if (d < near) near = d, q = c;
        }
      if (best_t > reach && near <= reach)
        best_t = (q[0] - p[0]) * dx + (q[1] - p[1]) * dy;
      if (best_t <= 0 || best_t > reach) return p;
      return Point2{p[0] + best_t * dx, p[1] + best_t * dy};
    };
    double ux = r[1][0] - r[0][0], uy = r[1][1] - r[0][1],
           l = std::hypot(ux, uy);
    ux /= l, uy /= l;
    Point2 m0{(r[0][0] + r[3][0]) / 2, (r[0][1] + r[3][1]) / 2},
        m1{(r[1][0] + r[2][0]) / 2, (r[1][1] + r[2][1]) / 2};
    double hx = (r[3][0] - r[0][0]) / 2, hy = (r[3][1] - r[0][1]) / 2;
    auto e0 = carry(m0, -ux, -uy), e1 = carry(m1, ux, uy);
    return Rectangle{
        Point2{e0[0] - hx, e0[1] - hy}, Point2{e1[0] - hx, e1[1] - hy},
        Point2{e1[0] + hx, e1[1] + hy}, Point2{e0[0] + hx, e0[1] + hy}};
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
      double best_t = 0, best_s = 0;
      size_t best_n = 0;
      for (int deg = 0; deg < 180; ++deg) {
        double t = deg * M_PI / 180, nx = -std::sin(t), ny = std::cos(t);
        double lo = INFINITY, hi = -INFINITY;
        for (auto i : left) {
          double s = nx * xyz[i][0] + ny * xyz[i][1];
          lo = std::min(lo, s), hi = std::max(hi, s);
        }
        bins.assign(size_t((hi - lo) / bin) + 1, 0);
        for (auto i : left) {
          size_t b = (nx * xyz[i][0] + ny * xyz[i][1] - lo) / bin;
          if (++bins[b] > best_n)
            best_n = bins[b], best_t = t, best_s = lo + (b + 0.5) * bin;
        }
      }
      if (best_n < 3) break;
      double ux = std::cos(best_t), uy = std::sin(best_t), nx = -uy, ny = ux;
      std::vector<double> in_bin;
      for (auto i : left) {
        double s = nx * xyz[i][0] + ny * xyz[i][1];
        if (std::abs(s - best_s) <= bin / 2) in_bin.push_back(s);
      }
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
    double best = INFINITY;
    for (const auto& ring : rings)
      for (size_t i = 0; i < ring.size(); ++i) {
        const auto &a = ring[i], &b = ring[(i + 1) % ring.size()];
        double ex = b[0] - a[0], ey = b[1] - a[1], l2 = ex * ex + ey * ey;
        double s = l2 > 0 ? std::clamp(((x - a[0]) * ex + (y - a[1]) * ey) / l2,
                                       0.0, 1.0)
                          : 0.0;
        best = std::min(best, std::hypot(a[0] + s * ex - x, a[1] + s * ey - y));
      }
    return best;
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
      for (long long dx = -1; dx <= 1; ++dx)
        for (long long dy = -1; dy <= 1; ++dy) {
          auto it = grid.find({cx + dx, cy + dy});
          if (it == grid.end()) continue;
          for (auto j : it->second) {
            if (labels[j] == labels[i] ||
                std::hypot(pts[j][0] - pts[i][0], pts[j][1] - pts[i][1]) >
                    beside)
              continue;
            double h = pts[i][2] - pts[j][2];
            lower |= h > epsilon && !steep.count(labels[j]);
            higher |= h < -epsilon;
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
    double ux = r[1][0] - r[0][0], uy = r[1][1] - r[0][1],
           l = std::hypot(ux, uy);
    ux /= l, uy /= l;
    double half =
        std::max(gap / 4, std::hypot(r[3][0] - r[0][0], r[3][1] - r[0][1]) / 2);
    Point2 c{(r[0][0] + r[2][0]) / 2, (r[0][1] + r[2][1]) / 2};
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
    double ux = r[1][0] - r[0][0], uy = r[1][1] - r[0][1],
           l = std::hypot(ux, uy), bin = alpha_radius / 2;
    ux /= l, uy /= l;
    Point2 c{(r[0][0] + r[3][0]) / 2, (r[0][1] + r[3][1]) / 2};
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

  /**
   * Whether a profile (`narrow_profile`) lies on the flank of a broader
   * ridge: across it, the roof rises on one side about as much as it falls
   * on the other (both by more than `step`, the larger less than twice the
   * smaller). A wall stands over the roof on one side at least.
   */
  inline bool narrow_on_flank(const std::array<double, 3>& m, double step) {
    double rise = m[1] - m[0], on = m[2] - m[1];
    if (std::isnan(rise) || std::isnan(on)) return false;
    return rise * on > 0 && std::min(std::abs(rise), std::abs(on)) > step &&
           std::max(std::abs(rise), std::abs(on)) <
               2 * std::min(std::abs(rise), std::abs(on));
  }

  // Whether (x, y) lies in `r`, within `margin`.
  inline bool in_rectangle(const Rectangle& r, double x, double y,
                           double margin) {
    double ux = r[1][0] - r[0][0], uy = r[1][1] - r[0][1],
           lu = std::hypot(ux, uy), vx = r[3][0] - r[0][0],
           vy = r[3][1] - r[0][1], lv = std::hypot(vx, vy);
    double u = ((x - r[0][0]) * ux + (y - r[0][1]) * uy) / lu,
           v = ((x - r[0][0]) * vx + (y - r[0][1]) * vy) / lv;
    return u >= -margin && u <= lu + margin && v >= -margin && v <= lv + margin;
  }

}  // namespace roofer::reconstruction
