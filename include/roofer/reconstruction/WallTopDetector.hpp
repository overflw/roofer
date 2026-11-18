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

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <unordered_map>
#include <vector>

namespace roofer::reconstruction {

  struct WallTopParams {
    float min_height = 0.3F;  // above the roof plane under it, in metres
    float max_height = 1.5F;
    float min_length = 5.0F;  // shortest wall, in metres
    float max_slope = 10.0F;  // steepest roof plane a wall is searched on,
                              // in degrees
    // the alpha shape's radius, the narrowest structure the reconstruction
    // resolves; the lengths below that are not a wall's own scale with it
    float alpha_radius = 0.5F;
  };

  struct WallTop {
    std::vector<size_t> points;
    // the plane a*x + b*y + c*z + d = 0: sloped along the wall, level across
    double a, b, c, d;
    // the wall's centre line, from (x0, y0) to (x1, y1)
    double x0, y0, x1, y1;
    double width;  // of the strip modelled on it
  };

  /**
   * Carry each end of `w` on to the roof's border where that lies within
   * `reach`: straight on along the wall, else to the nearest point of
   * it, a corner the wall runs into. A wall stops a little short of the
   * border in the points; unless it reaches it, the faces on either side of
   * it stay one.
   *
   * @param rings the footprint's rings, each a closed sequence of (x, y)
   */
  inline void snap_wall_top(
      WallTop& w, const std::vector<std::vector<std::array<double, 2>>>& rings,
      double reach) {
    auto carry = [&](double px, double py, double dx, double dy) {
      double best_t = INFINITY, near = INFINITY, nx = px, ny = py;
      for (const auto& r : rings) {
        for (size_t i = 0; i < r.size(); ++i) {
          const auto& a = r[i];
          const auto& b = r[(i + 1) % r.size()];
          double ex = b[0] - a[0], ey = b[1] - a[1];
          // the ray p + t * d against the edge a + s * e
          double den = dx * ey - dy * ex;
          if (std::abs(den) > 1e-12) {
            double t = ((a[0] - px) * ey - (a[1] - py) * ex) / den,
                   s = ((a[0] - px) * dy - (a[1] - py) * dx) / den;
            if (t >= 0 && s >= 0 && s <= 1) best_t = std::min(best_t, t);
          }
          double l2 = ex * ex + ey * ey,
                 s = l2 > 0 ? std::clamp(
                                  ((px - a[0]) * ex + (py - a[1]) * ey) / l2,
                                  0.0, 1.0)
                            : 0.0,
                 qx = a[0] + s * ex, qy = a[1] + s * ey,
                 d = std::hypot(qx - px, qy - py);
          if (d < near) near = d, nx = qx, ny = qy;
        }
      }
      // past the border by the strip's width, so that it cuts through
      auto past = [&](double qx, double qy) {
        double l = std::hypot(qx - px, qy - py);
        return l == 0 ? std::array<double, 2>{qx, qy}
                      : std::array<double, 2>{qx + (qx - px) / l * w.width,
                                              qy + (qy - py) / l * w.width};
      };
      if (best_t <= reach) return past(px + best_t * dx, py + best_t * dy);
      if (near <= reach) return past(nx, ny);
      return std::array<double, 2>{px, py};
    };
    double l = std::hypot(w.x1 - w.x0, w.y1 - w.y0);
    if (l == 0) return;
    double dx = (w.x1 - w.x0) / l, dy = (w.y1 - w.y0) / l;
    auto p0 = carry(w.x0, w.y0, -dx, -dy), p1 = carry(w.x1, w.y1, dx, dy);
    w.x0 = p0[0], w.y0 = p0[1], w.x1 = p1[0], w.y1 = p1[1];
  }

  /**
   * Walls standing on a roof: a firewall or parapet across a flat roof is a
   * few decimetres wide, too narrow for the region grower, so its points stay
   * unsegmented or fall into narrow planes of their own, and the roof plane
   * runs on underneath. Its points stand in a straight line above the roof
   * plane next to them.
   *
   * In units of the alpha radius r: takes the points unsegmented (`labels`
   * 0) or in a plane on top of walls (most of its points have another
   * plane's points below them within BESIDE = 2r) lying `min_height` to
   * `max_height` above the plane of the nearest point of another plane
   * within CELL = 3r, and more than EDGE = 2r inside the footprint, where
   * parapets run along the walls. Searches straight lines in them (offset
   * bins of r/2 over 1° steps) and keeps a line at least `min_length` long
   * with enough points in 70 % of its metres, no broad block (PV rows)
   * beside it and the roof not rising across it as much on one side as it
   * falls on the other (the flank of a ridge). Its strip is r wide, and an end within 5r of the footprint's
   * border runs on to it (`snap_wall_top`).
   *
   * @param planes plane coefficients (a, b, c, d) per label
   * @param rings the footprint's rings, each a closed sequence of (x, y)
   * @return per wall its points and a plane through them
   */
  template <typename Points>  // indexable (x, y, z), as std::vector<arr3f>
  std::vector<WallTop> detect_wall_tops(
      const Points& pts, const std::vector<int>& labels,
      const std::map<int, std::array<double, 4>>& planes,
      const std::vector<std::vector<std::array<double, 2>>>& rings,
      const WallTopParams& params = WallTopParams()) {
    constexpr double FILL = 0.7, THIN = 0.5;
    const double r = params.alpha_radius, BIN = r / 2, BESIDE = 2 * r,
                 EDGE = 2 * r, CELL = 3 * r;
    const double max_tilt = std::tan(params.max_slope * M_PI / 180.0);
    std::vector<WallTop> out;
    if (pts.empty()) return out;

    auto cell = [](long long cx, long long cy) {
      return (static_cast<unsigned long long>(cx) << 32) ^
             (static_cast<unsigned long long>(cy) & 0xffffffffULL);
    };
    auto key = [&](double x, double y) {
      return cell(std::floor(x / CELL), std::floor(y / CELL));
    };
    // A plane is a wall's top where most of its points have one of another
    // plane below them, as low as a wall is high, within BESIDE: a wall has
    // roof on both sides, a higher roof only along its edge.
    std::unordered_map<unsigned long long, std::vector<size_t>> near;
    auto near_key = [&](double x, double y) {
      return cell(std::floor(x / BESIDE), std::floor(y / BESIDE));
    };
    for (size_t i = 0; i < pts.size(); ++i)
      if (labels[i] > 0) near[near_key(pts[i][0], pts[i][1])].push_back(i);
    std::map<int, std::array<size_t, 2>> over;  // points over another, all
    for (size_t i = 0; i < pts.size(); ++i) {
      if (labels[i] <= 0) continue;
      long long cx = std::floor(pts[i][0] / BESIDE),
                cy = std::floor(pts[i][1] / BESIDE);
      bool above = false;
      for (long long dx = -1; dx <= 1 && !above; ++dx)
        for (long long dy = -1; dy <= 1 && !above; ++dy) {
          auto it = near.find(cell(cx + dx, cy + dy));
          if (it == near.end()) continue;
          for (auto j : it->second) {
            double h = pts[i][2] - pts[j][2];
            if (labels[j] != labels[i] && h > params.min_height &&
                h < params.max_height &&
                std::hypot(pts[j][0] - pts[i][0], pts[j][1] - pts[i][1]) <
                    BESIDE) {
              above = true;
              break;
            }
          }
        }
      auto& o = over[labels[i]];
      o[0] += above, o[1] += 1;
    }
    std::map<int, bool> narrow;
    for (auto& [l, o] : over) narrow[l] = o[0] > o[1] / 2;
    auto is_roof = [&](size_t i) {
      return labels[i] > 0 && !narrow[labels[i]];
    };
    auto to_border = [&](double x, double y) {
      double best = INFINITY;
      for (const auto& r : rings)
        for (size_t k = 0; k < r.size(); ++k) {
          const auto &a = r[k], &b = r[(k + 1) % r.size()];
          double ex = b[0] - a[0], ey = b[1] - a[1], l2 = ex * ex + ey * ey;
          double s = l2 > 0
                         ? std::clamp(((x - a[0]) * ex + (y - a[1]) * ey) / l2,
                                      0.0, 1.0)
                         : 0.0;
          best =
              std::min(best, std::hypot(a[0] + s * ex - x, a[1] + s * ey - y));
        }
      return best;
    };

    std::unordered_map<unsigned long long, std::vector<size_t>> grid;
    std::unordered_map<unsigned long long, size_t> occupied;
    for (size_t i = 0; i < pts.size(); ++i) {
      if (is_roof(i)) grid[key(pts[i][0], pts[i][1])].push_back(i);
      ++occupied[key(pts[i][0], pts[i][1])];
    }
    // points a metre of wall needs on its 2 * BIN wide band
    const double density = double(pts.size()) / (occupied.size() * CELL * CELL);
    const double dense = std::max(2.0, 0.3 * density * 2 * BIN);

    std::vector<size_t> cand, near_border;  // raised, inside EDGE or not
    for (size_t i = 0; i < pts.size(); ++i) {
      if (is_roof(i)) continue;
      const double x = pts[i][0], y = pts[i][1];
      long long cx = std::floor(x / CELL), cy = std::floor(y / CELL);
      double best = CELL * CELL;
      int label = 0;
      for (long long dx = -1; dx <= 1; ++dx)
        for (long long dy = -1; dy <= 1; ++dy) {
          auto it = grid.find(cell(cx + dx, cy + dy));
          if (it == grid.end()) continue;
          for (auto j : it->second) {
            double d2 = (pts[j][0] - x) * (pts[j][0] - x) +
                        (pts[j][1] - y) * (pts[j][1] - y);
            if (d2 < best) best = d2, label = labels[j];
          }
        }
      auto pl = planes.find(label);
      if (pl == planes.end()) continue;
      const auto& [a, b, c, d] = pl->second;
      if (std::abs(c) < 1e-9 || std::hypot(a, b) / std::abs(c) > max_tilt)
        continue;
      const double h = pts[i][2] + (a * x + b * y + d) / c;
      if (h > params.min_height && h < params.max_height)
        (to_border(x, y) > EDGE ? cand : near_border).push_back(i);
    }

    std::vector<double> dirs_x(180), dirs_y(180);
    for (int t = 0; t < 180; ++t) {
      dirs_x[t] = std::cos(t * M_PI / 180),
      dirs_y[t] = std::sin(t * M_PI / 180);
    }
    while (cand.size() >= params.min_length * FILL) {
      // the angle and offset bin holding most candidates
      int best_t = 0;
      long long best_bin = 0;
      size_t best_n = 0;
      for (int t = 0; t < 180; ++t) {
        std::unordered_map<long long, size_t> bins;
        for (auto i : cand) {
          double s = -dirs_y[t] * pts[i][0] + dirs_x[t] * pts[i][1];
          size_t n = ++bins[static_cast<long long>(std::floor(s / BIN))];
          if (n > best_n)
            best_n = n, best_t = t, best_bin = std::floor(s / BIN);
        }
      }
      const double nx = -dirs_y[best_t], ny = dirs_x[best_t],
                   ux = dirs_x[best_t], uy = dirs_y[best_t];
      std::vector<double> in_bin;
      for (auto i : cand) {
        double s = nx * pts[i][0] + ny * pts[i][1];
        if (std::floor(s / BIN) == best_bin) in_bin.push_back(s);
      }
      std::nth_element(in_bin.begin(), in_bin.begin() + in_bin.size() / 2,
                       in_bin.end());
      const double centre = in_bin[in_bin.size() / 2];

      std::vector<size_t> on;
      double umin = INFINITY, umax = -INFINITY;
      for (auto i : cand) {
        if (std::abs(nx * pts[i][0] + ny * pts[i][1] - centre) < BIN) {
          on.push_back(i);
          double u = ux * pts[i][0] + uy * pts[i][1];
          umin = std::min(umin, u), umax = std::max(umax, u);
        }
      }
      // the longest stretch with `dense` points in FILL of its metres
      std::vector<size_t> per_m(size_t(umax - umin) + 1, 0);
      for (auto i : on) ++per_m[size_t(ux * pts[i][0] + uy * pts[i][1] - umin)];
      std::vector<int> ok(per_m.size() + 1, 0);
      for (size_t m = 0; m < per_m.size(); ++m)
        ok[m + 1] = ok[m] + (per_m[m] >= dense);
      long best_i = -1, best_j = -1;
      for (size_t i = 0; i < per_m.size(); ++i) {
        if (per_m[i] < dense) continue;
        for (size_t j = per_m.size() - 1;
             j + 1 >= i + size_t(params.min_length) && j >= i; --j) {
          if (per_m[j] >= dense && ok[j + 1] - ok[i] >= FILL * (j - i + 1)) {
            if (long(j - i) > best_j - best_i) best_i = i, best_j = j;
            break;
          }
          if (j == 0) break;
        }
      }

      std::vector<size_t> at, rest;
      if (best_i >= 0) {
        const double lo = umin + best_i, hi = umin + best_j + 1;
        size_t beside = 0;
        for (auto i : cand) {
          double u = ux * pts[i][0] + uy * pts[i][1],
                 off = std::abs(nx * pts[i][0] + ny * pts[i][1] - centre);
          if (u < lo || u > hi)
            rest.push_back(i);
          else if (off < BIN)
            at.push_back(i);
          else {
            beside += off < BESIDE;
            rest.push_back(i);
          }
        }
        // a broad block cut short by EDGE is no wall either
        for (auto i : near_border) {
          double u = ux * pts[i][0] + uy * pts[i][1],
                 off = std::abs(nx * pts[i][0] + ny * pts[i][1] - centre);
          beside += u >= lo && u <= hi && off >= BIN && off < BESIDE;
        }
        // a wall stands over the roof on one side at least; a line on the
        // flank of a broader ridge has it rising across as much as falling
        std::array<std::vector<double>, 3> band;
        for (size_t i = 0; i < pts.size(); ++i) {
          double u = ux * pts[i][0] + uy * pts[i][1],
                 v = nx * pts[i][0] + ny * pts[i][1] - centre;
          if (u >= lo && u <= hi && std::abs(v) < 3 * BIN)
            band[v < -BIN ? 0 : v < BIN ? 1 : 2].push_back(pts[i][2]);
        }
        auto median = [](std::vector<double>& z) {
          std::nth_element(z.begin(), z.begin() + z.size() / 2, z.end());
          return z[z.size() / 2];
        };
        bool ramp = false;
        if (!band[0].empty() && !band[1].empty() && !band[2].empty()) {
          double top = median(band[1]), rise = top - median(band[0]),
                 on = median(band[2]) - top;
          ramp = rise * on > 0 &&
                 std::min(std::abs(rise), std::abs(on)) >
                     params.min_height / 2 &&
                 std::max(std::abs(rise), std::abs(on)) <
                     2 * std::min(std::abs(rise), std::abs(on));
        }
        if (beside <= THIN * at.size() && !ramp) {
          // z = alpha + g * u, least squares
          double su = 0, sz = 0, suu = 0, suz = 0, n = at.size();
          for (auto i : at) {
            double u = ux * pts[i][0] + uy * pts[i][1];
            su += u, sz += pts[i][2], suu += u * u, suz += u * pts[i][2];
          }
          double den = n * suu - su * su,
                 g = std::abs(den) > 1e-9 ? (n * suz - su * sz) / den : 0,
                 alpha = (sz - g * su) / n;
          out.push_back({at, -g * ux, -g * uy, 1.0, -alpha,
                         centre * nx + lo * ux, centre * ny + lo * uy,
                         centre * nx + hi * ux, centre * ny + hi * uy, r});
          snap_wall_top(out.back(), rings, 5 * r);
        }
      }
      if (at.empty()) break;
      cand = rest;
    }
    return out;
  }

  // The corners of the strip on top of `w`.
  inline std::array<std::array<double, 2>, 4> wall_top_corners(
      const WallTop& w) {
    double ux = w.x1 - w.x0, uy = w.y1 - w.y0, l = std::hypot(ux, uy);
    double nx = -uy / l * w.width / 2, ny = ux / l * w.width / 2;
    return {{{w.x0 + nx, w.y0 + ny},
             {w.x1 + nx, w.y1 + ny},
             {w.x1 - nx, w.y1 - ny},
             {w.x0 - nx, w.y0 - ny}}};
  }

  // The lines `w` puts in the arrangement: its strip's long edges, run on by
  // `extension` like the regularised lines, and its ends.
  inline std::vector<std::array<std::array<double, 2>, 2>> wall_top_edges(
      const WallTop& w, double extension) {
    auto c = wall_top_corners(w);
    double ux = w.x1 - w.x0, uy = w.y1 - w.y0, l = std::hypot(ux, uy);
    double ex = ux / l * extension, ey = uy / l * extension;
    return {{{{c[0][0] - ex, c[0][1] - ey}, {c[1][0] + ex, c[1][1] + ey}}},
            {{{c[3][0] - ex, c[3][1] - ey}, {c[2][0] + ex, c[2][1] + ey}}},
            {{c[0], c[3]}},
            {{c[1], c[2]}}};
  }

  // Whether (x, y) lies on the strip on top of `w`.
  inline bool in_wall_top(const WallTop& w, double x, double y) {
    double ux = w.x1 - w.x0, uy = w.y1 - w.y0, l = std::hypot(ux, uy);
    if (l == 0) return false;
    double u = ((x - w.x0) * ux + (y - w.y0) * uy) / l,
           v = ((x - w.x0) * -uy + (y - w.y0) * ux) / l;
    return u >= 0 && u <= l && std::abs(v) <= w.width;
  }

}  // namespace roofer::reconstruction
