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
#include <functional>
#include <memory>
#include <set>
#include <roofer/common/Raster.hpp>
#include <roofer/common/datastructures.hpp>
#include <roofer/reconstruction/cgal_shared_definitions.hpp>
#include <roofer/common/ConfigField.hpp>

namespace roofer::reconstruction {

#define ROOFER_ARRANGEMENT_OPTIMISER_FIELDS(X)                                \
  X(float, complexity_factor, 0.888F,                                         \
    "Balance between data fidelity (1) and smoothness (0).",                  \
    config::in_range(0.0F, 1.0F), public_)                                    \
  X(bool, preset_labels, false, "Preset graph-cut labels.",                   \
    config::no_validation<bool>(), internal)                                  \
  X(bool, normalise, false, "Normalise graph-cut costs.",                     \
    config::no_validation<bool>(), internal)                                  \
  X(int, graph_cut_impl, 0,                                                   \
    "Graph-cut backend (experimental). 0 uses Boost adjacency list "          \
    "(default). 1 uses "                                                      \
    "Boost compressed sparse row (slower, uses less memory). 2 uses MaxFlow " \
    "(faster but may consume a lot more memory).",                            \
    config::in_range(0, 2), public_)                                          \
  X(bool, label_ground_outside_footprint, true,                               \
    "Label ground outside the footprint.", config::no_validation<bool>(),     \
    internal)                                                                 \
  X(bool, use_ground, true, "Use ground labels (pipeline-derived).",          \
    config::no_validation<bool>(), internal)
  struct ArrangementOptimiserConfig {
    using Self = ArrangementOptimiserConfig;
    ROOFER_CONFIG_MEMBERS(ROOFER_ARRANGEMENT_OPTIMISER_FIELDS)

    [[nodiscard]] constexpr float data_weight() const {
      return complexity_factor;
    }
    [[nodiscard]] constexpr float smoothness_weight() const {
      return 1.0F - complexity_factor;
    }
  };
#undef ROOFER_ARRANGEMENT_OPTIMISER_FIELDS

  struct ArrangementOptimiserInterface {
    // add_input("arrangement", typeid(Arrangement_2));
    // add_input("heightfield", typeid(RasterTools::Raster));
    // add_input("pts_per_roofplane", typeid(IndexedPlanesWithPoints ));
    // add_input("ground_pts_per_roofplane", typeid(IndexedPlanesWithPoints ));
    // add_output("arrangement", typeid(Arrangement_2));
    // add_vector_output("groundparts", typeid(LinearRing));

    virtual ~ArrangementOptimiserInterface() = default;
    virtual void compute(
        Arrangement_2& arrangement, const RasterTools::Raster& heightfield,
        const IndexedPlanesWithPoints& roof_planes,
        const IndexedPlanesWithPoints& ground_planes,
        ArrangementOptimiserConfig config = ArrangementOptimiserConfig(),
        // whether a roof plane may label a face; one that may not costs it
        // more than any other label
        std::function<bool(Arrangement_2::Face_handle, int)> allowed = {}) = 0;
  };

  std::vector<LinearRing> arr2polygons(Arrangement_2& arr);

  /**
   * After the graph cut: a roof plane in `narrow` whose faces cover less than a
   * quarter of its outline (the faces of the footprint, ground ones included,
   * that `allowed` gives it) or 0.5 m², but at most half of it, goes, unless
   * that would leave one of its faces, or of the planes gone before it in id
   * order, without an allowed plane that stays (such small planes included).
   * Their faces take the cheapest allowed plane among their neighbours', in
   * rounds so a face can follow a relabelled neighbour, and then the cheapest
   * allowed plane of all. Equal costs go to the lower plane id.
   *
   * @return the planes that went
   */
  std::set<int> drop_small_labels(
      Arrangement_2& arr, const IndexedPlanesWithPoints& roof_planes,
      const std::set<int>& narrow,
      std::function<bool(Arrangement_2::Face_handle, int)> allowed);

  std::unique_ptr<ArrangementOptimiserInterface> createArrangementOptimiser();
}  // namespace roofer::reconstruction
