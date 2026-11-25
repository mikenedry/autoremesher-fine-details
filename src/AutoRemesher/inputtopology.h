// SPDX-License-Identifier: MIT
#ifndef AUTO_REMESHER_INPUT_TOPOLOGY_H
#define AUTO_REMESHER_INPUT_TOPOLOGY_H
#include <AutoRemesher/Vector3>
#include <vector>
namespace AutoRemesher {
// Original face IDs, in input order. Only isolated, negligible branches with
// an unambiguous continuous supporting sheet qualify; other input is untouched.
std::vector<size_t> redundantInputFins(const std::vector<Vector3>& vertices,
    const std::vector<std::vector<size_t>>& triangles);
}
#endif
