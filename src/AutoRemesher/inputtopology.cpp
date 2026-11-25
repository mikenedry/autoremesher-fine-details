// SPDX-License-Identifier: MIT
#include "inputtopology.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_map>

namespace AutoRemesher {
namespace {
    using Edge = std::pair<size_t, size_t>;
    struct EdgeHash {
        size_t operator()(const Edge& edge) const
        {
            const size_t a = std::hash<size_t> {}(edge.first);
            return a ^ (std::hash<size_t> {}(edge.second) + size_t(0x9e3779b9) + (a << 6) + (a >> 2));
        }
    };
    struct Incidence {
        std::array<size_t, 3> faces {};
        size_t count = 0;
    };
    Edge edgeKey(size_t a, size_t b) { return { std::min(a, b), std::max(a, b) }; }
    bool finite(const Vector3& p)
    {
        return std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z());
    }
    bool forward(const std::vector<size_t>& f, size_t a, size_t b)
    {
        for (size_t i = 0; i < 3; ++i)
            if (f[i] == a && f[(i + 1) % 3] == b)
                return true;
        return false;
    }
}

std::vector<size_t> redundantInputFins(const std::vector<Vector3>& vertices,
    const std::vector<std::vector<size_t>>& triangles)
{
    std::unordered_map<Edge, Incidence, EdgeHash> edges;
    std::vector<size_t> uses(vertices.size(), 0);
    for (size_t i = 0; i < triangles.size(); ++i) {
        const auto& f = triangles[i];
        // This helper does not reinterpret malformed indices, polygons or NaNs.
        if (f.size() != 3 || f[0] == f[1] || f[1] == f[2] || f[2] == f[0])
            return {};
        for (size_t v : f) {
            if (v >= vertices.size() || !finite(vertices[v]))
                return {};
            ++uses[v];
        }
        for (size_t j = 0; j < 3; ++j) {
            auto& incidence = edges[edgeKey(f[j], f[(j + 1) % 3])];
            if (incidence.count < 3)
                incidence.faces[incidence.count] = i;
            ++incidence.count;
        }
    }
    std::vector<size_t> removed;
    for (size_t i = 0; i < triangles.size(); ++i) {
        const auto& f = triangles[i];
        size_t boundary = 0, attachment = 3;
        for (size_t j = 0; j < 3; ++j) {
            const size_t count = edges.at(edgeKey(f[j], f[(j + 1) % 3])).count;
            boundary += count == 1;
            if (count == 3)
                attachment = j;
        }
        if (boundary != 2 || attachment == 3)
            continue;
        const size_t a = f[attachment], b = f[(attachment + 1) % 3], tip = f[(attachment + 2) % 3];
        if (uses[tip] != 1)
            continue;
        std::array<size_t, 2> support {};
        size_t next = 0;
        for (size_t face : edges.at(edgeKey(a, b)).faces)
            if (face != i)
                support[next++] = face;
        if (forward(triangles[support[0]], a, b) == forward(triangles[support[1]], a, b))
            continue;

        const Vector3 edge = vertices[b] - vertices[a];
        const double scale = std::max({ std::fabs(edge.x()), std::fabs(edge.y()), std::fabs(edge.z()) });
        if (!(scale > 0) || !std::isfinite(scale))
            continue;
        const double length = (edge / scale).length();
        const Vector3 axis = (edge / scale) / length;
        const Vector3 offset = ((vertices[tip] - vertices[a]) / scale) / length;
        const double along = Vector3::dotProduct(offset, axis);
        if (!finite(offset) || !(along > 0 && along < 1))
            continue;
        const double height = Vector3::crossProduct(axis, offset).length();
        // Relative to the attachment edge, not model units or requested quads.
        if (!(height <= .005))
            continue;
        std::array<Vector3, 2> normals;
        std::array<Vector3, 2> supportOffsets;
        std::array<double, 2> supportHeight {};
        for (size_t j = 0; j < 2; ++j) {
            size_t third = a;
            for (size_t v : triangles[support[j]])
                if (v != a && v != b)
                    third = v;
            supportOffsets[j] = ((vertices[third] - vertices[a]) / scale) / length;
            normals[j] = Vector3::crossProduct(axis, supportOffsets[j]);
            supportHeight[j] = normals[j].length();
        }
        if (!(supportHeight[0] > 0 && supportHeight[1] > 0)
            || !std::isfinite(supportHeight[0]) || !std::isfinite(supportHeight[1]))
            continue;
        // Opposite sides of the edge, with supporting face normals within 5 deg.
        if (Vector3::dotProduct(normals[0] / supportHeight[0], normals[1] / supportHeight[1])
            > -0.9961946980917455)
            continue;
        if (height > .002 * std::min(supportHeight[0], supportHeight[1]))
            continue;
        if (height > 0) {
            // Only discard an inward fold. An equally small outward branch can
            // carry real detail. Compare half-spaces without depending on the
            // mesh's global winding; a flat sheet has no resolved inward side.
            bool inward = true;
            for (size_t j = 0; j < 2; ++j) {
                const Vector3 normal = normals[j] / supportHeight[j];
                const double bend = Vector3::dotProduct(supportOffsets[1 - j] / supportHeight[1 - j], normal);
                const double side = Vector3::dotProduct(offset, normal);
                if (!std::isfinite(bend) || !std::isfinite(side)
                    || std::fabs(bend) <= 64 * std::numeric_limits<double>::epsilon()
                    || (side != 0 && std::signbit(side) != std::signbit(bend)))
                    inward = false;
            }
            if (!inward)
                continue;
        }
        removed.push_back(i);
    }
    return removed;
}
}
