#include "surfaceanalysis.h"
/*
 *  Copyright (c) 2026 Jeremy HU <jeremy-at-dust3d dot org>. All rights reserved.
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in all
 *  copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *  SOFTWARE.
 */
#include <AutoRemesher/AutoRemesher>
#include <AutoRemesher/IsotropicRemesher>
#include <AutoRemesher/MeshSeparator>
#include <AutoRemesher/Parameterizer>
#include <AutoRemesher/QuadExtractor>
#include <AutoRemesher/SurfaceMesh>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
// Qt defines `emit` as a macro, which collides with TBB profiling.h's `void emit()`.
// macOS `<mach/mach.h>` also defines `emit`. Undefine before including TBB headers.
#if defined(__APPLE__) || defined(emit)
#undef emit
#endif

// oneAPI TBB (2021+) moved headers under <oneapi/tbb/>. Use __has_include where
// available (clang + GCC) to pick the right path, falling back to the legacy path.
#if defined(__has_include)
#if __has_include(<oneapi/tbb/blocked_range.h>)
#include <oneapi/tbb/blocked_range.h>
#else
#include <tbb/blocked_range.h>
#endif
#if __has_include(<oneapi/tbb/mutex.h>)
#include <oneapi/tbb/mutex.h>
#else
#include <tbb/mutex.h>
#endif
#if __has_include(<oneapi/tbb/parallel_for.h>)
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_sort.h>
#else
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#endif
#else
#include <tbb/blocked_range.h>
#include <tbb/mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#endif
#include <cfloat>
#include <meshoptimizer.h>
#include <unordered_map>
#include <unordered_set>

namespace AutoRemesher {

namespace {
    const float parallelPhaseBegin = 0.03f;
    const float parallelPhaseEnd = 0.95f;

    // How an island's own 0..1 progress splits across its three stages, from the
    // measured cost of each on a typical model.  The phase report prints the
    // real accumulated times, so these can be re-checked against a run.
    const float islandResampleEnd = 0.17f;
    const float islandParameterizeEnd = 0.50f;
    // Quad extraction runs from islandParameterizeEnd to 1.0.

    const double decimateTriggerRatio = 8.0;
    const double decimateTargetRatio = 4.0;

    enum class PreparationRecovery { None,
        Source,
        Creases,
        Diagonals };
    struct LayoutQuality {
        double angle = 0, error = 0, area = 0, weight = 0, rimError = 0;
        size_t corners = 0, boundary = 0, nonmanifold = 0, parts = 0, nonquads = 0, collapsed = 0;
    };
    LayoutQuality measureLayout(QuadExtractor& mesh, const SurfaceAnalysis& reference, double sourceArea)
    {
        LayoutQuality q;
        const auto& points = mesh.remeshedVertices();
        const auto& faces = mesh.remeshedQuads();
        std::map<std::pair<size_t, size_t>, size_t> edges;
        for (const auto& f : faces) {
            Vector3 center;
            double area = 0;
            for (size_t v : f)
                center += points[v];
            center /= double(f.size());
            for (size_t k = 1; k + 1 < f.size(); ++k)
                area += Vector3::area(points[f[0]], points[f[k]], points[f[k + 1]]);
            q.area += area;
            q.error += area * reference.surfaceDistanceSquared(center);
            if (f.size() != 4) {
                bool trim = false;
                for (size_t k = 0; k < f.size(); ++k)
                    trim |= reference.onSourceBoundary((points[f[k]] + points[f[(k + 1) % f.size()]]) * .5);
                q.nonquads += !trim;
            }
            for (size_t k = 0; k < f.size(); ++k) {
                ++edges[std::minmax(f[k], f[(k + 1) % f.size()])];
                if (f.size() != 4)
                    continue;
                const Vector3 a = points[f[(k + 3) % 4]] - points[f[k]], b = points[f[(k + 1) % 4]] - points[f[k]];
                q.angle += area * (a.lengthSquared() * b.lengthSquared() > 0 ? std::asin(std::min(1.0, std::fabs(Vector3::dotProduct(a.normalized(), b.normalized())))) : M_PI / 2);
                ++q.corners;
                q.weight += area;
            }
        }
        // Authored holes and trimmed borders are geometry, not extraction defects.
        for (const auto& e : edges) {
            // Measure local rim displacement separately from holes inside the surface.
            if (e.second == 1 && reference.supportsRimConstraints()) {
                const Vector3 a = points[e.first.first], b = points[e.first.second];
                for (const Vector3 p : { a, b, (a + b) * .5 }) {
                    auto binding = reference.bindCurve(p, reference.length(), true);
                    if (binding.chain != SurfaceMesh::npos) {
                        binding.vertex = SurfaceMesh::npos;
                        q.rimError = std::max(q.rimError, (p - reference.projectCurve(binding, p)).lengthSquared());
                    }
                }
            }
            q.boundary += e.second == 1 && !reference.onSourceBoundary((points[e.first.first] + points[e.first.second]) * .5) && !reference.onSourceBoundary(points[e.first.first], points[e.first.second]);
            q.nonmanifold += e.second > 2;
            q.collapsed += (points[e.first.first] - points[e.first.second]).lengthSquared() <= std::pow(1e-10 * reference.length(), 2);
        }
        std::vector<std::vector<std::vector<size_t>>> parts;
        MeshSeparator::splitToIslands(faces, parts);
        // Numerical dust must not veto recovery of a complete source form.
        for (const auto& part : parts) {
            double area = 0;
            for (const auto& f : part)
                for (size_t k = 1; k + 1 < f.size(); ++k)
                    area += Vector3::area(points[f[0]], points[f[k]], points[f[k + 1]]);
            q.parts += area > 1e-6 * sourceArea;
        }
        q.nonquads -= std::min(q.nonquads, mesh.poleTriangles());
        q.angle = q.weight ? q.angle / q.weight : M_PI;
        q.error = q.area > 0 ? q.error / q.area : std::numeric_limits<double>::infinity();
        // A quad and its reversed copy enclose no surface; their angles
        // cannot veto recovery of the original component.
        if (faces.size() == 2 && faces[0].size() == 4 && faces[1].size() == 4) {
            const auto& a = faces[0];
            const auto& b = faces[1];
            for (size_t k = 0; k < 4; ++k)
                if (a[0] == b[k] && a[1] == b[(k + 3) % 4] && a[2] == b[(k + 2) % 4] && a[3] == b[(k + 1) % 4])
                    q.angle = M_PI;
        }
        q.error = .5 * (q.error + reference.missingSurfaceError(points, faces));
        q.area = q.area > 0 && sourceArea > 0 ? std::fabs(std::log(q.area / sourceArea)) : std::numeric_limits<double>::infinity();
        return q;
    }

    bool preferFeatureLayout(QuadExtractor& before, QuadExtractor& after,
        const SurfaceAnalysis& reference, double sourceArea, bool recoverConnectivity = false, PreparationRecovery preparation = PreparationRecovery::None)
    {
        const bool recoverPreparation = preparation == PreparationRecovery::Creases;
        const auto a = measureLayout(before, reference, sourceArea), b = measureLayout(after, reference, sourceArea);
        std::cerr << "Candidate quality (angle/error/area/boundary/parts/collapsed): "
                  << a.angle << '/' << a.error << '/' << a.area << '/' << a.boundary << '/' << a.parts << '/' << a.collapsed << " -> "
                  << b.angle << '/' << b.error << '/' << b.area << '/' << b.boundary << '/' << b.parts << '/' << b.collapsed << '\n';
        // Prefer jointly better boundary fitting and surface quality, even when
        // a different polygon subdivision changes the count of unmatched edges.
        const bool recoverRim = reference.supportsRimConstraints() && b.rimError < a.rimError && b.error <= a.error && b.angle <= a.angle && b.area <= a.area + .03 && b.parts <= a.parts;
        // Small curved forms can disappear without losing much total area.
        const bool recoverForm = a.area > .03 || (a.error > std::pow(.05 * reference.length(), 2) && b.error < a.error * .5 && b.boundary <= a.boundary);
        // A large area/angle recovery may slightly raise centroid-sampled error.
        // Bound this allowance to the diagonal proposal, with no new openings.
        const bool recoverArea = preparation == PreparationRecovery::Diagonals && b.parts == 1 && b.area < a.area * .5 && b.angle < a.angle && b.boundary <= a.boundary;
        const double fittingLimit = std::max(a.error * (recoverArea ? 1.05 : 1.), std::pow(.05 * reference.length(), 2));
        // Keep the same one-tenth-cell tolerance used to recognize authored rims.
        const bool keepsRim = b.rimError <= std::max(a.rimError, std::pow(.1 * reference.length(), 2));
        return b.corners && keepsRim && (!recoverPreparation || (b.parts == 1 && b.boundary <= a.boundary && b.error < a.error && b.area < a.area)) && b.nonmanifold <= a.nonmanifold && b.collapsed <= a.collapsed && (after.remeshedQuads().size() <= 2 * before.remeshedQuads().size() || (recoverForm && b.area < a.area * .5 && b.angle <= std::max(a.angle, .15))) && (recoverRim || (recoverForm && b.area < a.area && (b.area < a.area * .5 || b.error < a.error * .5 || (((recoverConnectivity && b.parts < a.parts) || (recoverPreparation && b.error < a.error)) && b.boundary <= a.boundary)) && b.angle <= a.angle + .15 && b.error < fittingLimit && b.parts <= std::max(size_t(1), a.parts)) ||
                   // Closing unintended holes permits small angle distortion, while
                   // retaining fitting and connectivity limits. Closing boundary edges
                   // may replace them with a small number of nonquad repair faces.
                   (a.boundary > 0 && b.boundary * 4 < a.boundary && b.angle < std::max(a.angle, .1) && b.area <= a.area + .03 && b.error < std::max(a.error, std::pow(.1 * reference.length(), 2)) && b.parts <= a.parts && b.nonquads + b.boundary <= a.nonquads + a.boundary) ||
                   // Avoid buying tiny fitting gains with a needlessly dense flat grid.
                   ((after.remeshedQuads().size() * 2 < before.remeshedQuads().size() || (after.remeshedQuads().size() * 4 < before.remeshedQuads().size() * 3 && b.angle < a.angle)) && b.nonquads == 0 && b.angle < .15 && b.area < .03 && b.error < std::pow(.05 * reference.length(), 2) && b.boundary <= a.boundary && b.parts <= a.parts) || (b.angle < a.angle && b.error <= a.error && b.area <= a.area + .03 && b.boundary <= a.boundary && b.nonmanifold <= a.nonmanifold && b.parts <= a.parts && b.nonquads <= a.nonquads));
    }

    bool isSmallConvexSource(const std::vector<Vector3>& vertices,
        const std::vector<std::vector<size_t>>& triangles)
    {
        // Bound the plane tests; dense or concave inputs retain full recovery.
        if (vertices.empty() || triangles.size() > 2048)
            return false;
        double extent = 0;
        for (const auto& p : vertices)
            extent = std::max(extent, (p - vertices.front()).length());
        const double tolerance = 1e-7 * extent;
        for (const auto& face : triangles) {
            const auto& origin = vertices[face[0]];
            const auto normal = Vector3::normal(origin, vertices[face[1]], vertices[face[2]]);
            if (normal.lengthSquared() == 0)
                return false;
            bool positive = false, negative = false;
            for (const auto& p : vertices) {
                const double distance = Vector3::dotProduct(p - origin, normal);
                positive |= distance > tolerance;
                negative |= distance < -tolerance;
                if (positive && negative)
                    return false;
            }
        }
        return true;
    }

    void markSharpEdgeVertices(const std::vector<Vector3>& vertices,
        const std::vector<unsigned int>& indices,
        double sharpEdgeRadians,
        std::vector<unsigned char>& vertexLock)
    {
        const size_t faceCount = indices.size() / 3;

        std::vector<Vector3> faceNormals(faceCount);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, faceCount),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    faceNormals[i] = Vector3::normal(vertices[indices[i * 3 + 0]],
                        vertices[indices[i * 3 + 1]],
                        vertices[indices[i * 3 + 2]]);
                }
            });

        std::vector<std::pair<uint64_t, unsigned int>> edges(faceCount * 3);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, faceCount),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    for (size_t j = 0; j < 3; ++j) {
                        unsigned int first = indices[i * 3 + j];
                        unsigned int second = indices[i * 3 + (j + 1) % 3];
                        if (first > second)
                            std::swap(first, second);
                        edges[i * 3 + j] = { ((uint64_t)first << 32) | second, (unsigned int)i };
                    }
                }
            });
        // Every (edge, face) pair is distinct, so the parallel sort produces the
        // same order the serial one did.
        tbb::parallel_sort(edges.begin(), edges.end());

        for (size_t i = 0; i + 1 < edges.size(); ++i) {
            if (edges[i].first != edges[i + 1].first)
                continue;
            if (Vector3::angle(faceNormals[edges[i].second],
                    faceNormals[edges[i + 1].second])
                < sharpEdgeRadians)
                continue;
            vertexLock[(size_t)(edges[i].first >> 32)] |= meshopt_SimplifyVertex_Priority;
            vertexLock[(size_t)(edges[i].first & 0xffffffffu)] |= meshopt_SimplifyVertex_Priority;
        }
    }
}

const double AutoRemesher::m_defaultSharpEdgeDegrees = 90;

double AutoRemesher::calculateAverageEdgeLength(const std::vector<Vector3>& vertices,
    const std::vector<std::vector<size_t>>& faces)
{
    double sumOfLength = 0.0;
    size_t edgeCount = 0;
    for (const auto& face : faces) {
        for (size_t i = 0; i < face.size(); ++i) {
            size_t j = (i + 1) % face.size();
            sumOfLength += (vertices[face[i]] - vertices[face[j]]).length();
            ++edgeCount;
        }
    }
    if (0 == edgeCount)
        return 0.0;
    return sumOfLength / edgeCount;
}

void AutoRemesher::initializeVoxelSize()
{
    double area = calculateMeshArea(m_vertices, m_triangles);
    double triangleArea = area / m_targetTriangleCount;
    m_voxelSize = std::sqrt(triangleArea / (0.86602540378 * 0.5));
#if AUTO_REMESHER_DEBUG
    std::cerr << "Area: " << area << " voxelSize: " << m_voxelSize << std::endl;
#endif
}

double AutoRemesher::calculateMeshArea(const std::vector<Vector3>& vertices,
    const std::vector<std::vector<size_t>>& triangles)
{
    double area = 0.0;
    for (const auto& it : triangles) {
        area += Vector3::area(vertices[it[0]], vertices[it[1]], vertices[it[2]]);
    }
    return area;
}

bool AutoRemesher::decimateIfTooDense(std::vector<Vector3>& vertices,
    std::vector<std::vector<size_t>>& triangles,
    double voxelSize,
    double sharpEdgeDegrees,
    size_t islandIndex,
    DecimationStats* stats)
{
    if (nullptr != stats)
        ++stats->islandsConsidered;

    if (vertices.empty() || triangles.empty() || voxelSize <= 0.0)
        return false;

    if (vertices.size() > (size_t)std::numeric_limits<unsigned int>::max())
        return false;

    const double targetTriangleArea = voxelSize * voxelSize * 0.86602540378 * 0.5;
    if (targetTriangleArea <= 0.0)
        return false;
    const double islandTargetTriangleCount = calculateMeshArea(vertices, triangles) / targetTriangleArea;
    if (islandTargetTriangleCount < 1.0)
        return false;

    if ((double)triangles.size() < islandTargetTriangleCount * decimateTriggerRatio)
        return false;

    const size_t decimateTriangleCount = (size_t)(islandTargetTriangleCount * decimateTargetRatio);

    std::vector<unsigned int> indices;
    indices.reserve(triangles.size() * 3);
    for (const auto& triangle : triangles) {
        if (3 != triangle.size())
            return false;
        for (size_t i = 0; i < 3; ++i)
            indices.push_back((unsigned int)triangle[i]);
    }

    Vector3 lowerBound = vertices.front();
    Vector3 upperBound = vertices.front();
    for (const auto& position : vertices) {
        for (size_t i = 0; i < 3; ++i) {
            lowerBound[i] = std::min(lowerBound[i], position[i]);
            upperBound[i] = std::max(upperBound[i], position[i]);
        }
    }
    const Vector3 center = (lowerBound + upperBound) * 0.5;

    std::vector<float> positions;
    positions.reserve(vertices.size() * 3);
    for (const auto& position : vertices) {
        positions.push_back((float)(position.x() - center.x()));
        positions.push_back((float)(position.y() - center.y()));
        positions.push_back((float)(position.z() - center.z()));
    }

    std::vector<unsigned int> remap(vertices.size());
    const size_t weldedVertexCount = meshopt_generateVertexRemap(remap.data(),
        indices.data(), indices.size(),
        positions.data(), vertices.size(), sizeof(float) * 3);
    std::vector<Vector3> weldedVertices(weldedVertexCount);
    std::vector<float> weldedPositions(weldedVertexCount * 3);
    meshopt_remapIndexBuffer(indices.data(), indices.data(), indices.size(), remap.data());
    meshopt_remapVertexBuffer(weldedPositions.data(), positions.data(),
        vertices.size(), sizeof(float) * 3, remap.data());
    for (size_t i = 0; i < vertices.size(); ++i) {
        if (~0u != remap[i])
            weldedVertices[remap[i]] = vertices[i];
    }

    std::vector<unsigned char> vertexLock;
    if (sharpEdgeDegrees > 0.0) {
        vertexLock.assign(weldedVertexCount, 0);
        markSharpEdgeVertices(weldedVertices, indices,
            sharpEdgeDegrees * (M_PI / 180.0), vertexLock);
    }

    std::vector<unsigned int> decimated(indices.size());
    float resultError = 0.0f;
    decimated.resize(meshopt_simplifyWithAttributes(decimated.data(),
        indices.data(), indices.size(),
        weldedPositions.data(), weldedVertexCount, sizeof(float) * 3,
        nullptr, 0, nullptr, 0,
        vertexLock.empty() ? nullptr : vertexLock.data(),
        decimateTriangleCount * 3, FLT_MAX, meshopt_SimplifyRegularize, &resultError));

    if (decimated.size() < 3 || decimated.size() >= indices.size())
        return false;

    std::vector<size_t> outputIndexOfWelded(weldedVertexCount, std::numeric_limits<size_t>::max());
    std::vector<Vector3> decimatedVertices;
    std::vector<std::vector<size_t>> decimatedTriangles;
    decimatedTriangles.reserve(decimated.size() / 3);
    for (size_t i = 0; i + 2 < decimated.size(); i += 3) {
        std::vector<size_t> triangle(3);
        for (size_t j = 0; j < 3; ++j) {
            const unsigned int weldedIndex = decimated[i + j];
            if (std::numeric_limits<size_t>::max() == outputIndexOfWelded[weldedIndex]) {
                outputIndexOfWelded[weldedIndex] = decimatedVertices.size();
                decimatedVertices.push_back(weldedVertices[weldedIndex]);
            }
            triangle[j] = outputIndexOfWelded[weldedIndex];
        }
        decimatedTriangles.push_back(triangle);
    }

    if (nullptr != stats) {
        ++stats->islandsDecimated;
        stats->trianglesBefore += triangles.size();
        stats->trianglesAfter += decimatedTriangles.size();
    }

#if AUTO_REMESHER_DEBUG
    std::cerr << "Island[" << islandIndex << "]: Decimated " << triangles.size()
              << " triangles to " << decimatedTriangles.size()
              << " (target " << decimateTriangleCount
              << ", island target " << (size_t)islandTargetTriangleCount
              << "), normalized error: " << resultError << std::endl;
#else
    (void)islandIndex;
    (void)resultError;
#endif

    vertices = std::move(decimatedVertices);
    triangles = std::move(decimatedTriangles);
    return true;
}

void AutoRemesher::resample(std::vector<Vector3>& vertices,
    std::vector<std::vector<size_t>>& triangles,
    double voxelSize,
    double adaptivity,
    double sharpEdgeDegrees,
    double smoothNormalDegrees,
    size_t islandIndex,
    DecimationStats* decimationStats,
    std::atomic<long long>* adaptiveFieldTimeUs,
    const ProgressHandler* progressHandler,
    std::vector<Vector3>* decimatedVerticesOut,
    std::vector<std::vector<size_t>>* decimatedTrianglesOut,
    const SurfaceAnalysis* analysis, bool balanceDiagonals, bool preserveCreases)
{
    auto t_decimateStart = std::chrono::high_resolution_clock::now();
    // Preserve preparation creases before curvature sampling (Remesher's 30-degree gate).
    if (!(analysis && analysis->featureLayout()) || vertices.size() > 50000)
        decimateIfTooDense(vertices, triangles, voxelSize,
            preserveCreases ? std::min(30., sharpEdgeDegrees) : sharpEdgeDegrees, islandIndex, decimationStats);
    if (nullptr != decimationStats) {
        decimationStats->timeUs += std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - t_decimateStart)
                                       .count();
    }

    if (nullptr != decimatedVerticesOut)
        *decimatedVerticesOut = vertices;
    if (nullptr != decimatedTrianglesOut)
        *decimatedTrianglesOut = triangles;

    auto t_fieldStart = std::chrono::high_resolution_clock::now();
    std::vector<double> vertexTargetLengths;
    if (adaptivity > 0.0 && analysis && !analysis->featureLayout()) {
        vertexTargetLengths.resize(vertices.size());
        tbb::parallel_for(size_t(0), vertices.size(), [&](size_t v) {
            vertexTargetLengths[v] = analysis->scalarSize(vertices[v]);
        });
        // Bound preparation work using the isotropic collapse spacing (4/5 h).
        // Remesher also separates its finite working budget from quad sizing.
        double samples = 0;
        for (const auto& t : triangles) {
            double density = 0;
            for (size_t v : t)
                density += 1 / std::pow(.8 * vertexTargetLengths[v], 2);
            samples += Vector3::area(vertices[t[0]], vertices[t[1]], vertices[t[2]]) * density / (3 * std::sqrt(3.) / 2);
            // Long, skinny input triangles also pay for splitting their edges.
            for (size_t k = 0; k < 3; ++k) {
                const size_t a = t[k], b = t[(k + 1) % 3];
                samples += .5 * std::max(0., std::ceil((vertices[a] - vertices[b]).length() / (4. / 3 * std::min(vertexTargetLengths[a], vertexTargetLengths[b]))) - 1);
            }
        }
        std::cerr << "Preparation estimate island " << islandIndex << " vertices " << samples << "\n";
        // Reserve the 25k feature trial from Remesher's 80k dense budget.
        if (samples > 80000 - 25000)
            vertexTargetLengths.clear();
    }
    if (nullptr != adaptiveFieldTimeUs) {
        *adaptiveFieldTimeUs += std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - t_fieldStart)
                                    .count();
    }

#if AUTO_REMESHER_DEBUG
    std::cerr << "Island[" << islandIndex << "]: Uniformly remeshing on target edge length: " << voxelSize << std::endl;
#endif
    IsotropicRemesher isotropicRemesher(vertices, triangles);
    if (nullptr != progressHandler && *progressHandler)
        isotropicRemesher.setProgressHandler(*progressHandler);
    isotropicRemesher.setRefineOnly(analysis && analysis->featureLayout());
    isotropicRemesher.setBalanceDiagonals(balanceDiagonals);
    isotropicRemesher.setTargetEdgeLength(voxelSize);
    if (!vertexTargetLengths.empty())
        isotropicRemesher.setVertexTargetEdgeLengths(&vertexTargetLengths);
    isotropicRemesher.setSharpEdgeDegrees(sharpEdgeDegrees);
    isotropicRemesher.setSmoothNormalDegrees(smoothNormalDegrees);
    isotropicRemesher.remesh();
    vertices = isotropicRemesher.remeshedVertices();
    triangles = isotropicRemesher.remeshedTriangles();
#if AUTO_REMESHER_DEBUG
    std::cerr << "Island[" << islandIndex << "]: Uniformly remesh done, vertex count: " << vertices.size() << " triangle count: " << triangles.size() << std::endl;
#endif
}

void AutoRemesher::updateProgress(size_t threadIndex, float progress, const char* status)
{
    if (nullptr == m_progressHandler)
        return;

    std::lock_guard<std::mutex> lock(m_progressMutex);
    if (threadIndex >= m_threadProgress.size())
        return;
    if (nullptr != status && '\0' != status[0])
        m_threadStatus[threadIndex] = status;
    if (progress > m_threadProgress[threadIndex]) {
        m_progressSum += (double)(progress - m_threadProgress[threadIndex])
            * m_threadProgressWeights[threadIndex];
        m_threadProgress[threadIndex] = progress;
    }

    const double overall = parallelPhaseBegin
        + (parallelPhaseEnd - parallelPhaseBegin) * std::min(1.0, std::max(0.0, m_progressSum));

    // Steps now report many times per island, so only wake the UI when the bar
    // would actually move or the status line would change.
    const int permille = (int)(overall * 1000.0);
    const char* islandStatus = m_threadStatus[threadIndex];
    if (permille == m_reportedPermille && islandStatus == m_reportedStatus)
        return;
    m_reportedPermille = permille;
    m_reportedStatus = islandStatus;

    // With several islands in flight, the run as a whole is only as far along as
    // its slowest island, so that is the step worth naming.
    size_t slowest = threadIndex;
    for (size_t i = 0; i < m_threadProgress.size(); ++i) {
        if (m_threadProgress[i] < m_threadProgress[slowest])
            slowest = i;
    }
    const char* name = m_threadStatus[slowest];
    m_progressHandler(m_tag, (float)overall, nullptr != name ? name : "");
}

ProgressHandler AutoRemesher::makeStageProgress(size_t islandIndex, float begin, float end, float stageOrder)
{
    return [this, islandIndex, begin, end, stageOrder,
               lastTime = std::chrono::high_resolution_clock::now(),
               lastName = (const char*)nullptr,
               lastOrder = 0.0f](float fraction, const char* name) mutable {
        const auto now = std::chrono::high_resolution_clock::now();
        if (nullptr != lastName) {
            accumulateStageTime(lastName, lastOrder,
                std::chrono::duration_cast<std::chrono::microseconds>(now - lastTime).count());
        }
        lastTime = now;
        lastName = name;
        lastOrder = stageOrder + fraction;
        updateProgress(islandIndex, begin + (end - begin) * fraction, name);
    };
}

void AutoRemesher::accumulateStageTime(const char* name, float order, long long microseconds)
{
    if (nullptr == name || '\0' == name[0])
        return;
    std::lock_guard<std::mutex> lock(m_stageTimingMutex);
    for (auto& it : m_stageTimes) {
        if (it.name == name) {
            it.microseconds += microseconds;
            return;
        }
    }
    m_stageTimes.push_back({ name, order, microseconds });
}

bool AutoRemesher::remesh()
{
    // Each run owns a fresh delivery, including when validation fails.
    m_preparedIslands.clear();
    m_remeshedVertices.clear();
    m_remeshedQuads.clear();
    m_decimatedVertices.clear();
    m_decimatedTriangles.clear();
    m_decimated = false;
    m_isotropicVertices.clear();
    m_isotropicTriangles.clear();
    m_isotropicTriangleUvs.clear();
    m_isotropicOriginalTriangleUvs.clear();
    m_isotropicSingularVertices.clear();
    m_isotropicExtractedConnections.clear();
    m_isotropicExtractedConnectionMoved.clear();
    m_phaseReport.clear();
    m_stageTimes.clear();
    m_reportedPermille = -1;
    m_reportedStatus = nullptr;
    // Validate inputs before any sizing math. In particular a zero target
    // triangle count would divide by zero in initializeVoxelSize().
    const char* invalidInputReason = nullptr;
    if (m_vertices.empty())
        invalidInputReason = "input mesh has no vertices";
    else if (m_triangles.empty())
        invalidInputReason = "input mesh has no triangles";
    else if (0 == m_targetTriangleCount)
        invalidInputReason = "target triangle count must be greater than zero";
    if (nullptr != invalidInputReason) {
        std::cerr << "Invalid remesh input: " << invalidInputReason << std::endl;
        if (nullptr != m_progressHandler)
            m_progressHandler(m_tag, 1.0, invalidInputReason);
        return false;
    }
    auto t_start = std::chrono::high_resolution_clock::now();

    // Each label names the step that is about to run, not the one that just
    // finished, so the status line matches what the process is actually doing.
    if (nullptr != m_progressHandler)
        m_progressHandler(m_tag, 0.0f, "Computing voxel size");
    auto t_voxelStart = std::chrono::high_resolution_clock::now();
    initializeVoxelSize();
    auto t_voxelEnd = std::chrono::high_resolution_clock::now();

    if (nullptr != m_progressHandler)
        m_progressHandler(m_tag, 0.01f, "Splitting mesh into islands");
    std::vector<std::vector<std::vector<size_t>>> trianglesIslands;
    std::vector<std::vector<size_t>> sourceTriangleIds;
    auto t_splitStart = std::chrono::high_resolution_clock::now();
    MeshSeparator::splitToIslands(m_triangles, trianglesIslands, &sourceTriangleIds);
    auto t_afterSplit = std::chrono::high_resolution_clock::now();

    if (trianglesIslands.empty()) {
        std::cerr << "Input mesh is empty" << std::endl;
        if (nullptr != m_progressHandler)
            m_progressHandler(m_tag, 1.0, "Input mesh is empty");
        return false;
    }

#if AUTO_REMESHER_DEBUG
    std::cerr << "Split to islands: " << trianglesIslands.size() << std::endl;
#endif

    struct IslandContext {
        std::shared_ptr<const SurfaceAnalysis> analysis;
        std::vector<Vector3> vertices;
        std::vector<std::vector<size_t>> triangles;
        double voxelSize;
        double samplingLength = 0;
        double targetQuads;
        double scaling;
        double adaptivity;
        double anisotropy;
        double sharpEdgeDegrees;
        double smoothNormalDegrees;
    };

    if (nullptr != m_progressHandler)
        m_progressHandler(m_tag, 0.02f, "Building island contexts");
    // Islands are compacted independently of each other, and writing into a
    // pre-sized vector by index keeps them in the original order.
    std::vector<IslandContext> islandContexes(trianglesIslands.size());
    m_preparedIslands.resize(trianglesIslands.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, trianglesIslands.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            for (size_t islandIndex = range.begin(); islandIndex != range.end(); ++islandIndex) {
                const auto& island = trianglesIslands[islandIndex];
                IslandContext& context = islandContexes[islandIndex];
                ReferenceSurface reference;
                reference.sourceTriangleIds = std::move(sourceTriangleIds[islandIndex]);
                context.triangles.reserve(island.size());
                std::unordered_map<size_t, size_t> oldToNewVertexMap;
                oldToNewVertexMap.reserve(island.size() * 2);
                for (const auto& face : island) {
                    std::vector<size_t> triangle;
                    triangle.reserve(3);
                    for (size_t i = 0; i < 3; ++i) {
                        auto insertResult = oldToNewVertexMap.insert({ face[i], context.vertices.size() });
                        if (insertResult.second) {
                            context.vertices.push_back(m_vertices[face[i]]);
                            reference.sourceVertexIds.push_back(face[i]);
                        }
                        triangle.push_back(insertResult.first->second);
                    }
                    context.triangles.push_back(std::move(triangle));
                }

                reference.vertices = context.vertices;
                reference.triangles = context.triangles;
                reference.sharpEdgeDegrees = m_sharpEdgeDegrees;
                const SurfaceMesh sourceMesh(reference.vertices, reference.triangles);
                reference.edgeFeatures.resize(sourceMesh.cornerCount(), ReferenceSurface::EdgeFeature::None);
                const double sharpRadians = m_sharpEdgeDegrees * M_PI / 180.0;
                for (size_t corner = 0; corner < sourceMesh.cornerCount(); ++corner) {
                    if (sourceMesh.isBoundaryCorner(corner))
                        reference.edgeFeatures[corner] = ReferenceSurface::EdgeFeature::Boundary;
                    else if (std::fabs(sourceMesh.normalAngle(corner)) > sharpRadians)
                        reference.edgeFeatures[corner] = ReferenceSurface::EdgeFeature::Sharp;
                }
                // No mutable owner remains once this snapshot is published.
                m_preparedIslands[islandIndex].reference = std::make_shared<const ReferenceSurface>(std::move(reference));

                context.analysis.reset(new SurfaceAnalysis(sourceMesh, m_voxelSize,
                    m_sharpEdgeDegrees, m_adaptivity, m_anisotropy));
                m_preparedIslands[islandIndex].analysis = context.analysis;
                // Zero means use the parameterizer's default scale of one.
                context.scaling = m_scaling > 0 ? m_scaling : 1.0;
                context.targetQuads = calculateMeshArea(context.vertices, context.triangles) / (.86602540378 * m_voxelSize * m_voxelSize * context.scaling * context.scaling);
                context.voxelSize = m_voxelSize;
                context.adaptivity = m_adaptivity;
                context.anisotropy = m_anisotropy;
                context.sharpEdgeDegrees = m_sharpEdgeDegrees;
                context.smoothNormalDegrees = m_smoothNormalDegrees;
            }
        });
    std::vector<IslandContext> featureContexts = islandContexes;
    auto t_buildEnd = std::chrono::high_resolution_clock::now();
    if (nullptr != m_progressHandler)
        m_progressHandler(m_tag, parallelPhaseBegin, "Remeshing uniformly");

    std::atomic<long long> resampleTime(0);
    std::atomic<long long> adaptiveFieldTime(0);
    DecimationStats decimationStats;

    {
        m_threadProgressWeights.assign(islandContexes.size(), 1.0f);
        for (size_t i = 0; i < islandContexes.size(); ++i) {
            if (!m_triangles.empty())
                m_threadProgressWeights[i] = (float)(((double)islandContexes[i].triangles.size() / m_triangles.size()));
        }
        m_threadProgress.assign(islandContexes.size(), 0.0f);
        m_threadStatus.assign(islandContexes.size(), nullptr);
        m_progressSum = 0.0;

        struct IsotropicPhase {
            IsotropicPhase(std::vector<IslandContext>* contexts,
                AutoRemesher* remesher,
                std::atomic<long long>* resampleTime,
                std::atomic<long long>* adaptiveFieldTime,
                DecimationStats* decimationStats,
                std::vector<std::vector<Vector3>>* islandVertices,
                std::vector<std::vector<std::vector<size_t>>>* islandTriangles,
                std::vector<std::vector<Vector3>>* decimatedIslandVertices,
                std::vector<std::vector<std::vector<size_t>>>* decimatedIslandTriangles)
                : m_contexts(contexts)
                , m_remesher(remesher)
                , m_resampleTime(resampleTime)
                , m_adaptiveFieldTime(adaptiveFieldTime)
                , m_decimationStats(decimationStats)
                , m_islandVertices(islandVertices)
                , m_islandTriangles(islandTriangles)
                , m_decimatedIslandVertices(decimatedIslandVertices)
                , m_decimatedIslandTriangles(decimatedIslandTriangles)
            {
            }

            void operator()(const tbb::blocked_range<size_t>& range) const
            {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    auto& ctx = (*m_contexts)[i];

                    m_remesher->updateProgress(i, 0.0f, "Remeshing uniformly");
                    const ProgressHandler isotropicProgress = m_remesher->makeStageProgress(i,
                        0.0f, islandResampleEnd, -1.0f);

                    auto t0 = std::chrono::high_resolution_clock::now();
                    resample(ctx.vertices, ctx.triangles, ctx.voxelSize, ctx.adaptivity, ctx.sharpEdgeDegrees, ctx.smoothNormalDegrees, i, m_decimationStats,
                        m_adaptiveFieldTime, &isotropicProgress,
                        &(*m_decimatedIslandVertices)[i], &(*m_decimatedIslandTriangles)[i], ctx.analysis.get());
                    auto t1 = std::chrono::high_resolution_clock::now();
                    *m_resampleTime += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

                    (*m_islandVertices)[i] = ctx.vertices;
                    (*m_islandTriangles)[i] = ctx.triangles;

                    m_remesher->updateProgress(i, islandResampleEnd);
                }
            }

        private:
            std::vector<IslandContext>* m_contexts = nullptr;
            AutoRemesher* m_remesher = nullptr;
            std::atomic<long long>* m_resampleTime = nullptr;
            std::atomic<long long>* m_adaptiveFieldTime = nullptr;
            DecimationStats* m_decimationStats = nullptr;
            std::vector<std::vector<Vector3>>* m_islandVertices = nullptr;
            std::vector<std::vector<std::vector<size_t>>>* m_islandTriangles = nullptr;
            std::vector<std::vector<Vector3>>* m_decimatedIslandVertices = nullptr;
            std::vector<std::vector<std::vector<size_t>>>* m_decimatedIslandTriangles = nullptr;
        };

        auto mergeIslands = [](const std::vector<std::vector<Vector3>>& islandVertices,
                                const std::vector<std::vector<std::vector<size_t>>>& islandTriangles,
                                std::vector<Vector3>& mergedVertices,
                                std::vector<std::vector<size_t>>& mergedTriangles) {
            for (size_t i = 0; i < islandVertices.size(); ++i) {
                const size_t vertexOffset = mergedVertices.size();
                mergedVertices.insert(mergedVertices.end(),
                    islandVertices[i].begin(), islandVertices[i].end());
                for (const auto& triangle : islandTriangles[i]) {
                    std::vector<size_t> offsetTriangle;
                    offsetTriangle.reserve(triangle.size());
                    for (const size_t index : triangle)
                        offsetTriangle.push_back(index + vertexOffset);
                    mergedTriangles.push_back(std::move(offsetTriangle));
                }
            }
        };

        m_isotropicVertices.clear();
        m_isotropicTriangles.clear();
        m_decimatedVertices.clear();
        m_decimatedTriangles.clear();
        std::vector<std::vector<Vector3>> isotropicIslandVertices(islandContexes.size());
        std::vector<std::vector<std::vector<size_t>>> isotropicIslandTriangles(islandContexes.size());
        std::vector<std::vector<Vector3>> decimatedIslandVertices(islandContexes.size());
        std::vector<std::vector<std::vector<size_t>>> decimatedIslandTriangles(islandContexes.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, islandContexes.size()),
            IsotropicPhase(&islandContexes, this, &resampleTime, &adaptiveFieldTime,
                &decimationStats,
                &isotropicIslandVertices, &isotropicIslandTriangles,
                &decimatedIslandVertices, &decimatedIslandTriangles));
        size_t vertexOffset = 0, triangleOffset = 0;
        for (size_t i = 0; i < m_preparedIslands.size(); ++i) {
            auto& prepared = m_preparedIslands[i];
            prepared.vertexOffset = vertexOffset;
            prepared.vertexCount = isotropicIslandVertices[i].size();
            prepared.triangleOffset = triangleOffset;
            prepared.triangleCount = isotropicIslandTriangles[i].size();
            vertexOffset += prepared.vertexCount;
            triangleOffset += prepared.triangleCount;
        }
        mergeIslands(isotropicIslandVertices, isotropicIslandTriangles,
            m_isotropicVertices, m_isotropicTriangles);
        m_decimated = decimationStats.islandsDecimated.load() > 0;
        if (m_decimated) {
            mergeIslands(decimatedIslandVertices, decimatedIslandTriangles,
                m_decimatedVertices, m_decimatedTriangles);
        }
    }
    // Remesher prepare/refinement_coordinator.cpp separates sampling density
    // from the output quota; even a small request needs a resolved work surface.
    const double requested = m_targetTriangleCount * .5;
    const double samplingCount = std::max(3000., requested < 8000 ? requested : requested < 15000 ? 8000 + (requested - 8000) * 4 / 7
                                                                                                  : 12000);
    const double samplingLength = std::sqrt(calculateMeshArea(m_vertices, m_triangles) / samplingCount);
    tbb::parallel_for(size_t(0), featureContexts.size(), [&](size_t i) {
        auto& c = featureContexts[i];
        // Remesher prepare/quota.cpp uses quad area, with a 20-cell spacing
        // floor for small pieces in large requests or multi-part meshes.
        c.voxelSize *= std::sqrt(std::sqrt(3.) / 2);
        if (m_targetTriangleCount > 800 || featureContexts.size() > 2)
            c.voxelSize = std::min(c.voxelSize, std::sqrt(calculateMeshArea(c.vertices, c.triangles) / 20));
        c.analysis.reset(new SurfaceAnalysis(SurfaceMesh(c.vertices, c.triangles), c.voxelSize,
            c.sharpEdgeDegrees, c.adaptivity, c.anisotropy, true, false));
        c.samplingLength = std::min(c.voxelSize, samplingLength);
        const auto start = std::chrono::high_resolution_clock::now();
        resample(c.vertices, c.triangles, c.samplingLength, c.adaptivity, c.sharpEdgeDegrees, c.smoothNormalDegrees,
            i, nullptr, &adaptiveFieldTime, nullptr, nullptr, nullptr, c.analysis.get());
        resampleTime += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count();
    });
    auto t_isotropicEnd = std::chrono::high_resolution_clock::now();

    class ParameterizationThread {
    public:
        size_t islandIndex = 0;
        size_t feedbackBefore = 0;
        bool feedbackRetried = false, feedbackAccepted = false;
        IslandContext* island = nullptr;
        std::shared_ptr<IslandContext> preparationProposal;
        std::unique_ptr<Parameterizer> parameterizer;
        std::unique_ptr<QuadExtractor> remesher, connectivityProposal;
        AutoRemesher* autoRemesher = nullptr;
        std::vector<std::vector<Vector2>> capturedUvs;
        std::vector<std::vector<Vector2>> capturedOriginalUvs;
        std::vector<uint8_t> capturedExtractedConnectionMoved;
        std::vector<Vector3> capturedSingularVertices;
        std::vector<size_t> capturedSingularVertexIndices;
        std::vector<std::pair<Vector3, Vector3>> capturedExtractedConnections;
    };

    std::vector<ParameterizationThread> parameterizationThreads(islandContexes.size());
    for (size_t i = 0; i < islandContexes.size(); ++i) {
        auto& thread = parameterizationThreads[i];
        auto& context = islandContexes[i];
        thread.islandIndex = i;
        thread.island = &context;
        thread.autoRemesher = this;
    }

    class SurfaceParameterizer {
    public:
        SurfaceParameterizer(std::vector<ParameterizationThread>* parameterizationThreads,
            std::atomic<long long>* parameterizeTime,
            std::atomic<long long>* extractTime, std::vector<IslandContext>* featureContexts)
            : m_parameterizationThreads(parameterizationThreads)
            , m_parameterizeTime(parameterizeTime)
            , m_extractTime(extractTime)
            , m_featureContexts(featureContexts)
        {
        }

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                auto& thread = (*m_parameterizationThreads)[i];

                bool layoutTrial = false, spacingRefinement = false;
                const auto run = [&](double scaling) {
                    const bool trial = layoutTrial;
                    auto t0 = std::chrono::high_resolution_clock::now();

                    const auto& vertices = thread.island->vertices;
                    const auto& triangles = thread.island->triangles;

                    if (vertices.empty() || triangles.empty()) {
                        // Still retire the island, otherwise its share of the bar
                        // is never filled in and the total stalls short of the end.
                        if (!trial)
                            thread.autoRemesher->updateProgress(thread.islandIndex, 1.0f);
                        return;
                    }

                    if (!trial)
                        thread.autoRemesher->updateProgress(thread.islandIndex, islandResampleEnd);
                    thread.parameterizer = std::make_unique<Parameterizer>(&vertices,
                        &triangles,
                        nullptr);
                    if (!trial)
                        thread.parameterizer->setProgressHandler(
                            thread.autoRemesher->makeStageProgress(thread.islandIndex,
                                islandResampleEnd, islandParameterizeEnd, 0.0f));
                    if (scaling > 0.0)
                        thread.parameterizer->setScaling(scaling);
                    thread.parameterizer->setSurfaceAnalysis(thread.island->analysis.get());
                    thread.parameterizer->setSpacingRefinement(spacingRefinement);
                    thread.parameterizer->setGradientAdaptivity(thread.island->adaptivity);
                    thread.parameterizer->setAnisotropy(thread.island->anisotropy);
                    thread.parameterizer->setSharpEdgeDegrees(thread.island->sharpEdgeDegrees);
                    bool parameterizeSucceeded = true;
                    try {
                        parameterizeSucceeded = thread.parameterizer->parameterize(layoutTrial);
                    } catch (const std::exception& e) {
                        // A pathological island must not abort the whole remesh,
                        // so log the parameterizer failure and skip its quads.
                        parameterizeSucceeded = false;
                        std::cerr << "Island " << (thread.islandIndex + 1)
                                  << ": parameterization failed (" << e.what()
                                  << "), skipping this island." << std::endl;
                    } catch (...) {
                        parameterizeSucceeded = false;
                        std::cerr << "Island " << (thread.islandIndex + 1)
                                  << ": parameterization failed (unknown error), skipping this island." << std::endl;
                    }

                    auto t1 = std::chrono::high_resolution_clock::now();
                    *m_parameterizeTime += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

                    if (parameterizeSucceeded) {
                        if (!trial)
                            thread.autoRemesher->updateProgress(thread.islandIndex, islandParameterizeEnd);
                        std::unique_ptr<std::vector<std::vector<Vector2>>> uvs = thread.parameterizer->takeTriangleUvs();
                        if (uvs) {
                            // Save a copy of UVs for the [param] preview overlay
                            thread.capturedUvs = *uvs;
                            thread.capturedOriginalUvs = thread.parameterizer->originalTriangleUvs();
                        }
                        // Capture singular vertex positions for the [param] preview
                        thread.capturedSingularVertices = thread.parameterizer->singularVertexPositions();
                        thread.capturedSingularVertexIndices = thread.parameterizer->singularVertexIndices();
                        thread.remesher = std::make_unique<QuadExtractor>(&vertices,
                            &triangles,
                            uvs.get());
                        thread.remesher->setSurfaceAnalysis(thread.island->analysis.get());
                        thread.remesher->setOriginalTriangleUvs(&thread.capturedOriginalUvs);
                        thread.remesher->setSingularVertices(&thread.capturedSingularVertexIndices);
                        thread.remesher->setFullTurnVertices(&thread.parameterizer->fullTurnVertices());
                        if (!trial)
                            thread.remesher->setProgressHandler(
                                thread.autoRemesher->makeStageProgress(thread.islandIndex,
                                    islandParameterizeEnd, 1.0f, 1.0f));
                        if (!thread.remesher->extract(false, bool(thread.preparationProposal))) {
                            thread.remesher.reset();
                        } else {
                            // Preserve the original source extraction. Only disconnected
                            // results need a source-supported proposal on the same UVs.
                            if (thread.island->analysis->featureLayout()) {
                                std::vector<std::vector<std::vector<size_t>>> parts;
                                MeshSeparator::splitToIslands(thread.remesher->remeshedQuads(), parts);
                                if (parts.size() > 1) {
                                    auto supported = std::make_unique<QuadExtractor>(&vertices, &triangles, uvs.get());
                                    supported->setSurfaceAnalysis(thread.island->analysis.get());
                                    supported->setOriginalTriangleUvs(&thread.capturedOriginalUvs);
                                    supported->setSingularVertices(&thread.capturedSingularVertexIndices);
                                    supported->setFullTurnVertices(&thread.parameterizer->fullTurnVertices());
                                    if (supported->extract(true, bool(thread.preparationProposal)))
                                        thread.connectivityProposal = std::move(supported);
                                }
                            }
                            thread.capturedExtractedConnections = thread.remesher->extractedConnections();
                            thread.capturedExtractedConnectionMoved = thread.remesher->extractedConnectionMoved();
                        }
                    }
                    if (!trial)
                        thread.autoRemesher->updateProgress(thread.islandIndex, 1.0f);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    *m_extractTime += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
                };
                run(thread.island->scaling);
                const auto feedback = [&]() {
                    if (!thread.remesher)
                        return;
                    const auto quadCount = [](const QuadExtractor& mesh) {
                        const auto& faces = mesh.remeshedQuads();
                        return std::count_if(faces.begin(), faces.end(), [](const std::vector<size_t>& f) { return f.size() == 4; });
                    };
                    const double target = thread.island->targetQuads;
                    const double before = quadCount(*thread.remesher);
                    if (target < 16 || before == 0 || std::fabs(std::log(before / target)) < .15)
                        return;
                    // Curvature imposes a physical resolution requirement; count
                    // feedback may refine it, but must not coarsen it away.
                    if (thread.island->adaptivity > 0 && before > target)
                        return;
                    // One bounded correction on the already prepared island. Keep
                    // the first delivery if the retry fails or worsens count/topology.
                    const auto defects = [](const QuadExtractor& mesh) {
                        std::map<std::pair<size_t, size_t>, size_t> edges;
                        for (const auto& f : mesh.remeshedQuads())
                            for (size_t k = 0; k < f.size(); ++k)
                                ++edges[std::minmax(f[k], f[(k + 1) % f.size()])];
                        std::vector<std::vector<std::vector<size_t>>> components;
                        MeshSeparator::splitToIslands(mesh.remeshedQuads(), components);
                        std::array<size_t, 3> result { { 0, 0, components.size() } };
                        for (const auto& e : edges) {
                            result[0] += e.second == 1;
                            result[1] += e.second > 2;
                        }
                        return result;
                    };
                    auto saved = std::move(thread);
                    thread = ParameterizationThread();
                    thread.island = saved.island;
                    thread.islandIndex = saved.islandIndex;
                    thread.autoRemesher = saved.autoRemesher;
                    run(thread.island->scaling * std::max(.8, std::min(1.25, std::sqrt(before / target))));
                    bool accepted = false;
                    if (thread.remesher && quadCount(*thread.remesher) > 0) {
                        const auto first = defects(*saved.remesher), second = defects(*thread.remesher);
                        accepted = std::fabs(std::log(quadCount(*thread.remesher) / target)) < std::fabs(std::log(before / target)) && second[0] <= first[0] && second[1] <= first[1] && second[2] == first[2];
                    }
                    if (!accepted)
                        thread = std::move(saved);
                    thread.feedbackBefore = before;
                    thread.feedbackRetried = true;
                    thread.feedbackAccepted = accepted;
                };
                feedback();
                // A small closed convex input can already have a regular, source-fitting
                // delivery. Keep that result instead of solving denser alternatives.
                const auto& source = *thread.autoRemesher->m_preparedIslands[i].reference;
                const bool smallClosedSource = m_parameterizationThreads->size() == 1
                    && isSmallConvexSource(source.vertices, source.triangles)
                    && std::find(source.edgeFeatures.begin(), source.edgeFeatures.end(), ReferenceSurface::EdgeFeature::Boundary) == source.edgeFeatures.end();
                const auto resolvedLayout = [&]() {
                    if (!smallClosedSource || !thread.remesher)
                        return false;
                    const auto& analysis = *thread.island->analysis;
                    const auto quality = measureLayout(*thread.remesher, analysis,
                        calculateMeshArea(source.vertices, source.triangles));
                    return quality.angle < .1 && quality.area < .03
                        && quality.error < std::pow(.1 * analysis.length(), 2)
                        && quality.boundary == 0 && quality.nonmanifold == 0
                        && quality.collapsed == 0 && quality.nonquads * 500 <= thread.remesher->remeshedQuads().size() && quality.parts == 1;
                };
                if (resolvedLayout()) {
                    std::cerr << "Resolved small source layout on island " << i << '\n';
                    continue;
                }
                layoutTrial = true;
                // Keep the baseline winner; one extra solve tests final-frame spacing
                // on the same prepared source under the unchanged quality checks.
                bool sourceLayoutAccepted = false;
                for (size_t variant = 0; variant < 5; ++variant) {
                    std::shared_ptr<IslandContext> crease;
                    if (variant == 3) {
                        // Preserve successful source layouts. Only failed source trials
                        // need a different preparation triangulation before the solve.
                        if (sourceLayoutAccepted)
                            continue;
                        auto& c = (*m_featureContexts)[i];
                        const auto& source = *thread.autoRemesher->m_preparedIslands[i].reference;
                        c.vertices = source.vertices;
                        c.triangles = source.triangles;
                        resample(c.vertices, c.triangles, c.samplingLength, c.adaptivity,
                            c.sharpEdgeDegrees, c.smoothNormalDegrees, i, nullptr, nullptr, nullptr,
                            nullptr, nullptr, c.analysis.get(), true);
                    }
                    if (variant == 4) {
                        const auto& source = *thread.autoRemesher->m_preparedIslands[i].reference;
                        if (source.vertices.size() <= 50000)
                            break;
                        // Keep an accepted source layout alive while testing preparation.
                        crease = std::make_shared<IslandContext>((*m_featureContexts)[i]);
                        crease->vertices = source.vertices;
                        crease->triangles = source.triangles;
                        resample(crease->vertices, crease->triangles, crease->samplingLength, crease->adaptivity,
                            crease->sharpEdgeDegrees, crease->smoothNormalDegrees, i, nullptr, nullptr, nullptr,
                            nullptr, nullptr, crease->analysis.get(), false, true);
                    }
                    spacingRefinement = variant == 2 || variant == 3;
                    auto saved = std::move(thread);
                    thread = ParameterizationThread();
                    thread.preparationProposal = std::move(crease);
                    thread.island = thread.preparationProposal ? thread.preparationProposal.get() : variant ? &(*m_featureContexts)[i]
                                                                                                            : saved.island;
                    thread.islandIndex = saved.islandIndex;
                    thread.autoRemesher = saved.autoRemesher;
                    run(thread.island->scaling);
                    const auto& reference = *thread.autoRemesher->m_preparedIslands[i].reference;
                    const double sourceArea = calculateMeshArea(reference.vertices, reference.triangles);
                    const auto recovery = variant == 4 ? PreparationRecovery::Creases : variant == 3 ? PreparationRecovery::Diagonals
                        : variant                                                                    ? PreparationRecovery::Source
                                                                                                     : PreparationRecovery::None;
                    bool accepted = thread.remesher && (!saved.remesher || preferFeatureLayout(*saved.remesher, *thread.remesher, *saved.island->analysis, sourceArea, false, recovery));
                    // Test both extractions against the retained winner: a local choice
                    // must not discard a source layout that beats the previous trial.
                    auto* winner = accepted ? thread.remesher.get() : saved.remesher.get();
                    if (thread.connectivityProposal && (!winner || preferFeatureLayout(*winner, *thread.connectivityProposal, *saved.island->analysis, sourceArea, true, recovery))) {
                        thread.remesher = std::move(thread.connectivityProposal);
                        accepted = true;
                        thread.capturedExtractedConnections = thread.remesher->extractedConnections();
                        thread.capturedExtractedConnectionMoved = thread.remesher->extractedConnectionMoved();
                    }
                    thread.connectivityProposal.reset();
                    std::cerr << "Feature layout island " << i << " variant " << variant << " accepted " << accepted << '\n';
                    sourceLayoutAccepted |= variant > 0 && accepted;
                    if (!accepted)
                        thread = std::move(saved);
                    else if (resolvedLayout())
                        break;
                }
            }
        }

    private:
        std::vector<ParameterizationThread>* m_parameterizationThreads = nullptr;
        std::vector<IslandContext>* m_featureContexts = nullptr;
        std::atomic<long long>* m_parameterizeTime = nullptr;
        std::atomic<long long>* m_extractTime = nullptr;
    };
    std::atomic<long long> parameterizeTimeAccumulated(0);
    std::atomic<long long> extractTimeAccumulated(0);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, parameterizationThreads.size()),
        SurfaceParameterizer(&parameterizationThreads,
            &parameterizeTimeAccumulated,
            &extractTimeAccumulated, &featureContexts));
    auto t_parallelEnd = std::chrono::high_resolution_clock::now();

    if (nullptr != m_progressHandler)
        m_progressHandler(m_tag, parallelPhaseEnd, "Merging mesh islands");

    m_isotropicVertices.clear();
    m_isotropicTriangles.clear();
    m_isotropicTriangleUvs.clear();
    m_isotropicOriginalTriangleUvs.clear();
    for (size_t i = 0; i < parameterizationThreads.size(); ++i) {
        const auto& thread = parameterizationThreads[i];
        const auto& c = *thread.island;
        auto& prepared = m_preparedIslands[i];
        prepared.analysis = c.analysis;
        prepared.vertexOffset = m_isotropicVertices.size();
        prepared.vertexCount = c.vertices.size();
        prepared.triangleOffset = m_isotropicTriangles.size();
        prepared.triangleCount = c.triangles.size();
        m_isotropicVertices.insert(m_isotropicVertices.end(), c.vertices.begin(), c.vertices.end());
        for (auto f : c.triangles) {
            for (auto& v : f)
                v += prepared.vertexOffset;
            m_isotropicTriangles.push_back(std::move(f));
        }
        m_isotropicTriangleUvs.resize(m_isotropicTriangles.size(), std::vector<Vector2>(3));
        m_isotropicOriginalTriangleUvs.resize(m_isotropicTriangles.size(), std::vector<Vector2>(3));
        if (thread.capturedUvs.size() == c.triangles.size())
            std::copy(thread.capturedUvs.begin(), thread.capturedUvs.end(), m_isotropicTriangleUvs.begin() + prepared.triangleOffset);
        if (thread.capturedOriginalUvs.size() == c.triangles.size())
            std::copy(thread.capturedOriginalUvs.begin(), thread.capturedOriginalUvs.end(), m_isotropicOriginalTriangleUvs.begin() + prepared.triangleOffset);
    }

    // Merge singular vertex positions from all islands (for [param] preview)
    m_isotropicSingularVertices.clear();
    for (size_t i = 0; i < parameterizationThreads.size(); ++i) {
        auto& thread = parameterizationThreads[i];
        if (thread.capturedSingularVertices.empty())
            continue;
        m_isotropicSingularVertices.insert(m_isotropicSingularVertices.end(),
            thread.capturedSingularVertices.begin(), thread.capturedSingularVertices.end());
    }

    // Merge the raw quad-extraction connections for the [param] preview.
    m_isotropicExtractedConnections.clear();
    m_isotropicExtractedConnectionMoved.clear();
    for (const auto& thread : parameterizationThreads) {
        m_isotropicExtractedConnections.insert(m_isotropicExtractedConnections.end(),
            thread.capturedExtractedConnections.begin(), thread.capturedExtractedConnections.end());
        m_isotropicExtractedConnectionMoved.resize(m_isotropicExtractedConnections.size()
                - thread.capturedExtractedConnections.size(),
            0);
        m_isotropicExtractedConnectionMoved.insert(m_isotropicExtractedConnectionMoved.end(),
            thread.capturedExtractedConnectionMoved.begin(),
            thread.capturedExtractedConnectionMoved.end());
        m_isotropicExtractedConnectionMoved.resize(m_isotropicExtractedConnections.size(), 0);
    }
    size_t curveVertices = 0;
    for (size_t i = 0; i < parameterizationThreads.size(); ++i) {
        auto& thread = parameterizationThreads[i];
        if (nullptr == thread.remesher)
            continue;
        if (thread.feedbackRetried) {
            const auto& faces = thread.remesher->remeshedQuads();
            const size_t count = std::count_if(faces.begin(), faces.end(), [](const std::vector<size_t>& f) { return f.size() == 4; });
            std::cerr << "Count feedback: " << thread.feedbackBefore << " -> " << count
                      << " target " << thread.island->targetQuads << " accepted " << thread.feedbackAccepted << '\n';
        }
        curveVertices += thread.remesher->constrainedCurveVertices();
        const auto& quads = thread.remesher->remeshedQuads();
        if (quads.empty())
            continue;
        const auto& vertices = thread.remesher->remeshedVertices();
        size_t vertexStartIndex = m_remeshedVertices.size();
        m_remeshedVertices.reserve(m_remeshedVertices.size() + vertices.size());
        for (const auto& it : vertices) {
            m_remeshedVertices.push_back(it);
        }
        for (const auto& it : quads) {
            std::vector<size_t> quad;
            quad.reserve(it.size());
            for (const auto& v : it)
                quad.push_back(vertexStartIndex + v);
            m_remeshedQuads.push_back(quad);
        }
    }

    std::cerr << "Source curve constrained vertices: " << curveVertices << '\n';
    auto t_mergeEnd = std::chrono::high_resolution_clock::now();

    const auto elapsedUs = [](const std::chrono::high_resolution_clock::time_point& from,
                               const std::chrono::high_resolution_clock::time_point& to) {
        return std::chrono::duration_cast<std::chrono::microseconds>(to - from).count();
    };
    const long long t_voxelUs = elapsedUs(t_voxelStart, t_voxelEnd);
    const long long t_splitUs = elapsedUs(t_splitStart, t_afterSplit);
    const long long t_buildUs = elapsedUs(t_afterSplit, t_buildEnd);
    const long long t_isotropicWallUs = elapsedUs(t_buildEnd, t_isotropicEnd);
    const long long t_parameterizeWallUs = elapsedUs(t_isotropicEnd, t_parallelEnd);
    const long long t_parallelWallUs = elapsedUs(t_buildEnd, t_parallelEnd);
    const long long t_mergeUs = elapsedUs(t_parallelEnd, t_mergeEnd);
    const long long t_totalUs = elapsedUs(t_start, t_mergeEnd);

    const long long t_decimateUs = decimationStats.timeUs.load();
    const long long t_adaptiveFieldUs = adaptiveFieldTime.load();
    const size_t decimatedIslands = decimationStats.islandsDecimated.load();

    m_phaseReport.clear();
    {
        std::ostringstream line;
        // Whole milliseconds hide the per-island steps on a mesh split into many
        // small islands, so keep one decimal place.
        const auto milliseconds = [](long long microseconds) {
            std::ostringstream value;
            value.setf(std::ios::fixed);
            value.precision(1);
            value << (double)microseconds / 1000.0 << " ms";
            return value.str();
        };
        auto phase = [&](const char* name, long long microseconds) {
            line.str(std::string());
            line << name << ": " << milliseconds(microseconds);
            m_phaseReport.push_back(line.str());
        };

        line.str(std::string());
        line << "Islands: " << islandContexes.size()
             << ", input triangles: " << m_triangles.size();
        m_phaseReport.push_back(line.str());

        phase("Compute voxel size", t_voxelUs);
        phase("Split into islands", t_splitUs);
        phase("Build island contexts", t_buildUs);

        line.str(std::string());
        if (decimatedIslands > 0) {
            line << "Mesh simplifier: RAN on " << decimatedIslands << " of "
                 << decimationStats.islandsConsidered.load() << " islands, "
                 << decimationStats.trianglesBefore.load() << " -> "
                 << decimationStats.trianglesAfter.load() << " triangles, "
                 << milliseconds(t_decimateUs);
        } else {
            line << "Mesh simplifier: SKIPPED (no island above "
                 << (long long)decimateTriggerRatio << "x target triangle count), "
                 << milliseconds(t_decimateUs);
        }
        m_phaseReport.push_back(line.str());

        // The accumulated figures sum the islands, so on a multi-island mesh they
        // add up to more than the wall clock next to them.  That gap is the point:
        // accumulated / wall is how many cores the phase actually kept busy.
        phase("Adaptive target length field (accumulated)", t_adaptiveFieldUs);
        phase("Isotropic remesh (accumulated)",
            resampleTime.load() - t_decimateUs - t_adaptiveFieldUs);
        phase("Parameterize (accumulated)", parameterizeTimeAccumulated.load());
        phase("Quad extract (accumulated)", extractTimeAccumulated.load());

        {
            std::lock_guard<std::mutex> lock(m_stageTimingMutex);
            std::sort(m_stageTimes.begin(), m_stageTimes.end(),
                [](const StageTime& first, const StageTime& second) {
                    return first.order < second.order;
                });
            for (const auto& it : m_stageTimes) {
                line.str(std::string());
                line << "    " << it.name << ": " << milliseconds(it.microseconds);
                m_phaseReport.push_back(line.str());
            }
        }

        phase("Isotropic phase wall clock", t_isotropicWallUs);
        phase("Parameterize phase wall clock", t_parameterizeWallUs);
        phase("Parallel phase wall clock", t_parallelWallUs);

        {
            const long long accumulated = resampleTime.load()
                + parameterizeTimeAccumulated.load() + extractTimeAccumulated.load();
            line.str(std::string());
            line.setf(std::ios::fixed);
            line.precision(2);
            line << "Cores kept busy across the parallel phase: "
                 << (t_parallelWallUs > 0 ? (double)accumulated / t_parallelWallUs : 0.0)
                 << " (islands are the unit of parallelism)";
            m_phaseReport.push_back(line.str());
        }

        phase("Merge islands", t_mergeUs);
        phase("Total", t_totalUs);
    }

    for (const auto& line : m_phaseReport)
        std::cerr << line << std::endl;

#if AUTO_REMESHER_DEBUG
    std::cerr << "Remesh done" << std::endl;
#endif

    if (nullptr != m_progressHandler)
        m_progressHandler(m_tag, 1.0, "Done");

    return true;
}

}
