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
#include <AutoRemesher/FrameField>

#include <AutoRemesher/ConstrainedLeastSquares>
#include <AutoRemesher/SurfaceMesh>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace AutoRemesher {
namespace {
    // Fourfold field representation: Ray et al., N-Symmetry Direction Field Design (2008).
    // https://doi.org/10.1145/1356682.1356683
    constexpr double kSymmetry = 4.0;

    struct FacetTangentBasis {
        Vector3 tangent;
        Vector3 perpendicularTangent;
        Vector3 normal;
    };

    Vector3 normalizedOrFallback(const Vector3& vector, const Vector3& fallback)
    {
        return vector.length() <= 1e-12 ? fallback.normalized() : vector.normalized();
    }

    FacetTangentBasis createFacetTangentBasis(const SurfaceMesh& mesh, size_t face)
    {
        const Vector3 normal = normalizedOrFallback(mesh.faceNormal(face), Vector3(0, 0, 1));
        Vector3 tangent = mesh.edgeVector(3 * face);
        tangent = tangent - Vector3::dotProduct(tangent, normal) * normal;
        if (tangent.length() <= 1e-12) {
            tangent = std::fabs(normal.x()) < .9 ? Vector3(1, 0, 0) : Vector3(0, 1, 0);
            tangent = tangent - Vector3::dotProduct(tangent, normal) * normal;
        }
        tangent = normalizedOrFallback(tangent, Vector3(1, 0, 0));
        return { tangent, Vector3::crossProduct(normal, tangent), normal };
    }

    double tangentAngle(const Vector3& vector, const FacetTangentBasis& basis)
    {
        return std::atan2(Vector3::dotProduct(vector, basis.perpendicularTangent), Vector3::dotProduct(vector, basis.tangent));
    }

}

bool FrameField::create(const SurfaceMesh& mesh, double sharpEdgeDegrees,
    std::vector<Vector3>* field, const SurfaceGuidance* guidance, bool fixedCurvature)
{
    if (nullptr == field || mesh.faceCount() == 0)
        return false;
    const size_t faces = mesh.faceCount();
    if (guidance && (guidance->faces.size() != faces || guidance->featureCorners.size() != mesh.cornerCount()))
        return false;
    std::vector<FacetTangentBasis> facetBases(faces);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, faces), [&](const tbb::blocked_range<size_t>& range) {
        for (size_t faceIndex = range.begin(); faceIndex != range.end(); ++faceIndex)
            facetBases[faceIndex] = createFacetTangentBasis(mesh, faceIndex);
    });

    std::vector<double> periodic(2 * faces, 0.0), certainty(faces, 0.0);
    std::vector<char> locked(faces, 0);
    SurfaceGuidance local;
    if (!guidance) {
        SurfaceAnalysis analysis(mesh, mesh.averageEdgeLength(), sharpEdgeDegrees, 0, 0);
        local = analysis.transfer(mesh);
        guidance = &local;
    }
    for (size_t f = 0; f < faces; ++f) {
        const auto& face = guidance->faces[f];
        const double angle = kSymmetry * tangentAngle(face.direction, facetBases[f]);
        periodic[2 * f] = std::cos(angle);
        periodic[2 * f + 1] = std::sin(angle);
        certainty[f] = face.confidence;
        locked[f] = face.confidence >= 1.0 && (fixedCurvature || guidance->featureCorners[3 * f] || guidance->featureCorners[3 * f + 1] || guidance->featureCorners[3 * f + 2]);
    }

    const auto normalizePeriodic = [&]() {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, faces), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t faceIndex = range.begin(); faceIndex != range.end(); ++faceIndex) {
                const double periodicLength = std::hypot(periodic[2 * faceIndex], periodic[2 * faceIndex + 1]);
                if (periodicLength > 1e-30) {
                    periodic[2 * faceIndex] /= periodicLength;
                    periodic[2 * faceIndex + 1] /= periodicLength;
                }
            }
        });
    };

    normalizePeriodic();
    ConstrainedLeastSquares system(2 * faces);
    for (size_t f = 0; f < faces; ++f)
        if (locked[f]) {
            system.addConstraint({ { 2 * f, 1.0 } }, periodic[2 * f]);
            system.addConstraint({ { 2 * f + 1, 1.0 } }, periodic[2 * f + 1]);
        }
    for (size_t c = 0; c < mesh.cornerCount(); ++c) {
        const size_t other = mesh.oppositeCorner(c);
        if (other == SurfaceMesh::npos)
            continue;
        const size_t f = mesh.cornerFace(c), g = mesh.cornerFace(other);
        if (f < g)
            continue;
        const double transport = -kSymmetry * (tangentAngle(mesh.edgeVector(c), facetBases[g]) - tangentAngle(mesh.edgeVector(c), facetBases[f]));
        const double co = std::cos(transport), si = std::sin(transport);
        system.addEnergy({ { 2 * f, co }, { 2 * f + 1, si }, { 2 * g, -1.0 } }, 0.0);
        system.addEnergy({ { 2 * f, -si }, { 2 * f + 1, co }, { 2 * g + 1, -1.0 } }, 0.0);
    }
    std::vector<std::pair<size_t, size_t>> certaintyRows;
    certaintyRows.reserve(faces);
    for (size_t f = 0; f < faces; ++f)
        if (certainty[f] > 0.0) {
            const double weight = certainty[f] * certainty[f];
            const size_t rowU = system.addEnergy({ { 2 * f, 1.0 } }, periodic[2 * f], weight);
            const size_t rowV = system.addEnergy({ { 2 * f + 1, 1.0 } }, periodic[2 * f + 1], weight);
            certaintyRows.push_back({ rowU, rowV });
        }

    std::vector<double> solved;
    for (size_t iteration = 0; iteration < 5; ++iteration) {
        size_t rowIndex = 0;
        for (size_t f = 0; f < faces; ++f)
            if (certainty[f] > 0.0) {
                system.setEnergyRightHandSide(certaintyRows[rowIndex].first, periodic[2 * f]);
                system.setEnergyRightHandSide(certaintyRows[rowIndex].second, periodic[2 * f + 1]);
                ++rowIndex;
            }
        if (!system.solve(&solved))
            return false;
        periodic.swap(solved);
        normalizePeriodic();
    }
    field->resize(faces);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, faces), [&](const tbb::blocked_range<size_t>& range) {
        for (size_t faceIndex = range.begin(); faceIndex != range.end(); ++faceIndex) {
            const double fieldAngle = std::atan2(periodic[2 * faceIndex + 1], periodic[2 * faceIndex]) / kSymmetry;
            (*field)[faceIndex] = std::cos(fieldAngle) * facetBases[faceIndex].tangent + std::sin(fieldAngle) * facetBases[faceIndex].perpendicularTangent;
        }
    });
    return true;
}

}
