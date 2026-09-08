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
#include <AutoRemesher/ConstrainedLeastSquares>
#include <AutoRemesher/FrameField>
#include <AutoRemesher/Parameterizer>
#include <AutoRemesher/QuadParameterizer>
#include <AutoRemesher/SingularitySimplifier>
#include <AutoRemesher/SurfaceMesh>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>

namespace AutoRemesher {

namespace {

    std::vector<double> computeConformalScaling(const SurfaceMesh& mesh,
        const std::vector<int>& vertexCharges,
        const std::vector<double>& desiredFaceScaling,
        double fitting)
    {
        const size_t vertexCount = mesh.vertexCount(), faceCount = mesh.faceCount();
        std::vector<double> cotanWeight(mesh.cornerCount(), 0.0);
        std::vector<double> angleSum(vertexCount, 0.0), faceArea(faceCount, 0.0);
        std::vector<char> onBoundary(vertexCount, 0);
        for (size_t f = 0; f < faceCount; ++f) {
            for (size_t l = 0; l < 3; ++l) {
                const size_t c = 3 * f + l;
                const Vector3 first = mesh.edgeVector(c);
                const Vector3 second = -mesh.edgeVector(mesh.previousCorner(c));
                const double cross = Vector3::crossProduct(first, second).length();
                const double dot = Vector3::dotProduct(first, second);
                angleSum[mesh.cornerVertex(c)] += std::atan2(cross, dot);
                if (cross > 1e-20)
                    cotanWeight[mesh.nextCorner(c)] = .5 * dot / cross;
                faceArea[f] += cross / 6.0;
                if (mesh.oppositeCorner(c) == SurfaceMesh::npos) {
                    onBoundary[mesh.cornerVertex(c)] = 1;
                    onBoundary[mesh.cornerVertex(mesh.nextCorner(c))] = 1;
                }
            }
        }

        std::vector<double> desiredVertexScaling(vertexCount, 0.0), vertexWeight(vertexCount, 0.0);
        for (size_t f = 0; f < faceCount; ++f) {
            const double logScaling = std::log(std::max(1e-6, desiredFaceScaling[f]));
            for (size_t l = 0; l < 3; ++l) {
                const size_t v = mesh.cornerVertex(3 * f + l);
                desiredVertexScaling[v] += faceArea[f] * logScaling;
                vertexWeight[v] += faceArea[f];
            }
        }
        for (size_t v = 0; v < vertexCount; ++v)
            desiredVertexScaling[v] = vertexWeight[v] > 0 ? desiredVertexScaling[v] / vertexWeight[v] : 0.0;

        std::vector<double> deficit(vertexCount, 0.0);
        for (size_t v = 0; v < vertexCount; ++v) {
            const double gaussian = (onBoundary[v] ? M_PI : 2.0 * M_PI) - angleSum[v];
            const int charge = v < vertexCharges.size() ? vertexCharges[v] : 0;
            const double cone = .5 * M_PI * double(charge > 2 ? charge - 4 : charge);
            deficit[v] = gaussian - cone;
        }

        ConstrainedLeastSquares system(vertexCount);
        std::vector<std::pair<size_t, double>> row;
        for (size_t v = 0; v < vertexCount; ++v) {
            if (onBoundary[v] || vertexWeight[v] <= 0.0)
                continue;
            row.clear();
            double diagonal = 0.0;
            for (const size_t c : mesh.cornersAroundVertex(v)) {
                const size_t opposite = mesh.oppositeCorner(c);
                if (opposite == SurfaceMesh::npos)
                    continue;
                const double weight = cotanWeight[c] + cotanWeight[opposite];
                if (0.0 == weight)
                    continue;
                row.push_back({ mesh.cornerVertex(mesh.nextCorner(c)), -weight });
                diagonal += weight;
            }
            if (row.empty())
                continue;
            row.push_back({ v, diagonal });
            system.addEnergy(row, deficit[v], 1.0);
        }
        for (size_t v = 0; v < vertexCount; ++v)
            system.addEnergy({ { v, 1.0 } }, desiredVertexScaling[v], fitting);

        std::vector<double> logScaling;
        if (!system.solve(&logScaling) || logScaling.size() != vertexCount)
            return desiredFaceScaling;

        std::vector<double> result(faceCount, 1.0);
        for (size_t f = 0; f < faceCount; ++f) {
            double total = 0.0;
            for (size_t l = 0; l < 3; ++l)
                total += logScaling[mesh.cornerVertex(3 * f + l)];
            result[f] = total / 3.0;
        }
        double solvedDensity = 0.0, desiredDensity = 0.0;
        for (size_t f = 0; f < faceCount; ++f) {
            solvedDensity += faceArea[f] * std::exp(-2.0 * result[f]);
            desiredDensity += faceArea[f] / (desiredFaceScaling[f] * desiredFaceScaling[f]);
        }
        if (solvedDensity > 0.0 && desiredDensity > 0.0) {
            const double shift = .5 * std::log(solvedDensity / desiredDensity);
            for (double& value : result)
                value -= shift;
        }
        for (double& value : result)
            value = std::max(.2, std::min(5.0, std::exp(value)));
        return result;
    }

}

bool Parameterizer::parameterize(bool featureLayout)
{
#if AUTO_REMESHER_DEV
    {
        FILE* fp = fopen("debug-input-for-parameterization.obj", "wb");
        for (size_t i = 0; i < m_vertices->size(); ++i) {
            const auto& vertex = (*m_vertices)[i];
            fprintf(fp, "v %f %f %f\n", vertex[0], vertex[1], vertex[2]);
        }
        for (size_t i = 0; i < m_triangles->size(); ++i) {
            const auto& indices = (*m_triangles)[i];
            fprintf(fp, "f %zu %zu %zu\n", indices[0] + 1, indices[1] + 1, indices[2] + 1);
        }
        fclose(fp);
    }
#endif

    // The fractions below are the measured share of parameterization each step
    // costs; "Solving quad cover" dominates, so the bar has to keep moving
    // through it rather than sitting still until it finishes.
    const auto report = [this](float fraction, const char* name) {
        if (m_progressHandler)
            m_progressHandler(fraction, name);
    };

    report(0.0f, "Computing vertex normals");
    report(0.02f, "Building surface topology");
    // The parameterization pipeline uses the triangle/corner mesh;
    // no attribute-backed interchange mesh is constructed.
    SurfaceMesh topology(*m_vertices, *m_triangles);
    if (topology.faceCount() != m_triangles->size()) {
        std::cerr << "Topology rejected a non-triangle face" << std::endl;
        return false;
    }

    std::unique_ptr<SurfaceAnalysis> localAnalysis;
    if (!m_analysis)
        localAnalysis.reset(new SurfaceAnalysis(topology, topology.averageEdgeLength(),
            m_sharpEdgeDegrees, m_adaptivity, m_anisotropy));
    const SurfaceAnalysis& analysis = m_analysis ? *m_analysis : *localAnalysis;
    SurfaceGuidance guidance;
    if (analysis.featureLayout()) {
        SurfaceAnalysis measured(topology, analysis.length(), m_sharpEdgeDegrees, m_adaptivity, m_anisotropy,
            true, true, m_spacingRefinement);
        guidance = measured.transfer(topology, true);
        for (size_t f = 0; f < topology.faceCount(); ++f)
            if (guidance.featureCorners[3 * f] || guidance.featureCorners[3 * f + 1] || guidance.featureCorners[3 * f + 2])
                guidance.faces[f].confidence = 1;
    } else
        guidance = analysis.transfer(topology);
    // Convert physical sizes to the cover's working-edge units exactly once.
    const double relativeLength = analysis.length() / std::max(1e-12, topology.averageEdgeLength());
    const double adaptive = std::max(0., std::min(2., m_adaptivity));
    guidance.spacingUpper = relativeLength;
    guidance.spacingLower = relativeLength / (adaptive <= 1 ? 1 + 3 * adaptive : 4 + 8 * (adaptive - 1));
    guidance.spacingAspect = std::pow(2.3, std::max(0., std::min(1., m_anisotropy)));
    std::vector<double> faceScalingField(topology.faceCount());
    for (size_t f = 0; f < topology.faceCount(); ++f)
        faceScalingField[f] = guidance.faces[f].scale * (m_adaptivity > 0 ? relativeLength : std::min(1.0, relativeLength));

    report(0.03f, "Solving frame field");
    // Topology, field, and quad cover form the complete active path.
    std::vector<Vector3> field;
    if (nullptr != m_triangleFieldVectors) {
        field = *m_triangleFieldVectors;
    } else if (!FrameField::create(topology, m_sharpEdgeDegrees,
                   &field, &guidance, analysis.featureLayout())) {
        std::cerr << "Frame field solve failed" << std::endl;
        return false;
    }
    if (field.size() != topology.faceCount()) {
        std::cerr << "Frame field has the wrong face count" << std::endl;
        return false;
    }

    SingularitySimplifier simplifier(topology, &field);
    if (m_singularitySimplification) {
        report(0.17f, "Simplifying singularities");
        simplifier.setSharpEdgeDegrees(m_sharpEdgeDegrees);
        simplifier.setFeatureCorners(&guidance.featureCorners);
        simplifier.setMaximumPairDistance(m_maximumSingularityPairDistance);
        simplifier.simplify();
    }

    //faceScalingField = computeConformalScaling(topology, simplifier.vertexCharges(),
    //    faceScalingField, std::max(1e-4, .05 * m_adaptivity));

    std::vector<double> faceScalingU(m_triangles->size(), 1.0);
    std::vector<double> faceScalingV(m_triangles->size(), 1.0);
    for (size_t f = 0; f < topology.faceCount(); ++f) {
        const double alignment = Vector3::dotProduct(field[f].normalized(), guidance.faces[f].direction);
        const double ratio = std::pow(guidance.faces[f].ratio, 2 * alignment * alignment - 1);
        faceScalingU[f] = std::sqrt(ratio);
        faceScalingV[f] = 1 / faceScalingU[f];
    }
    // The cover solve is the longest single step here, so it reports its own
    // sub-steps from 0.28 onwards rather than going quiet until it finishes.
    ProgressHandler coverProgress;
    if (m_progressHandler) {
        coverProgress = [this](float fraction, const char* name) {
            m_progressHandler(0.28f + (0.99f - 0.28f) * fraction, name);
        };
    }
    QuadParameterizer::Result cover;
    if (!QuadParameterizer::parameterize(*m_vertices, *m_triangles,
            &field, m_scaling, m_sharpEdgeDegrees, &cover,
            &faceScalingField, &faceScalingU, &faceScalingV,
            coverProgress ? &coverProgress : nullptr, &guidance.featureCorners, featureLayout || analysis.featureLayout(),
            m_spacingRefinement && analysis.featureLayout() && adaptive > 0 ? &guidance : nullptr)) {
        std::cerr << "Quad cover solve failed" << std::endl;
        return false;
    }

    report(0.99f, "Collecting singularities");
    m_fullTurnVertices = cover.fullTurnVertices;
    m_originalTriangleUvs = cover.triangleUvs;
    m_triangleUvs = std::make_unique<std::vector<std::vector<Vector2>>>(cover.triangleUvs);
    m_singularVertexPositions.clear();
    m_singularVertexIndices.clear();
    for (const size_t v : cover.singularVertices) {
        if (v >= m_vertices->size())
            continue;
        m_singularVertexPositions.push_back((*m_vertices)[v]);
        m_singularVertexIndices.push_back(v);
    }
    if (simplifier.cancelledPairCount() > 0) {
        std::cerr << "Simplified cross field singularities: "
                  << simplifier.singularityCountBefore() << " -> "
                  << simplifier.singularityCountAfter() << " ("
                  << simplifier.cancelledPairCount() << " pair(s) cancelled)"
                  << std::endl;
    }
    report(1.0f, "");
    return true;
}

}
