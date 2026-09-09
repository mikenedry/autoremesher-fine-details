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
#include <AutoRemesher/MixedIntegerLeastSquares>
#include <AutoRemesher/QuadParameterizer>
#include <AutoRemesher/SurfaceAnalysis>
#include <AutoRemesher/SurfaceMesh>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace AutoRemesher {
namespace {
    Vector3 unit(const Vector3& v, const Vector3& fallback)
    {
        return v.length() < 1e-12 ? fallback.normalized() : v.normalized();
    }

    int edgeQuarterTurn(const SurfaceMesh& mesh, size_t corner,
        const std::vector<Vector3>& field, const std::vector<Vector3>& normals)
    {
        const size_t opposite = mesh.oppositeCorner(corner);
        if (opposite == SurfaceMesh::npos)
            return 0;
        const size_t f = mesh.cornerFace(corner), g = mesh.cornerFace(opposite);
        if (f > g) {
            const int r = edgeQuarterTurn(mesh, opposite, field, normals);
            return (4 - r) % 4;
        }
        size_t v0 = mesh.cornerVertex(corner), v1 = mesh.cornerVertex(mesh.nextCorner(corner));
        if (v1 < v0)
            std::swap(v0, v1);
        const Vector3 e = unit(mesh.position(v1) - mesh.position(v0), Vector3(1, 0, 0));
        const Vector3 y0 = unit(Vector3::crossProduct(normals[f], e), Vector3(0, 1, 0));
        const Vector3 yg = unit(Vector3::crossProduct(normals[g], e), Vector3(0, 1, 0));
        const double a0 = std::atan2(Vector3::dotProduct(field[f], y0), Vector3::dotProduct(field[f], e));
        int best = 0;
        double bestError = 1e100;
        Vector3 candidate = field[g];
        for (int r = 0; r < 4; ++r) {
            const double ag = std::atan2(Vector3::dotProduct(candidate, yg), Vector3::dotProduct(candidate, e));
            double d = std::fabs(a0 - ag);
            while (d > M_PI)
                d = std::fabs(d - 2.0 * M_PI);
            if (d < bestError) {
                bestError = d;
                best = r;
            }
            candidate = Vector3::crossProduct(normals[g], candidate);
        }
        return best;
    }

    enum EdgeConstraint { ConstraintNone = 0,
        ConstraintU = 1,
        ConstraintV = 2 };
    int edgeConstraint(const SurfaceMesh& mesh, size_t c, const std::vector<Vector3>& field,
        const std::vector<Vector3>& normals, double hardEdgeDegrees,
        const std::vector<char>* featureCorners, bool nearestAxis = false, bool preserveBoundary = true)
    {
        if (!(featureCorners && (*featureCorners)[c]) && mesh.oppositeCorner(c) != SurfaceMesh::npos && std::fabs(mesh.normalAngle(c)) * 180.0 / M_PI < hardEdgeDegrees)
            return ConstraintNone;
        const size_t f = mesh.cornerFace(c);
        const Vector3 edge = unit(mesh.edgeVector(c), Vector3(1, 0, 0));
        const Vector3 b = unit(field[f], edge);
        const Vector3 br = unit(Vector3::crossProduct(normals[f], b), Vector3(0, 1, 0));
        // Spacing transports only along a feature, even when the cover cannot
        // impose a hard axis constraint on its current field alignment.
        // Authored boundaries must remain isolines despite imperfect field alignment.
        const bool boundary = preserveBoundary && mesh.isBoundaryCorner(c);
        if (nearestAxis || boundary) {
            const double u = Vector3::dotProduct(edge, b), v = Vector3::dotProduct(edge, br);
            // Diagonal tangents need a coherent choice for perpendicular edges.
            if (boundary && std::fabs(std::fabs(u) - std::fabs(v)) <= 1e-12)
                return u * v > 0 ? ConstraintV : ConstraintU;
            return std::fabs(u) > std::fabs(v) ? ConstraintV : ConstraintU;
        }
        const bool alongB = std::acos(std::max(-1.0, std::min(1.0, std::fabs(Vector3::dotProduct(edge, b))))) < 10.0 * M_PI / 180.0;
        const bool alongBr = std::acos(std::max(-1.0, std::min(1.0, std::fabs(Vector3::dotProduct(edge, br))))) < 10.0 * M_PI / 180.0;
        if (alongB == alongBr)
            return ConstraintNone;
        return alongB ? ConstraintV : ConstraintU;
    }

    struct CoverContext {
        const SurfaceMesh& mesh;
        const std::vector<Vector3>& field;
        const std::vector<Vector3>& normals;
        const std::vector<int>& rotation;
        const std::vector<char>& seam;
        const std::vector<signed char>& cornerConstraints;
        const std::vector<double>& scalingU;
        const std::vector<double>& scalingV;
        const std::vector<double>* faceScaling;
        double scale;
        bool featureLayout;
        const std::vector<char>& fullTurns;
        double hardEdgeDegrees;
    };

    void initializeFieldAndNormals(const SurfaceMesh& mesh, const std::vector<Vector3>* guidance,
        std::vector<Vector3>* normals, std::vector<Vector3>* field)
    {
        const bool hasGuidance = guidance && guidance->size() == mesh.faceCount();
        normals->assign(mesh.faceCount(), Vector3());
        field->resize(mesh.faceCount());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, mesh.faceCount()), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t faceIndex = range.begin(); faceIndex != range.end(); ++faceIndex) {
                (*normals)[faceIndex] = unit(mesh.faceNormal(faceIndex), Vector3(0, 0, 1));
                Vector3 tangentAxis(1, 0, 0);
                if (std::fabs(Vector3::dotProduct(tangentAxis, (*normals)[faceIndex])) > .8)
                    tangentAxis = Vector3(0, 1, 0);
                Vector3 fieldDirection = hasGuidance ? (*guidance)[faceIndex] : tangentAxis;
                if (!hasGuidance)
                    fieldDirection = fieldDirection - (*normals)[faceIndex] * Vector3::dotProduct(fieldDirection, (*normals)[faceIndex]);
                (*field)[faceIndex] = unit(fieldDirection, mesh.edgeVector(3 * faceIndex));
            }
        });
    }

    void smoothCrossField(const SurfaceMesh& mesh, const std::vector<Vector3>& normals,
        double hardEdgeDegrees, std::vector<Vector3>* field)
    {
        std::vector<double> alpha(2 * mesh.faceCount(), 0.0);
        std::vector<char> locked(mesh.faceCount(), 0);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, mesh.faceCount()), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t faceIndex = range.begin(); faceIndex != range.end(); ++faceIndex) {
                alpha[2 * faceIndex] = 1.0;
                for (size_t localCorner = 0; localCorner < 3; ++localCorner) {
                    const size_t cornerIndex = 3 * faceIndex + localCorner, oppositeCornerIndex = mesh.oppositeCorner(cornerIndex);
                    if (oppositeCornerIndex != SurfaceMesh::npos && std::fabs(mesh.normalAngle(cornerIndex)) * 180.0 / M_PI < hardEdgeDegrees)
                        continue;
                    const Vector3 edge = unit(mesh.edgeVector(cornerIndex), Vector3(1, 0, 0));
                    const Vector3 fieldDirection = (*field)[faceIndex];
                    const Vector3 perpendicular = unit(Vector3::crossProduct(normals[faceIndex], fieldDirection), Vector3(0, 1, 0));
                    const double fieldAngle = std::atan2(Vector3::dotProduct(edge, perpendicular), Vector3::dotProduct(edge, fieldDirection));
                    alpha[2 * faceIndex] = std::cos(4.0 * fieldAngle);
                    alpha[2 * faceIndex + 1] = std::sin(4.0 * fieldAngle);
                    locked[faceIndex] = 1;
                }
            }
        });
        for (int iteration = 0; iteration < 40; ++iteration) {
            std::vector<double> next = alpha;
            tbb::parallel_for(tbb::blocked_range<size_t>(0, mesh.faceCount()), [&](const tbb::blocked_range<size_t>& range) {
                for (size_t f = range.begin(); f != range.end(); ++f)
                    if (!locked[f]) {
                        double x = alpha[2 * f], y = alpha[2 * f + 1];
                        const Vector3 bf = (*field)[f];
                        const Vector3 btf = unit(Vector3::crossProduct(normals[f], bf), Vector3(0, 1, 0));
                        for (size_t l = 0; l < 3; ++l) {
                            const size_t oc = mesh.oppositeCorner(3 * f + l);
                            if (oc == SurfaceMesh::npos)
                                continue;
                            const size_t g = mesh.cornerFace(oc);
                            Vector3 bg = (*field)[g];
                            bg = unit(bg - normals[f] * Vector3::dotProduct(bg, normals[f]), bf);
                            const double d = std::atan2(Vector3::dotProduct(bg, btf), Vector3::dotProduct(bg, bf));
                            const double cs = std::cos(4.0 * d), sn = std::sin(4.0 * d);
                            x += cs * alpha[2 * g] - sn * alpha[2 * g + 1];
                            y += sn * alpha[2 * g] + cs * alpha[2 * g + 1];
                        }
                        const double length = std::hypot(x, y);
                        if (length > 1e-12) {
                            next[2 * f] = x / length;
                            next[2 * f + 1] = y / length;
                        }
                    }
            });
            alpha.swap(next);
        }
        tbb::parallel_for(tbb::blocked_range<size_t>(0, mesh.faceCount()), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t faceIndex = range.begin(); faceIndex != range.end(); ++faceIndex) {
                const double fieldAngle = .25 * std::atan2(alpha[2 * faceIndex + 1], alpha[2 * faceIndex]);
                const Vector3 fieldDirection = (*field)[faceIndex];
                const Vector3 perpendicular = unit(Vector3::crossProduct(normals[faceIndex], fieldDirection), Vector3(0, 1, 0));
                (*field)[faceIndex] = unit(fieldDirection * std::cos(fieldAngle) + perpendicular * std::sin(fieldAngle), fieldDirection);
            }
        });
    }

    void brushFieldAlongSpanningTree(const SurfaceMesh& mesh, const std::vector<Vector3>& normals,
        std::vector<Vector3>* field)
    {
        std::vector<char> seen(mesh.faceCount(), 0);
        if (mesh.faceCount() != 0) {
            std::queue<size_t> q;
            q.push(0);
            seen[0] = 1;
            while (!q.empty()) {
                size_t f = q.front();
                q.pop();
                for (size_t l = 0; l < 3; ++l) {
                    size_t c = 3 * f + l, oc = mesh.oppositeCorner(c);
                    if (oc == SurfaceMesh::npos)
                        continue;
                    size_t g = mesh.cornerFace(oc);
                    if (!seen[g]) {
                        int turns = edgeQuarterTurn(mesh, c, *field, normals);
                        Vector3 brushed = (*field)[g];
                        for (int k = 0; k < turns; ++k)
                            brushed = Vector3::crossProduct(normals[g], brushed);
                        (*field)[g] = unit(brushed, mesh.edgeVector(3 * g));
                        seen[g] = 1;
                        q.push(g);
                    }
                }
            }
        }
    }

    std::vector<int> computeCornerRotations(const SurfaceMesh& mesh,
        const std::vector<Vector3>& field, const std::vector<Vector3>& normals)
    {
        const size_t corners = mesh.cornerCount();
        std::vector<int> rotation(corners, 0);
        for (size_t c = 0; c < corners; ++c) {
            const size_t oc = mesh.oppositeCorner(c);
            if (oc == SurfaceMesh::npos)
                continue;
            rotation[c] = edgeQuarterTurn(mesh, c, field, normals);
        }
        return rotation;
    }

    std::vector<signed char> computeCornerConstraints(const SurfaceMesh& mesh,
        const std::vector<Vector3>& field, const std::vector<Vector3>& normals, double hardEdgeDegrees,
        const std::vector<char>* featureCorners, bool nearestAxis = false, bool preserveBoundary = true)
    {
        const size_t corners = mesh.cornerCount();
        std::vector<signed char> cornerConstraints(corners, ConstraintNone);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, corners), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t c = range.begin(); c != range.end(); ++c)
                cornerConstraints[c] = static_cast<signed char>(edgeConstraint(mesh, c, field, normals, hardEdgeDegrees, featureCorners, nearestAxis, preserveBoundary));
        });
        return cornerConstraints;
    }

    void applyDirectionalSwaps(const SurfaceMesh& mesh, const std::vector<Vector3>& fieldBeforeBrush,
        const std::vector<Vector3>& field, const std::vector<Vector3>& normals,
        std::vector<double>* scalingU, std::vector<double>* scalingV)
    {
        size_t directionalSwaps = 0;
        for (size_t f = 0; f < mesh.faceCount(); ++f) {
            const Vector3 before = unit(fieldBeforeBrush[f], field[f]);
            const Vector3 after = unit(field[f], before);
            const Vector3 perpendicular = unit(Vector3::crossProduct(normals[f], after), before);
            if (std::fabs(Vector3::dotProduct(before, after)) < std::fabs(Vector3::dotProduct(before, perpendicular)))
                std::swap((*scalingU)[f], (*scalingV)[f]), ++directionalSwaps;
        }
    }

    const double quarterTurn[4][2][2] = {
        { { 1, 0 }, { 0, 1 } },
        { { 0, 1 }, { -1, 0 } },
        { { -1, 0 }, { 0, -1 } },
        { { 0, -1 }, { 1, 0 } }
    };

    void regularizeSpacing(const SurfaceMesh& mesh, const std::vector<int>& rotation,
        const std::vector<signed char>& constraints, const std::vector<double>& faceScale,
        const SurfaceGuidance& sizing, std::vector<double>& u, std::vector<double>& v)
    {
        // Remesher spacing.cpp/spacing_regularizer.cpp: final-frame stencil,
        // then one area-weighted solve that holds tight curvature targets.
        const size_t count = mesh.faceCount();
        std::vector<double> lengths(2 * count), areas(count), weights(2 * count);
        double area = 0;
        for (size_t f = 0; f < count; ++f) {
            lengths[2 * f] = faceScale[f] * u[f];
            lengths[2 * f + 1] = faceScale[f] * v[f];
            areas[f] = Vector3::crossProduct(mesh.edgeVector(3 * f), -mesh.edgeVector(3 * f + 2)).length();
            area += areas[f];
        }
        if (!(area > 0) || !(sizing.spacingLower > 0) || sizing.spacingUpper < sizing.spacingLower || std::any_of(lengths.begin(), lengths.end(), [](double x) { return !std::isfinite(x) || x <= 0; }))
            return;
        const auto allowed = [&](size_t c, size_t k) { return constraints[c] != (k ? ConstraintV : ConstraintU); };
        for (size_t pass = 0; pass < 11; ++pass) {
            auto next = lengths;
            const bool minimum = pass != 0 && pass != 3 && pass != 7 && pass != 10;
            for (size_t f = 0; f < count; ++f)
                for (size_t k = 0; k < 2; ++k) {
                    double sum = lengths[2 * f + k], weight = 1;
                    for (size_t c = 3 * f; c < 3 * f + 3; ++c) {
                        const size_t g = mesh.adjacentFace(c);
                        if (g == SurfaceMesh::npos || !allowed(c, k))
                            continue;
                        const double x = lengths[2 * g + (k ^ (rotation[c] & 1))];
                        if (minimum)
                            sum = std::min(sum, x);
                        else {
                            sum += x;
                            ++weight;
                        }
                    }
                    next[2 * f + k] = sum / weight;
                }
            lengths.swap(next);
        }
        for (size_t f = 0; f < count; ++f) {
            double& a = lengths[2 * f];
            double& b = lengths[2 * f + 1];
            a = std::min(a, b * sizing.spacingAspect);
            b = std::min(b, a * sizing.spacingAspect);
        }
        ConstrainedLeastSquares system(2 * count);
        for (size_t f = 0; f < count; ++f)
            for (size_t k = 0; k < 2; ++k) {
                const size_t i = 2 * f + k;
                const double t = sizing.spacingUpper > sizing.spacingLower ? std::max(0., std::min(1., (sizing.spacingUpper - lengths[i]) / (sizing.spacingUpper - sizing.spacingLower))) : 0;
                weights[i] = count * areas[f] / area * (sizing.faces[f].adaptiveWeight ? 10 * (1 + 100 * t) : 1);
                system.addEnergy({ { i, 1 } }, lengths[i], weights[i]);
            }
        for (size_t c = 0; c < mesh.cornerCount(); ++c) {
            const size_t o = mesh.oppositeCorner(c);
            if (o == SurfaceMesh::npos || o < c)
                continue;
            for (size_t k = 0; k < 2; ++k)
                if (allowed(c, k)) {
                    const size_t i = 2 * mesh.cornerFace(c) + k, j = 2 * mesh.cornerFace(o) + (k ^ (rotation[c] & 1));
                    system.addEnergy({ { i, 1 }, { j, -1 } }, 0, std::max(20., 3 * (weights[i] + weights[j])));
                }
        }
        std::vector<double> solution;
        if (!system.solve(&solution) || solution.size() != lengths.size() || std::any_of(solution.begin(), solution.end(), [](double x) { return !std::isfinite(x) || x <= 0; }))
            return;
        for (size_t f = 0; f < count; ++f) {
            u[f] = solution[2 * f] / faceScale[f];
            v[f] = solution[2 * f + 1] / faceScale[f];
        }
    }

    void applyCurlCorrection(const SurfaceMesh& mesh, const std::vector<Vector3>& normals,
        const std::vector<int>& rotation, const std::vector<signed char>& cornerConstraints,
        const std::vector<double>* faceScaling, double scale, double regularization,
        std::vector<double>* scalingU, std::vector<double>* scalingV,
        std::vector<Vector3>* field)
    {
        const size_t faceCount = mesh.faceCount();
        const bool hasFaceScaling = faceScaling && faceScaling->size() == faceCount;
        std::vector<double> su(faceCount), sv(faceCount);
        std::vector<Vector3> perpendicular(faceCount);
        for (size_t f = 0; f < faceCount; ++f) {
            const double faceScale = hasFaceScaling ? std::max(1e-12, (*faceScaling)[f]) : 1.0;
            su[f] = scale * faceScale * std::max(1e-12, (*scalingU)[f]);
            sv[f] = scale * faceScale * std::max(1e-12, (*scalingV)[f]);
            perpendicular[f] = unit(Vector3::crossProduct(normals[f], (*field)[f]),
                mesh.edgeVector(3 * f));
        }

        std::vector<char> anchored(faceCount, 0);
        for (size_t c = 0; c < mesh.cornerCount(); ++c)
            if (cornerConstraints[c] != ConstraintNone)
                anchored[mesh.cornerFace(c)] = 1;

        const auto request = [&](size_t corner, size_t face, double* constant,
                                 double* turnGradient, double* alongGradient, double* acrossGradient) {
            const Vector3 e = mesh.edgeVector(corner);
            const double along = Vector3::dotProduct((*field)[face], e);
            const double across = Vector3::dotProduct(perpendicular[face], e);
            constant[0] = along / su[face];
            constant[1] = across / sv[face];
            turnGradient[0] = across / su[face];
            turnGradient[1] = -along / sv[face];
            alongGradient[0] = -along / su[face];
            alongGradient[1] = 0.0;
            acrossGradient[0] = 0.0;
            acrossGradient[1] = -across / sv[face];
        };

        ConstrainedLeastSquares system(3 * faceCount);
        size_t edgeCount = 0;
        for (size_t c = 0; c < mesh.cornerCount(); ++c) {
            const size_t oc = mesh.oppositeCorner(c);
            if (oc == SurfaceMesh::npos || oc < c)
                continue;
            const size_t f = mesh.cornerFace(c), g = mesh.cornerFace(oc);
            const int r = (rotation[c] % 4 + 4) % 4;
            double here[2], hereTurn[2], hereAlong[2], hereAcross[2];
            double across[2], acrossTurn[2], acrossAlong[2], acrossAcross[2];
            request(c, f, here, hereTurn, hereAlong, hereAcross);
            request(oc, g, across, acrossTurn, acrossAlong, acrossAcross);
            const auto turned = [&](const double* value, size_t k) {
                return quarterTurn[r][k][0] * value[0] + quarterTurn[r][k][1] * value[1];
            };
            for (size_t k = 0; k < 2; ++k) {
                const double residual = here[k] + turned(across, k);
                system.addEnergy({ { 3 * f, hereTurn[k] },
                                     { 3 * f + 1, hereAlong[k] },
                                     { 3 * f + 2, hereAcross[k] },
                                     { 3 * g, turned(acrossTurn, k) },
                                     { 3 * g + 1, turned(acrossAlong, k) },
                                     { 3 * g + 2, turned(acrossAcross, k) } },
                    -residual, 1.0);
            }
            ++edgeCount;
        }
        if (0 == edgeCount)
            return;
        for (size_t f = 0; f < faceCount; ++f) {
            if (anchored[f])
                system.addConstraint({ { 3 * f, 1.0 } }, 0.0);
            else
                system.addEnergy({ { 3 * f, 1.0 } }, 0.0, regularization);
            system.addEnergy({ { 3 * f + 1, 1.0 } }, 0.0, regularization);
            system.addEnergy({ { 3 * f + 2, 1.0 } }, 0.0, regularization);
        }

        std::vector<double> correction;
        if (!system.solve(&correction) || correction.size() != 3 * faceCount)
            return;

        const double limit = .3, scaleLimit = std::log(1.5);
        for (size_t f = 0; f < faceCount; ++f) {
            const double angle = std::max(-limit, std::min(limit, correction[3 * f]));
            (*field)[f] = unit(std::cos(angle) * (*field)[f] + std::sin(angle) * perpendicular[f],
                (*field)[f]);
            (*scalingU)[f] *= std::exp(std::max(-scaleLimit, std::min(scaleLimit, correction[3 * f + 1])));
            (*scalingV)[f] *= std::exp(std::max(-scaleLimit, std::min(scaleLimit, correction[3 * f + 2])));
        }
        double areaBefore = 0.0, areaAfter = 0.0;
        for (size_t f = 0; f < faceCount; ++f) {
            const double faceScale = hasFaceScaling ? std::max(1e-12, (*faceScaling)[f]) : 1.0;
            const double area = Vector3::crossProduct(mesh.edgeVector(3 * f),
                -mesh.edgeVector(3 * f + 2))
                                    .length();
            areaBefore += area / (su[f] * sv[f]);
            su[f] = scale * faceScale * std::max(1e-12, (*scalingU)[f]);
            sv[f] = scale * faceScale * std::max(1e-12, (*scalingV)[f]);
            areaAfter += area / (su[f] * sv[f]);
        }
        if (areaBefore > 0.0 && areaAfter > 0.0) {
            const double factor = std::sqrt(areaAfter / areaBefore);
            for (size_t f = 0; f < faceCount; ++f) {
                (*scalingU)[f] *= factor;
                (*scalingV)[f] *= factor;
            }
        }
    }

    std::vector<char> computeSeam(const SurfaceMesh& mesh, const std::vector<int>& rotation, const std::vector<char>& fullTurns)
    {
        const size_t corners = mesh.cornerCount();
        std::vector<char> insideBall(corners, 0), ballSeen(mesh.faceCount(), 0);
        if (mesh.faceCount() != 0) {
            std::queue<size_t> q;
            q.push(0);
            ballSeen[0] = 1;
            while (!q.empty()) {
                const size_t f = q.front();
                q.pop();
                for (size_t l = 0; l < 3; ++l) {
                    const size_t c = 3 * f + l, oc = mesh.oppositeCorner(c);
                    if (oc == SurfaceMesh::npos || rotation[c] != 0)
                        continue;
                    const size_t g = mesh.cornerFace(oc);
                    if (!ballSeen[g]) {
                        ballSeen[g] = 1;
                        q.push(g);
                        insideBall[c] = insideBall[oc] = 1;
                    }
                }
            }
        }
        std::vector<char> seam(corners, 0);
        for (size_t c = 0; c < corners; ++c)
            seam[c] = (mesh.oppositeCorner(c) == SurfaceMesh::npos || !insideBall[c]);
        std::vector<size_t> borderDegree(mesh.vertexCount(), 0);
        for (size_t c = 0; c < corners; ++c)
            if (seam[c])
                ++borderDegree[mesh.cornerVertex(c)];
        bool changed = true;
        while (changed) {
            changed = false;
            for (size_t c = 0; c < corners; ++c) {
                const size_t oc = mesh.oppositeCorner(c);
                if (oc == SurfaceMesh::npos || !seam[c] || rotation[c] != 0)
                    continue;
                const size_t v0 = mesh.cornerVertex(c);
                if (borderDegree[v0] != 1 || fullTurns[v0])
                    continue;
                const size_t v1 = mesh.cornerVertex(mesh.nextCorner(c));
                seam[c] = seam[oc] = 0;
                insideBall[c] = insideBall[oc] = 1;
                if (borderDegree[v0] > 0)
                    --borderDegree[v0];
                if (borderDegree[v1] > 0)
                    --borderDegree[v1];
                changed = true;
            }
        }
        for (size_t c = 0; c < corners; ++c)
            if (fullTurns[mesh.cornerVertex(c)] || fullTurns[mesh.cornerVertex(mesh.nextCorner(c))])
                seam[c] = 1;
        return seam;
    }

    void addRotationConstraints(MixedIntegerLeastSquares& s, size_t ax, size_t bx, int r, double sign)
    {
        r = (r % 4 + 4) % 4;
        if (r == 0) {
            s.addConstraint(ax, 1, bx, sign);
            s.addConstraint(ax + 1, 1, bx + 1, sign);
        } else if (r == 1) {
            s.addConstraint(ax, 1, bx + 1, sign);
            s.addConstraint(ax + 1, 1, bx, -sign);
        } else if (r == 2) {
            s.addConstraint(ax, 1, bx, -sign);
            s.addConstraint(ax + 1, 1, bx + 1, -sign);
        } else {
            s.addConstraint(ax, 1, bx + 1, -sign);
            s.addConstraint(ax + 1, 1, bx, sign);
        }
    }

    void separateFeatureRegions(const CoverContext& ctx, MixedIntegerLeastSquares& system)
    {
        // Remesher param/speedup_regions.cpp + speedup_policy.cpp: connect
        // nearest feature-coordinate regions, retaining one cell across thin bands.
        // Charts stop at seams; their integer translations remain in the cover solve.
        const auto& mesh = ctx.mesh;
        const size_t count = mesh.cornerCount(), none = SurfaceMesh::npos;
        const auto value = [&](size_t i) { return double(float(system.value(i))); };
        struct Link {
            size_t to;
            double primary, orthogonal, direction;
        };
        struct Connection {
            size_t first, second;
            double primary, orthogonal, direction;
        };
        for (size_t axis = 0; axis < 2; ++axis) {
            std::vector<size_t> parent(count);
            std::iota(parent.begin(), parent.end(), 0);
            const auto root = [&](size_t x) {while (parent[x]!=x) {parent[x]=parent[parent[x]];x=parent[x];}return x; };
            const auto join = [&](size_t a, size_t b) { parent[root(a)] = root(b); };
            std::vector<std::vector<Link>> graph(count);
            std::vector<char> seed(count, 0);
            for (size_t c = 0; c < count; ++c) {
                const size_t n = mesh.nextCorner(c), o = mesh.oppositeCorner(c), f = mesh.cornerFace(c);
                if (o != none && !ctx.seam[c]) {
                    const size_t other = mesh.nextCorner(o);
                    join(c, other);
                    graph[c].push_back({ other, 0, 0, 0 });
                    graph[other].push_back({ c, 0, 0, 0 });
                }
                if (ctx.cornerConstraints[c] == (axis ? ConstraintV : ConstraintU)) {
                    join(c, n);
                    seed[c] = seed[n] = 1;
                }
                const double p = value(2 * n + axis) - value(2 * c + axis);
                const double q = value(2 * n + 1 - axis) - value(2 * c + 1 - axis);
                const Vector3 d = axis ? Vector3::crossProduct(ctx.normals[f], ctx.field[f]) : ctx.field[f];
                const double direction = Vector3::dotProduct(mesh.edgeVector(c), d);
                graph[c].push_back({ n, p, q, direction });
                graph[n].push_back({ c, -p, -q, -direction });
            }
            std::vector<size_t> owner(count, none), origin(count, none);
            std::vector<double> distance(count, 1e100), direction(count, 0);
            using Item = std::pair<double, size_t>;
            std::priority_queue<Item, std::vector<Item>, std::greater<Item>> queue;
            for (size_t c = 0; c < count; ++c)
                if (seed[c]) {
                    owner[c] = root(c);
                    origin[c] = c;
                    distance[c] = 0;
                    queue.push({ 0, c });
                }
            while (!queue.empty()) {
                const auto item = queue.top();
                queue.pop();
                const size_t c = item.second;
                if (item.first != distance[c])
                    continue;
                for (const auto& e : graph[c]) {
                    const double d = item.first + std::hypot(e.primary, e.orthogonal);
                    if (d < distance[e.to]) {
                        distance[e.to] = d;
                        owner[e.to] = owner[c];
                        origin[e.to] = origin[c];
                        direction[e.to] = direction[c] + e.direction;
                        queue.push({ d, e.to });
                    }
                }
            }
            std::map<std::pair<size_t, size_t>, Connection> connections;
            for (size_t c = 0; c < count; ++c)
                if (owner[c] != none)
                    for (const auto& e : graph[c])
                        if (owner[e.to] != none && owner[c] < owner[e.to]) {
                            const size_t a = origin[c], b = origin[e.to];
                            const double p = value(2 * b + axis) - value(2 * a + axis);
                            const double q = value(2 * b + 1 - axis) - value(2 * a + 1 - axis);
                            const auto key = std::make_pair(owner[c], owner[e.to]);
                            const auto old = connections.find(key);
                            if (old == connections.end() || std::fabs(p) + std::fabs(q) < std::fabs(old->second.primary) + std::fabs(old->second.orthogonal))
                                connections[key] = { a, b, p, q, direction[c] + e.direction - direction[e.to] };
                        }
            std::vector<Connection> ordered;
            for (const auto& item : connections)
                ordered.push_back(item.second);
            std::sort(ordered.begin(), ordered.end(), [](const Connection& a, const Connection& b) { return std::fabs(a.primary) < std::fabs(b.primary); });
            for (const auto& c : ordered)
                if (std::fabs(c.primary) < 1.35 && std::fabs(c.orthogonal) < 2 * std::fabs(c.primary) && std::fabs(c.direction) > 1e-12)
                    system.separateIntegerCoordinates(2 * c.first + axis, 2 * c.second + axis, c.direction > 0 ? 1 : -1);
        }
    }

    // Seamless integer-grid parameterization background: Bommes et al., MIQ (2009).
    // https://doi.org/10.1145/1531326.1531383
    // Directional weights and lattice-separation rules here are implementation-specific.
    bool solveQuadCover(const CoverContext& ctx, std::vector<double>* values,
        const ProgressHandler* progressHandler)
    {
        const SurfaceMesh& mesh = ctx.mesh;
        const std::vector<Vector3>& field = ctx.field;
        const std::vector<Vector3>& normals = ctx.normals;
        const std::vector<int>& rotation = ctx.rotation;
        const std::vector<char>& seam = ctx.seam;
        const std::vector<signed char>& cornerConstraints = ctx.cornerConstraints;
        const std::vector<double>& activeScalingU = ctx.scalingU;
        const std::vector<double>& activeScalingV = ctx.scalingV;
        const std::vector<double>* faceScaling = ctx.faceScaling;
        const double scale = ctx.scale;
        const size_t corners = mesh.cornerCount();
        const size_t uvVariables = 2 * corners;
        const size_t variables = 2 * uvVariables;
        const auto report = [progressHandler](float fraction, const char* name) {
            if (nullptr != progressHandler && *progressHandler)
                (*progressHandler)(fraction, name);
        };
        report(0.0f, "Building cover system");
        MixedIntegerLeastSquares s(variables);
        // Remesher uses unit chart offsets. Apply them to smooth open pieces;
        // retain the established lattice where interior hard-feature rows meet.
        bool open = false, interiorFeature = false;
        for (size_t c = 0; c < corners; ++c) {
            if (mesh.isBoundaryCorner(c))
                open = true;
            else
                interiorFeature |= std::fabs(mesh.normalAngle(c)) * 180 / M_PI >= ctx.hardEdgeDegrees;
        }
        for (size_t t = 0; t < 2 * corners; ++t)
            s.setVariablePeriod(uvVariables + t, ctx.featureLayout && open && !interiorFeature ? 1 : 2);
        std::vector<std::array<double, 4>> weights(mesh.faceCount(), { 1, 1, 1, 1 });
        const auto energies = [&]() {
            for (size_t f = 0; f < mesh.faceCount(); ++f) {
                const Vector3 u = field[f], v = unit(Vector3::crossProduct(normals[f], u), mesh.edgeVector(3 * f));
                const double faceScale = faceScaling && faceScaling->size() == mesh.faceCount()
                    ? std::max(1e-12, (*faceScaling)[f])
                    : 1.0;
                const double directionalU = std::max(1e-12, activeScalingU[f]);
                const double directionalV = std::max(1e-12, activeScalingV[f]);
                const double su = scale * faceScale * directionalU;
                const double sv = scale * faceScale * directionalV;
                // Remesher param/system.cpp: area-weighted derivatives in the
                // field frame, with 0.55 on the two nonzero-target rows.
                const auto& t = mesh.triangle(f);
                const double area2 = Vector3::crossProduct(mesh.position(t[1]) - mesh.position(t[0]),
                    mesh.position(t[2]) - mesh.position(t[0]))
                                         .length();
                if (area2 <= 1e-30)
                    continue;
                for (size_t coordinate = 0; coordinate < 2; ++coordinate)
                    for (size_t axis = 0; axis < 2; ++axis) {
                        std::vector<std::pair<size_t, double>> row;
                        for (size_t k = 0; k < 3; ++k) {
                            const Vector3 opposite = mesh.position(t[(k + 2) % 3]) - mesh.position(t[(k + 1) % 3]);
                            const Vector3 gradient = Vector3::crossProduct(normals[f], opposite) / area2;
                            row.push_back({ 2 * (3 * f + k) + coordinate, (coordinate ? sv : su) * Vector3::dotProduct(axis ? v : u, gradient) });
                        }
                        s.addEnergy(row, coordinate == axis ? 1 : 0, area2 * (coordinate == axis ? .55 : 1) * weights[f][2 * coordinate + axis]);
                    }
            }
        };
        energies();
        s.addConstraint(0, 1);
        s.addConstraint(1, 1);
        size_t hardCoordinateCount = 0;
        for (size_t c = 0; c < corners; ++c) {
            size_t oc = mesh.oppositeCorner(c);
            if (oc == SurfaceMesh::npos)
                continue;
            const size_t tc = uvVariables + 2 * c;
            const size_t toc = uvVariables + 2 * oc;
            int r = rotation[c];
            if (seam[c]) {
                addRotationConstraints(s, tc, toc, r, 1.0);
            } else {
                s.addConstraint(tc, 1);
                s.addConstraint(tc + 1, 1);
            }
        }
        for (size_t c = 0; c < corners; ++c) {
            const size_t oc = mesh.oppositeCorner(c);
            if (oc == SurfaceMesh::npos)
                continue;
            // Full-turn spokes meet at an unsplit far endpoint. Remesher's
            // singleton-chain rule leaves their pole copies independent.
            if (ctx.fullTurns[mesh.cornerVertex(c)])
                continue;
            const size_t other = mesh.nextCorner(oc);
            const size_t tc = uvVariables + 2 * c;
            const int r = (rotation[c] % 4 + 4) % 4;
            if (r == 0) {
                s.addConstraint(2 * c, 1, 2 * other, -1, tc, -1);
                s.addConstraint(2 * c + 1, 1, 2 * other + 1, -1, tc + 1, -1);
            } else if (r == 1) {
                s.addConstraint(2 * c, 1, 2 * other + 1, -1, tc, -1);
                s.addConstraint(2 * c + 1, 1, 2 * other, 1, tc + 1, -1);
            } else if (r == 2) {
                s.addConstraint(2 * c, 1, 2 * other, 1, tc, -1);
                s.addConstraint(2 * c + 1, 1, 2 * other + 1, 1, tc + 1, -1);
            } else {
                s.addConstraint(2 * c, 1, 2 * other + 1, 1, tc, -1);
                s.addConstraint(2 * c + 1, 1, 2 * other, -1, tc + 1, -1);
            }
        }
        for (size_t vertex = 0; vertex < mesh.vertexCount(); ++vertex) {
            const auto& incident = mesh.cornersAroundVertex(vertex);
            if (incident.empty() || ctx.fullTurns[vertex])
                continue;
            size_t start = incident.front(), c = start;
            int accumulated = 0;
            bool closed = true;
            std::vector<std::pair<size_t, int>> wheel;
            do {
                if (mesh.oppositeCorner(c) == SurfaceMesh::npos) {
                    closed = false;
                    break;
                }
                wheel.push_back({ c, accumulated });
                accumulated = (accumulated + rotation[c]) % 4;
                c = mesh.nextCorner(mesh.oppositeCorner(c));
            } while (c != start && wheel.size() <= incident.size() + 1);
            if (!closed || c != start || accumulated != 0)
                continue;
            for (int coord = 0; coord < 2; ++coord) {
                std::vector<std::pair<size_t, double>> row;
                for (const auto& item : wheel) {
                    const size_t t = uvVariables + 2 * item.first;
                    const int r = item.second;
                    if (coord == 0) {
                        if (r == 0)
                            row.push_back({ t, 1 });
                        else if (r == 1)
                            row.push_back({ t + 1, 1 });
                        else if (r == 2)
                            row.push_back({ t, -1 });
                        else
                            row.push_back({ t + 1, -1 });
                    } else {
                        if (r == 0)
                            row.push_back({ t + 1, 1 });
                        else if (r == 1)
                            row.push_back({ t, -1 });
                        else if (r == 2)
                            row.push_back({ t + 1, -1 });
                        else
                            row.push_back({ t, 1 });
                    }
                }
                s.addConstraint(row);
            }
        }
        for (size_t c = 0; c < corners; ++c) {
            const size_t n = mesh.nextCorner(c);
            const int constraint = cornerConstraints[c];
            if (constraint == ConstraintV) {
                s.setVariablePeriod(2 * c + 1, 1);
                s.setVariablePeriod(2 * n + 1, 1);
                hardCoordinateCount += 2;
                s.addConstraint(2 * c + 1, 1, 2 * n + 1, -1);
            } else if (constraint == ConstraintU) {
                s.setVariablePeriod(2 * c, 1);
                s.setVariablePeriod(2 * n, 1);
                hardCoordinateCount += 2;
                s.addConstraint(2 * c, 1, 2 * n, -1);
            }
        }
        report(0.2f, "Eliminating cover constraints");
        s.finalizeConstraints();
        if (ctx.featureLayout) {
            if (!s.solveIteration(false))
                return false;
            // Remesher speedup_reweight.cpp: spacing/flip penalties, transported
            // across ten adjacency passes, before discrete layout decisions.
            for (size_t f = 0; f < mesh.faceCount(); ++f) {
                const auto& t = mesh.triangle(f);
                const auto n = normals[f];
                const double area2 = Vector3::crossProduct(mesh.position(t[1]) - mesh.position(t[0]), mesh.position(t[2]) - mesh.position(t[0])).length();
                if (area2 <= 1e-30)
                    continue;
                // The original reweight stage reads stored float UVs. In
                // double precision a near-zero determinant can change sign.
                Vector2 uv[3];
                for (size_t k = 0; k < 3; ++k)
                    uv[k] = { float(s.value(2 * (3 * f + k))), float(s.value(2 * (3 * f + k) + 1)) };
                Vector3 gradient[2];
                for (size_t k = 0; k < 3; ++k)
                    for (size_t axis = 0; axis < 2; ++axis)
                        gradient[axis] += Vector3::crossProduct(n, mesh.position(t[(k + 2) % 3]) - mesh.position(t[(k + 1) % 3])) * (uv[k][axis] / area2);
                const auto a = uv[1] - uv[0], b = uv[2] - uv[0];
                const double flip = a.x() * b.y() - a.y() * b.x() < 0 ? 5 : 0;
                for (size_t axis = 0; axis < 2; ++axis) {
                    const double target = scale * (faceScaling ? (*faceScaling)[f] : 1) * (axis ? activeScalingV[f] : activeScalingU[f]);
                    const double stretch = 1 / std::max(1e-30, gradient[axis].length() * target);
                    weights[f][3 * axis] += std::min(3., std::max(0., (stretch - 1.8) * 3 / 2.2));
                }
                for (double& w : weights[f])
                    w = std::min(10., w + flip);
            }
            for (size_t pass = 0; pass < 10; ++pass) {
                auto next = weights;
                for (size_t f = 0; f < mesh.faceCount(); ++f)
                    for (size_t axis = 0; axis < 2; ++axis) {
                        double sum[2] = { weights[f][3 * axis], weights[f][1 + axis] }, count = 1;
                        for (size_t c = 3 * f; c < 3 * f + 3; ++c) {
                            const size_t g = mesh.adjacentFace(c);
                            if (g == SurfaceMesh::npos)
                                continue;
                            if (cornerConstraints[c] == (axis ? ConstraintV : ConstraintU))
                                continue;
                            const size_t other = axis ^ (rotation[c] & 1);
                            sum[0] += weights[g][3 * other];
                            sum[1] += weights[g][1 + other];
                            ++count;
                        }
                        next[f][3 * axis] = sum[0] / count;
                        next[f][1 + axis] = sum[1] / count;
                    }
                weights.swap(next);
            }
            s.clearEnergy();
            energies();
        }
        if (ctx.featureLayout) {
            if (!s.solveIteration(false))
                return false;
            separateFeatureRegions(ctx, s);
            if (!s.solveIteration(false))
                return false;
        }
        // A budgeted pass fixes at least one remaining integer variable. Allow
        // that worst case, plus the initial continuous solve, for every layout.
        const size_t maximumIterations = ctx.featureLayout ? s.integerKernelVariableCount() + 1 : 100;
        const size_t expectedIterations = 4;
        for (size_t iteration = 0; iteration < maximumIterations; ++iteration) {
            report(0.3f + 0.65f * std::min(1.0f, (float)iteration / expectedIterations),
                "Rounding cover to integers");
            if (!s.solveIteration(true, ctx.featureLayout))
                return false;
            if (s.converged())
                break;
        }
        report(0.98f, "Rounding cover to integers");
        values->resize(variables);
        for (size_t i = 0; i < variables; ++i)
            (*values)[i] = s.value(i);
        return s.converged();
    }

    void buildResultUv(const SurfaceMesh& mesh, const std::vector<int>& rotation,
        const std::vector<double>& allValues,
        size_t uvVariables, QuadParameterizer::Result* result)
    {
        std::vector<double> uv(allValues.begin(), allValues.begin() + uvVariables);
        for (double& coordinate : uv) {
            const double integer = std::round(coordinate);
            if (std::fabs(coordinate - integer) < 0.01)
                coordinate = integer;
        }
        result->singularVertices.clear();
        for (size_t vertex = 0; vertex < mesh.vertexCount(); ++vertex) {
            const auto& fan = mesh.cornersAroundVertex(vertex);
            int sum = 0;
            bool boundary = false;
            for (size_t c : fan) {
                sum = (sum + rotation[c]) % 4;
                if (mesh.oppositeCorner(c) == SurfaceMesh::npos)
                    boundary = true;
            }
            if (!boundary && sum != 0)
                result->singularVertices.push_back(vertex);
        }
        result->triangleUvs.assign(mesh.faceCount(), std::vector<Vector2>(3));
        for (size_t f = 0; f < mesh.faceCount(); ++f)
            for (size_t l = 0; l < 3; ++l) {
                size_t c = 3 * f + l;
                result->triangleUvs[f][l] = Vector2(uv[2 * c], uv[2 * c + 1]);
            }
    }

}

bool QuadParameterizer::parameterize(const std::vector<Vector3>& vertices,
    const std::vector<std::vector<size_t>>& triangles,
    const std::vector<Vector3>* guidance, double scaling,
    double hardEdgeDegrees, Result* result,
    const std::vector<double>* faceScaling,
    const std::vector<double>* faceScalingU,
    const std::vector<double>* faceScalingV,
    const ProgressHandler* progressHandler,
    const std::vector<char>* featureCorners, bool featureLayout, const SurfaceGuidance* sizing, bool preserveBoundary)
{
    const auto report = [progressHandler](float fraction, const char* name) {
        if (nullptr != progressHandler && *progressHandler)
            (*progressHandler)(fraction, name);
    };

    if ((faceScaling && faceScaling->size() != triangles.size()) || (featureCorners && featureCorners->size() != 3 * triangles.size()))
        return false;
    if (vertices.empty() || triangles.empty() || scaling <= 0.0)
        return false;
    report(0.0f, "Initializing cover field");
    SurfaceMesh mesh(vertices, triangles);
    if (mesh.faceCount() != triangles.size())
        return false;
    const size_t corners = mesh.cornerCount();
    const size_t uvVariables = 2 * corners;
    const double scale = std::max(1e-12, scaling * mesh.averageEdgeLength());

    std::vector<Vector3> normals;
    initializeFieldAndNormals(mesh, guidance, &normals, &result->field);
    const std::vector<Vector3> fieldBeforeBrush = result->field;

    std::vector<double> activeScalingU(mesh.faceCount(), 1.0), activeScalingV(mesh.faceCount(), 1.0);
    const bool trackDirectionalScale = faceScalingU && faceScalingV
        && faceScalingU->size() == mesh.faceCount() && faceScalingV->size() == mesh.faceCount();
    if (trackDirectionalScale) {
        activeScalingU = *faceScalingU;
        activeScalingV = *faceScalingV;
    }

    report(0.04f, "Smoothing cross field");
    if (!(guidance && guidance->size() == mesh.faceCount()))
        smoothCrossField(mesh, normals, hardEdgeDegrees, &result->field);
    brushFieldAlongSpanningTree(mesh, normals, &result->field);

    report(0.14f, "Computing corner rotations");
    const std::vector<int> rotation = computeCornerRotations(mesh, result->field, normals);
    result->cornerRotations = rotation;
    const std::vector<signed char> cornerConstraints = computeCornerConstraints(mesh, result->field,
        normals, hardEdgeDegrees, featureCorners, false, preserveBoundary);
    if (trackDirectionalScale)
        applyDirectionalSwaps(mesh, fieldBeforeBrush, result->field, normals,
            &activeScalingU, &activeScalingV);
    if (featureLayout && sizing && sizing->faces.size() == mesh.faceCount() && faceScaling && faceScaling->size() == mesh.faceCount())
        regularizeSpacing(mesh, rotation, computeCornerConstraints(mesh, result->field, normals, hardEdgeDegrees, featureCorners, true, preserveBoundary), *faceScaling, *sizing, activeScalingU, activeScalingV);
    report(0.20f, "Correcting field curl");

    if (!featureLayout)
        applyCurlCorrection(mesh, normals, rotation, cornerConstraints,
            faceScaling, scale, 1e-4, &activeScalingU, &activeScalingV, &result->field);
    // Turning-number/index background: Ray et al. (2008), Sections 2.5-2.6.
    // https://doi.org/10.1145/1356682.1356683 ; the full-turn tolerance below is local.
    // Remesher field/singularities.cpp: modulo-four charge misses radial
    // full-turn fans. Compare transported field holonomy with the corner sum.
    std::vector<char> fullTurns(mesh.vertexCount(), 0);
    if (featureLayout)
        for (size_t vertex = 0; vertex < mesh.vertexCount(); ++vertex) {
            double holonomy = 0, geometry = 0;
            bool boundary = false;
            for (size_t c : mesh.cornersAroundVertex(vertex)) {
                const size_t p = mesh.previousCorner(c), o = mesh.oppositeCorner(p), f = mesh.cornerFace(c);
                const Vector3 a = -mesh.edgeVector(p).normalized(), b = mesh.edgeVector(c).normalized();
                geometry += std::acos(std::max(-1., std::min(1., Vector3::dotProduct(a, b))));
                if (o == SurfaceMesh::npos) {
                    boundary = true;
                    continue;
                }
                const size_t g = mesh.cornerFace(o);
                const Vector3 e = mesh.edgeVector(p).normalized();
                const auto angle = [&](size_t h) { return std::atan2(Vector3::dotProduct(result->field[h], Vector3::crossProduct(normals[h], e)), Vector3::dotProduct(result->field[h], e)); };
                holonomy += std::remainder(angle(g) - angle(f), M_PI / 2);
            }
            fullTurns[vertex] = !boundary && std::fabs(holonomy - geometry) < .05;
        }
    result->fullTurnVertices.clear();
    for (size_t v = 0; v < fullTurns.size(); ++v)
        if (fullTurns[v])
            result->fullTurnVertices.push_back(v);
    const std::vector<char> seam = computeSeam(mesh, rotation, fullTurns);

    const CoverContext ctx { mesh, result->field, normals, rotation, seam, cornerConstraints,
        activeScalingU, activeScalingV, faceScaling, scale, featureLayout, fullTurns, hardEdgeDegrees };

    // The cover solve reports on its own 0..1, so remap it into the tail of this
    // function's range and keep the fractions monotonic end to end.
    ProgressHandler coverProgress;
    if (nullptr != progressHandler && *progressHandler) {
        coverProgress = [progressHandler](float fraction, const char* name) {
            (*progressHandler)(0.35f + 0.63f * fraction, name);
        };
    }
    std::vector<double> allValues;
    if (!solveQuadCover(ctx, &allValues, coverProgress ? &coverProgress : nullptr))
        return false;

    report(0.99f, "Building cover uvs");
    buildResultUv(mesh, rotation, allValues, uvVariables, result);
    return true;
}
}
