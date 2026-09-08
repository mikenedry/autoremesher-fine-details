#include "surfaceanalysis.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <numeric>
#include <queue>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <unordered_map>

namespace AutoRemesher {
namespace {
    double dot(const Vector3& a, const Vector3& b) { return Vector3::dotProduct(a, b); }
    Vector3 center(const SurfaceMesh& m, size_t f)
    {
        const auto& t = m.triangle(f);
        return (m.position(t[0]) + m.position(t[1]) + m.position(t[2])) / 3.0;
    }
    Vector3 segmentPoint(const Vector3& p, const Vector3& a, const Vector3& b)
    {
        const Vector3 d = b - a;
        return a + d * std::max(0.0, std::min(1.0, dot(p - a, d) / std::max(1e-30, d.lengthSquared())));
    }
    Vector3 trianglePoint(const SurfaceMesh& m, size_t f, const Vector3& p)
    {
        const auto& t = m.triangle(f);
        const Vector3 a = m.position(t[0]), b = m.position(t[1]), c = m.position(t[2]);
        const Vector3 n = Vector3::crossProduct(b - a, c - a);
        if (n.lengthSquared() > 1e-30) {
            const Vector3 q = p - n * (dot(p - a, n) / n.lengthSquared());
            if (dot(Vector3::crossProduct(b - a, q - a), n) >= 0 && dot(Vector3::crossProduct(c - b, q - b), n) >= 0 && dot(Vector3::crossProduct(a - c, q - c), n) >= 0)
                return q;
        }
        Vector3 q = segmentPoint(p, a, b);
        for (const auto& v : { segmentPoint(p, b, c), segmentPoint(p, c, a) })
            if ((p - v).lengthSquared() < (p - q).lengthSquared())
                q = v;
        return q;
    }
    double triangleDistance2(const SurfaceMesh& m, size_t f, const Vector3& p)
    {
        return (p - trianglePoint(m, f, p)).lengthSquared();
    }
    bool safeMove(const std::vector<Vector3>& vertices, const std::vector<std::vector<size_t>>& faces,
        const std::vector<size_t>& incident, size_t v, const Vector3& q)
    {
        for (size_t f : incident)
            for (size_t k = 1; k + 1 < faces[f].size(); ++k) {
                const size_t ids[3] = { faces[f][0], faces[f][k], faces[f][k + 1] };
                const Vector3 a = vertices[ids[0]], b = vertices[ids[1]], c = vertices[ids[2]];
                const Vector3 before = Vector3::crossProduct(b - a, c - a);
                const Vector3 x = ids[0] == v ? q : a, y = ids[1] == v ? q : b, z = ids[2] == v ? q : c;
                const Vector3 after = Vector3::crossProduct(y - x, z - x);
                if (after.lengthSquared() < .0025 * before.lengthSquared() || dot(before, after) <= .1 * before.length() * after.length())
                    return false;
            }
        return true;
    }

    constexpr size_t scaleCount = 5;

    size_t oppositeVertex(const SurfaceMesh& mesh, size_t corner, size_t vertex)
    {
        const size_t first = mesh.cornerVertex(corner);
        const size_t second = mesh.cornerVertex(mesh.nextCorner(corner));
        return first == vertex ? second : first;
    }

    void assignDirectionalMetric(SurfaceGuidance::Face& f, double u, double v, double aspect)
    {
        // Cap anisotropy by shortening the long axis, preserving the tight one.
        u = std::min(u, v * aspect);
        v = std::min(v, u * aspect);
        f.scale = std::sqrt(u * v);
        f.ratio = u / v;
    }

    std::array<double, scaleCount> curvatureSupportRadii(double referenceLength, double totalArea)
    {
        const double support = std::max(referenceLength, std::sqrt(totalArea / 10000.0));
        const double cap = totalArea > 0 ? .5 * std::sqrt(totalArea / 24) : support;
        double low = std::min(.33 * support, cap), high = std::min(1.75 * support, cap);
        std::array<double, scaleCount> radii;
        if (low == high)
            low = high * .33 / 1.75;
        // Remesher curvature_radii.cpp: cap the endpoints independently and use
        // the mean of the linear and geometric ladders.
        for (size_t s = 0; s < scaleCount; ++s) {
            const double t = double(s) / (scaleCount - 1);
            radii[s] = .5 * ((1 - t) * low + t * high + low * std::pow(high / low, t));
        }
        return radii;
    }

    SurfaceGuidance::Face evaluateCurvatureScales(const SurfaceMesh& mesh, size_t faceIndex,
        const std::array<Eigen::Matrix3d, scaleCount>& tensorsByRadius, const double* supportWeights,
        const double* supportRadii, bool useVertexTensor)
    {
        SurfaceGuidance::Face selected;
        const Vector3 u = mesh.edgeVector(3 * faceIndex).normalized(), v = Vector3::crossProduct(mesh.faceNormal(faceIndex), u);
        Eigen::Matrix<double, 3, 2> basis;
        basis << u.x(), v.x(), u.y(), v.y(), u.z(), v.z();
        SurfaceGuidance::Face samples[scaleCount];
        double normalAlignment[scaleCount] = {};
        for (size_t s = 0; s < scaleCount; ++s) {
            samples[s].radius = supportRadii[s];
            const Eigen::Matrix2d tensor = basis.transpose() * tensorsByRadius[s] * basis / (useVertexTensor ? std::max(1.0, supportWeights[s]) : std::max(1e-30, supportWeights[s]));
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(tensor);
            if (eig.info() != Eigen::Success)
                continue;
            const int major = std::fabs(eig.eigenvalues()[0]) > std::fabs(eig.eigenvalues()[1]) ? 0 : 1;
            samples[s].major = std::fabs(eig.eigenvalues()[major]);
            samples[s].minor = std::fabs(eig.eigenvalues()[1 - major]);
            const auto d = eig.eigenvectors().col(major);
            samples[s].direction = u * d[0] + v * d[1];
            if (useVertexTensor) {
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(tensorsByRadius[s] / std::max(1.0, supportWeights[s]));
                if (eigen.info() != Eigen::Success)
                    continue;
                std::array<int, 3> order { { 0, 1, 2 } };
                std::sort(order.begin(), order.end(), [&](int a, int b) { return std::fabs(eigen.eigenvalues()[a]) > std::fabs(eigen.eigenvalues()[b]); });
                samples[s].major = std::fabs(eigen.eigenvalues()[order[0]]);
                samples[s].minor = std::fabs(eigen.eigenvalues()[order[1]]);
                const auto axis = eigen.eigenvectors().col(order[0]);
                normalAlignment[s] = std::fabs(dot(mesh.faceNormal(faceIndex), Vector3(axis[0], axis[1], axis[2])));
                samples[s].direction = Vector3::crossProduct(mesh.faceNormal(faceIndex), Vector3(axis[0], axis[1], axis[2])).normalized();
            }
            samples[s].confidence = (samples[s].major - samples[s].minor) / std::max(1e-30, samples[s].major);
        }
        double bestScore = -1;
        for (size_t s = 0; s < scaleCount; ++s) {
            double agreement = 0, stability = 0, neighborCount = 0, directionDisagreement = 0;
            for (size_t j = 0; j < scaleCount; ++j)
                if (j != s && (j + 1 == s || s + 1 == j)) {
                    const double d = dot(samples[s].direction, samples[j].direction);
                    agreement += d * d;
                    directionDisagreement += std::acos(std::min(1.0, std::fabs(d)));
                    const double a = samples[s].major + samples[s].minor, b = samples[j].major + samples[j].minor;
                    stability += 1 - std::fabs(a - b) / std::max(1e-30, a + b);
                    ++neighborCount;
                }
            // Require enough spatial support and stable magnitudes as well as
            // direction agreement; isotropic curvature still contributes density.
            const double coverage = useVertexTensor ? (supportWeights[s] > 0 ? 1.0 : 0.0) : std::min(1.0, supportWeights[s] / (M_PI * supportRadii[s] * supportRadii[s]));
            const double confidence = samples[s].confidence * agreement / neighborCount * stability / neighborCount;
            // Remesher's strict selector takes precedence over the common
            // fallback; magnitude gating avoids directions from flat-face noise.
            bool strict = useVertexTensor;
            for (size_t j = s ? s - 1 : 0; j <= std::min(scaleCount - 1, s + 1); ++j)
                strict &= samples[j].confidence >= .7 && samples[j].major * supportRadii[j] >= 1 / 3.7979 && normalAlignment[j] <= .5;
            const double score = strict ? 3 - directionDisagreement / neighborCount : confidence * std::sqrt(coverage) + .1 * stability / neighborCount + .002 * s;
            if (score > bestScore) {
                bestScore = score;
                selected = samples[s];
                selected.confidence = confidence;
                selected.adaptiveWeight = strict || (samples[s].confidence > .5 && samples[s].major * supportRadii[s] < 1 / (1.5 * 3.7979));
                // Common attachment is capped at 60/500; strict winners are fixed.
                if (useVertexTensor)
                    selected.confidence = strict ? 1 : std::min(.12, confidence) * std::min(1., samples[s].major * supportRadii[s] * 3.7979);
            }
        }
        return selected;
    }

}

SurfaceAnalysis::SurfaceAnalysis(const SurfaceMesh& mesh, double length,
    double sharpDegrees, double adaptivity, double anisotropy, bool featureLayout, bool measureCurvature, bool deferSpacing)
    : m_mesh(mesh)
    , m_length(std::max(1e-12, length))
    , m_features(mesh.cornerCount(), 0)
    , m_faces(mesh.faceCount())
    , m_featureLayout(featureLayout)
{
    traceFeatureChains(sharpDegrees);
    resolveFeatureJunctions();
    if (measureCurvature) {
        if (m_featureLayout)
            measureVertexCurvature();
        else
            measureFaceCurvature();
        const double aspect = initializeDirectionalSizes(adaptivity, anisotropy);
        regularizeDirectionalSizes(aspect, deferSpacing);
    }
    buildSourceIndex();
}

void SurfaceAnalysis::traceFeatureChains(double sharpDegrees)
{
    const auto& mesh = m_mesh;
    const size_t none = SurfaceMesh::npos;
    // Decompose angle candidates at junctions and turns; preserve boundaries and
    // the user's hard-angle edges, reject short/incoherent automatic fragments.
    std::vector<std::vector<size_t>> incident(mesh.vertexCount());
    std::vector<size_t> edges;
    bool open = false;
    const double hard = sharpDegrees * M_PI / 180.0, candidate = std::min(hard, M_PI / 6);
    for (size_t c = 0; c < mesh.cornerCount(); ++c) {
        const size_t o = mesh.oppositeCorner(c);
        open |= o == none;
        if (o != none && (o < c || std::fabs(mesh.normalAngle(c)) <= candidate))
            continue;
        edges.push_back(c);
        incident[mesh.cornerVertex(c)].push_back(c);
        incident[mesh.cornerVertex(mesh.nextCorner(c))].push_back(c);
    }
    const auto joint = [&](size_t v) {
        const auto& e = incident[v];
        return e.size() != 2 || dot((mesh.position(oppositeVertex(mesh, e[0], v)) - mesh.position(v)).normalized(), (mesh.position(oppositeVertex(mesh, e[1], v)) - mesh.position(v)).normalized()) > -std::sqrt(.5);
    };
    std::vector<char> visited(mesh.cornerCount(), false);
    const auto walk = [&](size_t first, size_t start) {
        Chain chain;
        chain.first = start;
        size_t c = first, v = start;
        bool protectedEdge = false;
        double maximumDihedral = 0, maximumTurn = 0;
        Vector3 previous;
        while (!visited[c]) {
            const Vector3 tangent = (mesh.position(oppositeVertex(mesh, c, v)) - mesh.position(v)).normalized();
            if (!chain.corners.empty())
                maximumTurn = std::max(maximumTurn, std::acos(std::max(-1., std::min(1., dot(previous, tangent)))));
            previous = tangent;
            maximumDihedral = std::max(maximumDihedral, std::fabs(mesh.normalAngle(c)));
            visited[c] = true;
            chain.corners.push_back(c);
            chain.length += mesh.edgeVector(c).length();
            protectedEdge = protectedEdge || mesh.isBoundaryCorner(c) || std::fabs(mesh.normalAngle(c)) > hard;
            v = oppositeVertex(mesh, c, v);
            if (v == start || joint(v))
                break;
            c = incident[v][0] == c ? incident[v][1] : incident[v][0];
        }
        chain.last = v;
        chain.closed = v == start;
        // Compact automatic_classifier_rounds.ipp rule: shallow curved chains
        // lack the straight/strong evidence needed to force a grid axis.
        if (m_featureLayout && open && !protectedEdge && !chain.closed && maximumDihedral < 44.5 * M_PI / 180 && maximumTurn >= 5 * M_PI / 180)
            chain.directional = false;
        const double coherence = chain.closed ? 1.0 : (mesh.position(v) - mesh.position(start)).length() / std::max(1e-30, chain.length);
        chain.strength = protectedEdge ? 1.0 : (chain.length >= 2 * m_length && coherence >= .7 ? .5 * coherence : 0.0);
        for (size_t e : chain.corners) {
            m_features[e] = chain.strength;
            if (mesh.oppositeCorner(e) != none)
                m_features[mesh.oppositeCorner(e)] = chain.strength;
        }
        m_chains.push_back(std::move(chain));
    };
    for (size_t c : edges)
        for (size_t v : { mesh.cornerVertex(c), mesh.cornerVertex(mesh.nextCorner(c)) })
            if (!visited[c] && joint(v))
                walk(c, v);
    for (size_t c : edges)
        if (!visited[c])
            walk(c, mesh.cornerVertex(c));
}

void SurfaceAnalysis::resolveFeatureJunctions()
{
    const auto& mesh = m_mesh;
    const size_t none = SurfaceMesh::npos;
    // Resolve endpoint competition against a frozen round, so chain order cannot
    // choose the winner. Protect authored/hard chains; remove weak acute branches.
    std::vector<std::vector<size_t>> atEndpoint(mesh.vertexCount());
    for (size_t i = 0; i < m_chains.size(); ++i)
        if (!m_chains[i].closed) {
            atEndpoint[m_chains[i].first].push_back(i);
            atEndpoint[m_chains[i].last].push_back(i);
        }
    const auto outgoing = [&](const Chain& chain, size_t v) {
        const size_t c = v == chain.first ? chain.corners.front() : chain.corners.back();
        return (mesh.position(oppositeVertex(mesh, c, v)) - mesh.position(v)).normalized();
    };
    for (size_t round = 0; round < 2; ++round) {
        std::vector<double> next;
        for (const auto& chain : m_chains)
            next.push_back(chain.strength);
        for (size_t i = 0; i < m_chains.size(); ++i) {
            auto& chain = m_chains[i];
            if (chain.closed || chain.strength >= 1)
                continue;
            double support = 0;
            size_t attached = 0;
            bool crowded = false;
            for (size_t v : { chain.first, chain.last }) {
                const Vector3 d = outgoing(chain, v);
                double evidence = 0;
                for (size_t j : atEndpoint[v])
                    if (j != i && m_chains[j].strength > 0) {
                        const auto& other = m_chains[j];
                        const double alignment = dot(d, outgoing(other, v));
                        if (alignment < -.866 || std::fabs(alignment) < .259) {
                            support += other.strength;
                            evidence = std::max(evidence, other.strength);
                        }
                        const double own = chain.strength * std::min(chain.length, 4 * m_length);
                        const double rival = other.strength * std::min(other.length, 4 * m_length);
                        if (alignment > .866 && (rival > own || (rival == own && j < i)))
                            crowded = true;
                    }
                attached += evidence >= .5;
            }
            const double coherence = (mesh.position(chain.last) - mesh.position(chain.first)).length() / std::max(1e-30, chain.length);
            double base = chain.strength;
            // Rescue a short coherent link only when both ends join strong chains.
            if (base == 0 && attached == 2 && chain.length >= m_length && coherence >= .7)
                base = .25 * coherence;
            // Remesher's high-incidence classifier separates protected geometry
            // from directional chains when neither endpoint supports a grid axis.
            if (m_featureLayout && !attached && atEndpoint[chain.first].size() > 2 && atEndpoint[chain.last].size() > 2)
                chain.directional = false;
            next[i] = crowded ? 0 : std::min(.9, base * (1 + .15 * std::min(2.0, support)));
        }
        for (size_t i = 0; i < next.size(); ++i)
            m_chains[i].strength = next[i];
    }
    m_cornerChain.assign(mesh.cornerCount(), none);
    for (size_t i = 0; i < m_chains.size(); ++i)
        for (size_t c : m_chains[i].corners) {
            m_features[c] = m_chains[i].directional ? m_chains[i].strength : 0;
            m_cornerChain[c] = i;
            const size_t o = mesh.oppositeCorner(c);
            if (o != none) {
                m_features[o] = m_features[c];
                m_cornerChain[o] = i;
            }
        }
}

void SurfaceAnalysis::measureFaceCurvature()
{
    const auto& mesh = m_mesh;
    const size_t count = mesh.faceCount(), none = SurfaceMesh::npos;
    std::vector<double> areas(count);
    std::vector<Vector3> centers(count);
    std::vector<Eigen::Matrix3d> tensors(count, Eigen::Matrix3d::Zero());
    double totalArea = 0;
    for (size_t f = 0; f < count; ++f) {
        const auto& t = mesh.triangle(f);
        centers[f] = center(mesh, f);
        areas[f] = Vector3::area(mesh.position(t[0]), mesh.position(t[1]), mesh.position(t[2]));
        totalArea += areas[f];
    }
    // Ignore crease bending only on resolved planar patches. Narrow bands still
    // need its density safeguard; 2 * area / perimeter estimates their width.
    std::vector<size_t> patchOfFace(count, none);
    std::vector<char> resolvedPlanar(count, false);
    const double planeTolerance = 1e-8 * std::sqrt(totalArea);
    std::vector<size_t> patch;
    for (size_t seed = 0; seed < count; ++seed) {
        if (patchOfFace[seed] != none)
            continue;
        patch.assign(1, seed);
        patchOfFace[seed] = seed;
        const Vector3 normal = mesh.faceNormal(seed);
        double area = 0, perimeter = 0;
        for (size_t i = 0; i < patch.size(); ++i) {
            const size_t f = patch[i];
            area += areas[f];
            for (size_t c = 3 * f; c < 3 * f + 3; ++c) {
                const size_t g = mesh.adjacentFace(c);
                if (g == none || dot(normal, mesh.faceNormal(g)) < 1 - 1e-10 || std::fabs(dot(normal, centers[g] - centers[seed])) > planeTolerance || (patchOfFace[g] != none && patchOfFace[g] != seed)) {
                    perimeter += mesh.edgeVector(c).length();
                } else if (patchOfFace[g] == none) {
                    patchOfFace[g] = seed;
                    patch.push_back(g);
                }
            }
        }
        for (size_t f : patch)
            resolvedPlanar[f] = perimeter > 0 && 2 * area > m_length * perimeter;
    }
    for (size_t f = 0; f < count; ++f) {
        for (size_t c = 3 * f; c < 3 * f + 3; ++c)
            if (!mesh.isBoundaryCorner(c) && (!resolvedPlanar[f] || m_features[c] == 0)) {
                const Vector3 across = Vector3::crossProduct(mesh.faceNormal(f), mesh.edgeVector(c).normalized());
                const Eigen::Vector3d d(across.x(), across.y(), across.z());
                const Vector3 normal = mesh.faceNormal(f), acrossNormal = mesh.faceNormal(mesh.adjacentFace(c));
                const Vector3 turn = Vector3::crossProduct(normal, acrossNormal);
                // atan2 avoids acos noise on nearly coplanar faces at the finest scale.
                double angle = std::atan2(turn.length(), dot(normal, acrossNormal));
                if (dot(turn, mesh.edgeVector(c)) > 0)
                    angle = -angle;
                tensors[f] += .5 * mesh.edgeVector(c).length() * angle * (d * d.transpose());
            }
    }
    // Physical support ladder from remesher's radius schedule. Connected
    // face-distance neighborhoods do not jump across folds or feature chains.
    const double support = std::max(m_length, std::sqrt(totalArea / 10000.0));
    const double radii[scaleCount] = { .33 * support, .50 * support, .76 * support, 1.15 * support, 1.75 * support };
    tbb::parallel_for(size_t(0), count, [&](size_t f) {
        std::array<Eigen::Matrix3d, scaleCount> sum;
        for (auto& tensor : sum)
            tensor.setZero();
        double weights[scaleCount] = {};
        std::priority_queue<std::pair<double, size_t>, std::vector<std::pair<double, size_t>>, std::greater<std::pair<double, size_t>>> queue;
        std::unordered_map<size_t, double> distance;
        distance[f] = 0;
        queue.push({ 0, f });
        while (!queue.empty()) {
            const auto item = queue.top();
            queue.pop();
            const size_t g = item.second;
            if (item.first > distance[g])
                continue;
            for (size_t s = 0; s < scaleCount; ++s)
                if (item.first <= radii[s]) {
                    sum[s] += tensors[g];
                    weights[s] += areas[g];
                }
            for (size_t c = 3 * g; c < 3 * g + 3; ++c) {
                const size_t h = mesh.adjacentFace(c);
                if (h == none || m_features[c] > 0)
                    continue;
                const double d = item.first + (centers[g] - centers[h]).length();
                if (d > radii[scaleCount - 1])
                    continue;
                const auto old = distance.find(h);
                if (old == distance.end() || d < old->second) {
                    distance[h] = d;
                    queue.push({ d, h });
                }
            }
        }
        m_faces[f] = evaluateCurvatureScales(mesh, f, sum, weights, radii, false);
    });
}

// Normal-cycle edge tensors: Cohen-Steiner & Morvan (2003).
// https://doi.org/10.1145/777792.777839
// The support-radius schedule and scale-confidence rules here are separate policies.
void SurfaceAnalysis::measureVertexCurvature()
{
    const auto& mesh = m_mesh;
    const size_t count = mesh.faceCount(), none = SurfaceMesh::npos;
    // Remesher field/curvature.cpp: integrate edge bending at non-crease
    // vertices, then average only those vertices into each face. A crease's
    // tensor may be reached, but its endpoints never seed or continue a walk.
    double totalArea = 0;
    std::vector<std::vector<size_t>> vertexEdges(mesh.vertexCount());
    std::vector<char> excluded(mesh.vertexCount(), false);
    std::vector<Eigen::Matrix3d> tensors(mesh.cornerCount(), Eigen::Matrix3d::Zero());
    for (size_t f = 0; f < count; ++f) {
        const auto& t = mesh.triangle(f);
        totalArea += Vector3::area(mesh.position(t[0]), mesh.position(t[1]), mesh.position(t[2]));
        for (size_t c = 3 * f; c < 3 * f + 3; ++c) {
            const size_t o = mesh.oppositeCorner(c);
            if (o != none && o < c)
                continue;
            const size_t a = mesh.cornerVertex(c), b = mesh.cornerVertex(mesh.nextCorner(c));
            vertexEdges[a].push_back(c);
            vertexEdges[b].push_back(c);
            if (m_features[c] > 0)
                excluded[a] = excluded[b] = true;
            if (o == none)
                continue;
            const Vector3 edge = mesh.edgeVector(c), normal = mesh.faceNormal(f), other = mesh.faceNormal(mesh.cornerFace(o));
            const Vector3 turn = Vector3::crossProduct(normal, other);
            double angle = std::atan2(turn.length(), dot(normal, other));
            if (dot(turn, edge) > 0)
                angle = -angle;
            const Vector3 direction = edge.normalized();
            const Eigen::Vector3d d(direction.x(), direction.y(), direction.z());
            tensors[c] = edge.length() * angle * (d * d.transpose());
        }
    }
    const auto radii = curvatureSupportRadii(m_length, totalArea);
    std::vector<std::array<Eigen::Matrix3d, scaleCount>> vertexTensors(mesh.vertexCount());
    // Remesher curvature.cpp uses indexed edge/vertex stamps, not a hash
    // allocation at every visited edge. Keep the identical radius walk order.
    struct Walk {
        std::vector<size_t> edges, vertices, stack;
    };
    tbb::enumerable_thread_specific<Walk> walks([&] { return Walk { std::vector<size_t>(mesh.cornerCount()), std::vector<size_t>(mesh.vertexCount()), {} }; });
    tbb::parallel_for(size_t(0), mesh.vertexCount(), [&](size_t seed) {
        auto& scratch = walks.local();
        auto& walk = scratch.stack;
        for (size_t s = 0; s < scaleCount; ++s) {
            auto& sum = vertexTensors[seed][s];
            sum.setZero();
            if (excluded[seed])
                continue;
            const size_t stamp = seed * scaleCount + s + 1;
            scratch.vertices[seed] = stamp;
            walk.assign(1, seed);
            while (!walk.empty()) {
                const size_t near = walk.back();
                walk.pop_back();
                const double a = (mesh.position(near) - mesh.position(seed)).length();
                for (size_t c : vertexEdges[near])
                    if (scratch.edges[c] != stamp) {
                        scratch.edges[c] = stamp;
                        const size_t far = oppositeVertex(mesh, c, near);
                        const double b = (mesh.position(far) - mesh.position(seed)).length();
                        if (b < radii[s]) {
                            sum += tensors[c];
                            if (!excluded[far] && scratch.vertices[far] != stamp) {
                                scratch.vertices[far] = stamp;
                                walk.push_back(far);
                            }
                        } else if (b > a)
                            sum += tensors[c] * std::max(0.0, (radii[s] - a) / (b - a));
                    }
            }
            sum /= M_PI * radii[s] * radii[s];
        }
    });
    tbb::parallel_for(size_t(0), count, [&](size_t f) {
        std::array<Eigen::Matrix3d, scaleCount> sum;
        for (auto& tensor : sum)
            tensor.setZero();
        double contributors = 0;
        for (size_t v : mesh.triangle(f))
            if (!excluded[v]) {
                ++contributors;
                for (size_t s = 0; s < scaleCount; ++s)
                    sum[s] += vertexTensors[v][s];
            }
        const double weights[scaleCount] = { contributors, contributors, contributors, contributors, contributors };
        m_faces[f] = evaluateCurvatureScales(mesh, f, sum, weights, radii.data(), true);
    });
}

// Curvature-aligned anisotropic sizing background: Alliez et al. (2003).
// https://doi.org/10.1145/882262.882296
// The chord-spacing formula and limits below are this implementation's choices.
double SurfaceAnalysis::initializeDirectionalSizes(double adaptivity, double anisotropy)
{
    // Absolute radius sizing: adaptivity 1 matches remesher's 50 preset.
    // Tight radii may pass the ordinary floor (the five-sample safeguard).
    const double a = std::max(0.0, std::min(2.0, adaptivity));
    const double densityRatio = a <= 1 ? 1 + 3 * a : 4 + 8 * (a - 1);
    const double samples = a <= 1 ? 4 + 12 * a : 16 + 16 * (a - 1);
    const double chord = 2 * std::sin(M_PI / samples), safeguard = 2 * std::sin(M_PI / 5);
    const double aspect = std::pow(2.3, std::max(0.0, std::min(1.0, anisotropy)));
    const auto size = [&](double curvature) {
        const double k = curvature * m_length;
        if (a == 0 || k < 1e-12)
            return 1.0;
        return std::max(std::min(1.0, chord / k), std::min(.9 / densityRatio, safeguard / k));
    };
    for (auto& f : m_faces) {
        const double blended = (2 * f.major + f.minor) / 3;
        assignDirectionalMetric(f, size(f.confidence > .01 ? f.major : blended), size(f.confidence > .01 ? f.minor : blended), aspect);
    }
    return aspect;
}

void SurfaceAnalysis::regularizeDirectionalSizes(double aspect, bool deferSpacing)
{
    const auto& mesh = m_mesh;
    const size_t count = mesh.faceCount(), none = SurfaceMesh::npos;
    if (m_featureLayout && !deferSpacing) {
        // Remesher spacing.cpp: mean, min x2, mean, min x3, mean, min x2, mean.
        // Transport each sizing axis before propagating a neighbour's requirement.
        for (size_t pass = 0; pass < 11; ++pass) {
            auto next = m_faces;
            const bool minimum = pass != 0 && pass != 3 && pass != 7 && pass != 10;
            for (size_t f = 0; f < count; ++f) {
                double u = m_faces[f].scale * std::sqrt(m_faces[f].ratio), v = m_faces[f].scale / std::sqrt(m_faces[f].ratio), weight = 1;
                for (size_t c = 3 * f; c < 3 * f + 3; ++c) {
                    const size_t g = mesh.adjacentFace(c);
                    if (g == none || m_features[c] > 0)
                        continue;
                    const double d = dot(m_faces[f].direction, m_faces[g].direction);
                    double a = m_faces[g].scale * std::sqrt(m_faces[g].ratio), b = m_faces[g].scale / std::sqrt(m_faces[g].ratio);
                    if (d * d < .5)
                        std::swap(a, b);
                    if (minimum) {
                        u = std::min(u, a);
                        v = std::min(v, b);
                    } else {
                        u += a;
                        v += b;
                        ++weight;
                    }
                }
                assignDirectionalMetric(next[f], u / weight, v / weight, aspect);
            }
            m_faces.swap(next);
        }
    } else if (!m_featureLayout) {
        const auto limits = m_faces;
        // Regularize log sizes in transported axes; never smooth across a chain.
        for (size_t pass = 0; pass < 4; ++pass) {
            auto next = m_faces;
            for (size_t f = 0; f < count; ++f) {
                double scale = 2 * std::log(m_faces[f].scale), ratio = 2 * std::log(m_faces[f].ratio), weight = 2;
                for (size_t c = 3 * f; c < 3 * f + 3; ++c) {
                    const size_t g = mesh.adjacentFace(c);
                    if (g == none || m_features[c] > 0)
                        continue;
                    const double d = dot(m_faces[f].direction, m_faces[g].direction);
                    scale += std::log(m_faces[g].scale);
                    ratio += (2 * d * d - 1) * std::log(m_faces[g].ratio);
                    ++weight;
                }
                const double s = std::exp(scale / weight), r = std::exp(.5 * ratio / weight);
                const double bound = std::sqrt(limits[f].ratio);
                assignDirectionalMetric(next[f], std::min(s * r, limits[f].scale * bound),
                    std::min(s / r, limits[f].scale / bound), aspect);
            }
            m_faces.swap(next);
        }
    }
    // No per-island density normalization: it would cancel a thin tube's sizing.
}

void SurfaceAnalysis::buildSourceIndex()
{
    const auto& mesh = m_mesh;
    const size_t count = mesh.faceCount();
    m_boxes.resize(count);
    std::vector<size_t> indices(count);
    AxisAlignedBoudingBox bounds;
    for (size_t f = 0; f < count; ++f) {
        indices[f] = f;
        for (size_t v : mesh.triangle(f)) {
            const auto& p = mesh.position(v);
            const ::Vector3 q(p.x(), p.y(), p.z());
            m_boxes[f].update(q);
            bounds.update(q);
        }
        m_boxes[f].updateCenter();
    }
    bounds.updateCenter();
    if (count)
        m_tree.reset(new AxisAlignedBoudingBoxTree(&m_boxes, indices, bounds));
}

size_t SurfaceAnalysis::nearestFace(const Vector3& p) const
{
    size_t found = SurfaceMesh::npos;
    double best = std::numeric_limits<double>::max();
    const auto distance = [&](const AxisAlignedBoudingBox& box) {
        double d = 0;
        for (size_t k = 0; k < 3; ++k) {
            const double x = std::max({ box.lowerBound()[k] - p[k], 0.0, p[k] - box.upperBound()[k] });
            d += x * x;
        }
        return d;
    };
    std::function<void(const AxisAlignedBoudingBoxTree::Node*)> visit = [&](const AxisAlignedBoudingBoxTree::Node* node) {
        if (!node || distance(node->boundingBox) > best)
            return;
        if (node->isLeaf()) {
            for (size_t f : node->boxIndices) {
                const double d = triangleDistance2(m_mesh, f, p);
                if (d < best || (d == best && f < found)) {
                    best = d;
                    found = f;
                }
            }
        } else {
            const bool left = distance(node->left->boundingBox) < distance(node->right->boundingBox);
            visit(left ? node->left : node->right);
            visit(left ? node->right : node->left);
        }
    };
    if (m_tree)
        visit(m_tree->root());
    return found;
}

bool SurfaceAnalysis::onSourceBoundary(const Vector3& p) const
{
    // Boundary ownership must not be hidden by a closer interior feature.
    return bindCurve(p, .1 * m_length, true).chain != SurfaceMesh::npos;
}

bool SurfaceAnalysis::onSourceBoundary(const Vector3& a, const Vector3& b) const
{
    const auto first = bindCurve(a, .1 * m_length, true), second = bindCurve(b, .1 * m_length, true);
    return first.chain != SurfaceMesh::npos && first.chain == second.chain;
}

double SurfaceAnalysis::surfaceDistanceSquared(const Vector3& p, Vector3* normal) const
{
    const size_t f = nearestFace(p);
    if (normal)
        *normal = f == SurfaceMesh::npos ? Vector3() : m_mesh.faceNormal(f);
    return f == SurfaceMesh::npos ? std::numeric_limits<double>::infinity() : triangleDistance2(m_mesh, f, p);
}

// Measure the reference too: a fragment lying on the source is not a good
// approximation of a whole component. Reuse the existing triangle BVH.
double SurfaceAnalysis::missingSurfaceError(const std::vector<Vector3>& vertices,
    const std::vector<std::vector<size_t>>& faces) const
{
    std::vector<std::vector<size_t>> triangles;
    for (const auto& f : faces)
        for (size_t k = 1; k + 1 < f.size(); ++k)
            triangles.push_back({ f[0], f[k], f[k + 1] });
    if (triangles.empty())
        return std::numeric_limits<double>::infinity();
    SurfaceAnalysis output(SurfaceMesh(vertices, triangles), m_length, 180, 0, 0, false, false);
    double error = 0, area = 0;
    for (size_t f = 0; f < m_mesh.faceCount(); ++f) {
        const auto& t = m_mesh.triangle(f);
        const double a = Vector3::area(m_mesh.position(t[0]), m_mesh.position(t[1]), m_mesh.position(t[2]));
        error += a * output.surfaceDistanceSquared(center(m_mesh, f));
        area += a;
    }
    return area > 0 ? error / area : 0;
}

SurfaceAnalysis::CurveBinding SurfaceAnalysis::bindCurve(const Vector3& p, double radius, bool boundaryOnly) const
{
    CurveBinding result;
    if (!m_tree || !(radius > 0))
        return result;
    std::vector<AxisAlignedBoudingBox> query(1);
    query[0].update(::Vector3(p.x() - radius, p.y() - radius, p.z() - radius));
    query[0].update(::Vector3(p.x() + radius, p.y() + radius, p.z() + radius));
    query[0].updateCenter();
    AxisAlignedBoudingBoxTree tree(&query, { 0 }, query[0]);
    std::vector<std::pair<size_t, size_t>> hits;
    m_tree->test(m_tree->root(), tree.root(), &query, &hits);
    double best = radius * radius;
    for (const auto& hit : hits)
        for (size_t c = 3 * hit.first; c < 3 * hit.first + 3; ++c)
            if (m_cornerChain[c] != SurfaceMesh::npos && m_chains[m_cornerChain[c]].strength > 0) {
                if (boundaryOnly && !m_mesh.isBoundaryCorner(c))
                    continue;
                const Vector3 a = m_mesh.position(m_mesh.cornerVertex(c)), b = a + m_mesh.edgeVector(c);
                const double distance = (p - segmentPoint(p, a, b)).lengthSquared();
                if (distance < best || (distance == best && m_cornerChain[c] < result.chain)) {
                    best = distance;
                    result.chain = m_cornerChain[c];
                }
            }
    if (result.chain != SurfaceMesh::npos && !m_chains[result.chain].closed) {
        const auto& chain = m_chains[result.chain];
        double closest = .25 * radius * radius;
        for (size_t v : { chain.first, chain.last }) {
            const double distance = (p - m_mesh.position(v)).lengthSquared();
            if (distance < closest) {
                closest = distance;
                result.vertex = v;
            }
        }
    }
    return result;
}

Vector3 SurfaceAnalysis::projectCurve(const CurveBinding& binding, const Vector3& p) const
{
    if (binding.chain >= m_chains.size())
        return p;
    if (binding.vertex != SurfaceMesh::npos)
        return m_mesh.position(binding.vertex);
    Vector3 result = p;
    double best = std::numeric_limits<double>::max();
    for (size_t c : m_chains[binding.chain].corners) {
        const Vector3 a = m_mesh.position(m_mesh.cornerVertex(c));
        const Vector3 q = segmentPoint(p, a, a + m_mesh.edgeVector(c));
        const double distance = (q - p).lengthSquared();
        if (distance < best) {
            best = distance;
            result = q;
        }
    }
    return result;
}

// Related surface-relocation background: Botsch & Kobbelt (2004), Section 4.
// https://graphics.rwth-aachen.de/media/papers/remeshing1.pdf
// This variant uses uniform damped averaging and connected-sheet/curve constraints.
void SurfaceAnalysis::relaxSurface(std::vector<Vector3>& vertices,
    const std::vector<std::unordered_set<size_t>>& neighbors,
    const std::vector<bool>& locked, const std::vector<std::vector<size_t>>& polygons, size_t iterations) const
{
    if (!m_tree || neighbors.size() != vertices.size() || locked.size() != vertices.size())
        return;
    std::vector<std::vector<size_t>> incident(vertices.size());
    for (size_t f = 0; f < polygons.size(); ++f)
        for (size_t v : polygons[f])
            incident[v].push_back(f);
    std::vector<size_t> faces(vertices.size());
    std::vector<CurveBinding> curves(vertices.size());
    tbb::parallel_for(size_t(0), vertices.size(), [&](size_t i) {
        faces[i] = nearestFace(vertices[i]);
        double radius = m_length;
        for (size_t j : neighbors[i])
            radius = std::min(radius, (vertices[i] - vertices[j]).length());
        curves[i] = bindCurve(vertices[i], .1 * radius);
    });
    for (size_t pass = 0; pass < iterations; ++pass) {
        auto next = vertices;
        auto nextFaces = faces;
        tbb::parallel_for(size_t(0), vertices.size(), [&](size_t i) {
            if (locked[i] || neighbors[i].empty() || faces[i] == SurfaceMesh::npos)
                return;
            Vector3 target;
            for (size_t j : neighbors[i])
                target += vertices[j];
            target = (vertices[i] + target / double(neighbors[i].size())) * .5;
            if (curves[i].chain != SurfaceMesh::npos) {
                Vector3 along;
                size_t count = 0;
                for (size_t j : neighbors[i])
                    if (curves[j].chain == curves[i].chain) {
                        along += vertices[j];
                        ++count;
                    }
                target = count == 2 ? .5 * vertices[i] + .25 * along : vertices[i];
                next[i] = projectCurve(curves[i], target);
                return;
            }
            // Stay on the connected source sheet, bounded to 32 adjacent faces.
            // Feature edges are barriers; a nearby opposite sheet is never a candidate.
            size_t face = faces[i];
            Vector3 q = trianglePoint(m_mesh, face, target);
            constexpr size_t visitLimit = 32;
            std::array<size_t, 1 + 3 * visitLimit> pending, visited;
            pending[0] = visited[0] = face;
            size_t pendingCount = 1, visitedCount = 1;
            for (size_t step = 0; step < visitLimit && pendingCount; ++step) {
                const size_t current = pending[--pendingCount];
                double distance = (q - target).lengthSquared();
                if (triangleDistance2(m_mesh, current, target) > distance)
                    continue;
                for (size_t c = 3 * current; c < 3 * current + 3; ++c) {
                    const size_t f = m_mesh.adjacentFace(c);
                    if (f == SurfaceMesh::npos || m_features[c] > 0 || std::find(visited.begin(), visited.begin() + visitedCount, f) != visited.begin() + visitedCount)
                        continue;
                    const Vector3 p = trianglePoint(m_mesh, f, target);
                    const double d = (p - target).lengthSquared();
                    // Several faces can project to the same source vertex.
                    // Traverse those ties once instead of pinning the output there.
                    if (d <= distance) {
                        pending[pendingCount++] = f;
                        visited[visitedCount++] = f;
                    }
                    if (d < distance) {
                        distance = d;
                        face = f;
                        q = p;
                    }
                }
            }
            nextFaces[i] = face;
            next[i] = q;
        });
        for (size_t i = 0; i < vertices.size(); ++i)
            if (!locked[i] && safeMove(vertices, polygons, incident[i], i, next[i])) {
                vertices[i] = next[i];
                faces[i] = nextFaces[i];
            }
    }
}

size_t SurfaceAnalysis::finishCurves(std::vector<Vector3>& vertices,
    const std::vector<std::vector<size_t>>& faces, size_t iterations) const
{
    std::vector<std::vector<size_t>> incident(vertices.size());
    std::vector<std::unordered_set<size_t>> neighbors(vertices.size());
    std::vector<double> radius(vertices.size(), m_length);
    for (size_t f = 0; f < faces.size(); ++f)
        for (size_t k = 0; k < faces[f].size(); ++k) {
            const size_t a = faces[f][k], b = faces[f][(k + 1) % faces[f].size()];
            incident[a].push_back(f);
            neighbors[a].insert(b);
            neighbors[b].insert(a);
            const double length = (vertices[a] - vertices[b]).length();
            radius[a] = std::min(radius[a], length);
            radius[b] = std::min(radius[b], length);
        }
    // Fix chain identity once. A slide cannot jump to another nearby feature.
    std::vector<CurveBinding> bindings(vertices.size());
    for (size_t v = 0; v < vertices.size(); ++v)
        if (!incident[v].empty())
            bindings[v] = bindCurve(vertices[v], .25 * radius[v]);
    for (size_t pass = 0; pass <= iterations; ++pass) {
        const auto previous = vertices;
        for (size_t v = 0; v < vertices.size(); ++v) {
            if (bindings[v].chain == SurfaceMesh::npos)
                continue;
            Vector3 target = previous[v];
            if (pass && !neighbors[v].empty()) {
                Vector3 sum;
                size_t count = 0;
                for (size_t n : neighbors[v])
                    if (bindings[n].chain == bindings[v].chain) {
                        sum += previous[n];
                        ++count;
                    }
                if (count == 2)
                    target = .5 * target + .5 * sum / double(count);
            }
            const Vector3 q = projectCurve(bindings[v], target);
            if (safeMove(vertices, faces, incident[v], v, q))
                vertices[v] = q;
            else if (!pass)
                bindings[v] = CurveBinding();
        }
    }
    return std::count_if(bindings.begin(), bindings.end(), [](const CurveBinding& b) { return b.chain != SurfaceMesh::npos; });
}

double SurfaceAnalysis::scalarSize(const Vector3& p) const
{
    const size_t f = nearestFace(p);
    // Bound triangle sampling to 16x base density (~50k points at 3000 quads).
    // The quad field retains the uncapped physical metric, including tight radii.
    return f == SurfaceMesh::npos ? m_length : m_length * std::max(.25, m_faces[f].scale);
}

SurfaceGuidance SurfaceAnalysis::transfer(const SurfaceMesh& mesh, bool sameTopology) const
{
    SurfaceGuidance result;
    result.faces.resize(mesh.faceCount());
    result.featureCorners.assign(mesh.cornerCount(), 0);
    for (size_t f = 0; f < mesh.faceCount(); ++f) {
        const size_t source = sameTopology ? f : nearestFace(center(mesh, f));
        if (source == SurfaceMesh::npos)
            continue;
        result.faces[f] = m_faces[source];
        auto& face = result.faces[f];
        double strongest = 0;
        for (size_t c = 3 * f; c < 3 * f + 3; ++c) {
            const Vector3 a = mesh.position(mesh.cornerVertex(c)), edge = mesh.edgeVector(c), mid = a + edge * .5;
            const size_t near = sameTopology ? f : nearestFace(mid);
            if (near == SurfaceMesh::npos)
                continue;
            for (size_t s = sameTopology ? c : 3 * near; s < (sameTopology ? c + 1 : 3 * near + 3); ++s)
                if (m_features[s] > 0) {
                    const Vector3 p = m_mesh.position(m_mesh.cornerVertex(s)), q = p + m_mesh.edgeVector(s);
                    if (!sameTopology && ((mid - segmentPoint(mid, p, q)).length() > .1 * std::min(edge.length(), m_length) || std::fabs(dot(edge.normalized(), (q - p).normalized())) < .866))
                        continue;
                    result.featureCorners[c] = 1;
                    const double strength = m_features[s] * edge.length();
                    if (strength > strongest) {
                        strongest = strength;
                        face.direction = edge.normalized();
                        face.confidence = m_features[s];
                    }
                }
        }
        const Vector3 normal = mesh.faceNormal(f);
        face.direction = (face.direction - normal * dot(face.direction, normal)).normalized();
        // Sizes follow the original curvature axes even when a chain wins direction.
        const double alignment = dot(face.direction, m_faces[source].direction);
        face.ratio = std::pow(m_faces[source].ratio, 2 * alignment * alignment - 1);
    }
    return result;
}
}
