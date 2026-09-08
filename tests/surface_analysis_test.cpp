#include <AutoRemesher/SurfaceAnalysis>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>
static_assert(!std::is_move_constructible<AutoRemesher::SurfaceAnalysis>::value, "analysis tree storage must stay put");
using namespace AutoRemesher;
using V = AutoRemesher::Vector3;
static void require(bool ok, const char* message)
{
    if (!ok)
        throw std::runtime_error(message);
}

static void checkMetric(const SurfaceAnalysis& analysis)
{
    for (const auto& metric : analysis.faces()) {
        require(std::isfinite(metric.scale) && metric.scale > 0 && std::isfinite(metric.ratio), "nonfinite sizing");
        require(std::isfinite(metric.confidence) && metric.confidence >= 0 && metric.confidence <= 1 + 1e-8, "invalid curvature confidence");
        require(metric.ratio >= 1 / 2.3 - 1e-8 && metric.ratio <= 2.3 + 1e-8, "aspect ratio exceeded limit");
        require(metric.scale * std::sqrt(metric.ratio) <= 1 + 1e-12 && metric.scale / std::sqrt(metric.ratio) <= 1 + 1e-12, "curvature enlarged the base spacing");
    }
}

static void cylinder(bool featureLayout = false)
{
    std::vector<V> p;
    std::vector<std::vector<size_t>> t;
    const size_t sides = 48, rings = 13;
    for (size_t z = 0; z < rings; ++z)
        for (size_t x = 0; x < sides; ++x) {
            const double a = 2 * M_PI * x / sides;
            p.emplace_back(std::cos(a), std::sin(a), double(z) / 6);
        }
    for (size_t z = 0; z + 1 < rings; ++z)
        for (size_t x = 0; x < sides; ++x) {
            size_t a = z * sides + x, b = z * sides + (x + 1) % sides, c = a + sides, d = b + sides;
            t.push_back({ a, b, d });
            t.push_back({ a, d, c });
        }
    SurfaceMesh mesh(p, t);
    SurfaceAnalysis analysis(mesh, .25, 90, 1, 1, featureLayout);
    checkMetric(analysis);
    size_t reliable = 0;
    for (size_t i = 0; i < t.size(); ++i) {
        const auto& f = analysis.faces()[i];
        const auto& v = p[t[i][0]];
        const V tangent = V(-v.y(), v.x(), 0).normalized();
        if (f.major > .5 && f.major < 2 && f.minor < .1 && std::fabs(V::dotProduct(f.direction, tangent)) > .95)
            ++reliable;
    }
    require(reliable > t.size() / 2, "cylinder curvature/directional sizing is incorrect");
    SurfaceAnalysis coarse(mesh, 1, 90, 1, 1, featureLayout);
    size_t directional = 0;
    for (const auto& f : coarse.faces())
        directional += f.ratio < .9 && f.scale < .8;
    require(directional > t.size() / 2, "physical curvature did not contract the circumference direction");
    // A uniformly thin island must gain resolution even at a fixed global budget.
    auto thin = p;
    for (auto& v : thin)
        v = v * .01;
    SurfaceAnalysis wire(SurfaceMesh(thin, t), .25, 90, 1, 1, featureLayout);
    size_t resolved = 0;
    for (const auto& f : wire.faces())
        resolved += f.scale < .15 && .25 * f.scale * std::sqrt(f.ratio) < .025;
    require(wire.scalarSize(thin[0]) >= .25 * .25, "triangle sampling exceeded its density budget");
    require(resolved > t.size() / 2, "thin island refinement was normalized away or stopped at the ordinary floor");
    // Physical scale changes curvature inversely, and leave relative sizes alone.
    for (auto& v : p)
        v = v * 7;
    SurfaceAnalysis scaled(SurfaceMesh(p, t), 1.75, 90, 1, 1, featureLayout);
    for (size_t f = 0; f < t.size(); ++f) {
        require(std::fabs(analysis.faces()[f].major - 7 * scaled.faces()[f].major) < 1e-7, "curvature is not scale covariant");
        require(std::fabs(analysis.faces()[f].scale - scaled.faces()[f].scale) < 1e-7, "sizes changed under uniform scaling");
    }
    // A coarse rim chord belongs to its source boundary even when its midpoint
    // is inside the circle; endpoints on distinct rims must not be conflated.
    require(scaled.onSourceBoundary(p[0], p[2]), "coarse rim chord lost boundary ownership");
    require(!scaled.onSourceBoundary(p[0], p.back()), "separate source rims were conflated");
    // Every round boundary is one closed chain, not many edge fragments.
    require(analysis.chains().size() == 2 && analysis.chains()[0].closed && analysis.chains()[1].closed,
        "closed boundary chains were not recognized");
}

static void creaseAndTransfer()
{
    // 45-degree crease: automatic at a 90-degree hard threshold.
    std::vector<V> p;
    std::vector<std::vector<size_t>> t;
    for (size_t i = 0; i < 9; ++i) {
        const double x = double(i) / 4;
        p.emplace_back(x, -1, 0);
        p.emplace_back(x, 0, 0);
        p.emplace_back(x, 1, 1);
    }
    for (size_t i = 0; i < 8; ++i)
        for (size_t j = 0; j < 2; ++j) {
            const size_t a = 3 * i + j, b = a + 3;
            t.push_back({ a, b, b + 1 });
            t.push_back({ a, b + 1, a + 1 });
        }
    SurfaceMesh mesh(p, t);
    SurfaceAnalysis analysis(mesh, .2, 90, 1, 1);
    require(analysis.onSourceBoundary(V(.005, 0, 0)), "interior crease hid a nearby authored boundary");
    require(!analysis.onSourceBoundary(V(1, 0, 0)), "interior crease was mistaken for a boundary");
    const auto exact = analysis.transfer(mesh, true);
    std::vector<char> expected(mesh.cornerCount(), 0);
    for (const auto& chain : analysis.chains())
        if (chain.strength > 0)
            for (size_t c : chain.corners) {
                expected[c] = 1;
                if (mesh.oppositeCorner(c) != SurfaceMesh::npos)
                    expected[mesh.oppositeCorner(c)] = 1;
            }
    require(exact.featureCorners == expected, "identity transfer changed the feature corner IDs");
    checkMetric(analysis);
    bool supported = false;
    for (const auto& c : analysis.chains())
        if (c.strength > .5 && c.strength < 1)
            supported = true;
    require(supported, "orthogonal endpoint context did not strengthen the crease");
    const auto binding = analysis.bindCurve(V(1, .01, 0), .05);
    require(binding.chain != SurfaceMesh::npos, "source crease binding failed");
    const auto slide = analysis.projectCurve(binding, V(1, -.99, 0));
    require(std::fabs(slide.y()) < 1e-12 && std::fabs(slide.z()) < 1e-12, "bound vertex hopped to another curve");
    auto working = p;
    for (size_t i = 0; i < 9; ++i)
        working[3 * i + 1] = working[3 * i + 1] + V(0, .01, 0);
    const auto before = working;
    require(analysis.finishCurves(working, t) > 0, "final curve constraints were not exercised");
    for (size_t i = 0; i < 9; ++i)
        require(std::fabs(working[3 * i + 1].y()) < 1e-10 && std::fabs(working[3 * i + 1].z()) < 1e-10, "finishing drifted away from the original crease");
    require((working[1] - p[1]).length() < 1e-10 && (working[25] - p[25]).length() < 1e-10, "junction endpoints moved");
    for (const auto& f : t) {
        const V a = V::crossProduct(before[f[1]] - before[f[0]], before[f[2]] - before[f[0]]);
        const V b = V::crossProduct(working[f[1]] - working[f[0]], working[f[2]] - working[f[0]]);
        require(V::dotProduct(a, b) > 0, "curve finishing flipped an incident triangle");
    }

    SurfaceAnalysis supportedLink(mesh, 1.2, 90, 1, 1);
    bool rescued = false;
    for (const auto& c : supportedLink.chains())
        rescued = rescued || (c.strength > 0 && c.strength < 1);
    require(rescued, "short link supported at both endpoints was not rescued");
    SurfaceAnalysis shortChain(mesh, 3, 90, 1, 1);
    for (const auto& c : shortChain.chains())
        require(c.strength == 0 || c.strength == 1, "short automatic fragment was admitted");
    const auto guidance = analysis.transfer(mesh);
    size_t selected = 0;
    for (size_t c = 0; c < mesh.cornerCount(); ++c)
        if (!mesh.isBoundaryCorner(c) && guidance.featureCorners[c])
            ++selected;
    require(selected == 16, "crease provenance did not transfer to working corners");
    // With both sizing controls off, the shared analysis gives a uniform metric.
    SurfaceAnalysis uniform(mesh, .2, 90, 0, 0);
    for (const auto& f : uniform.faces())
        require(std::fabs(f.scale - 1) < 1e-12 && f.ratio == 1, "disabled sizing still adapted");
    SurfaceAnalysis empty(SurfaceMesh({}, {}), .2, 90, 1, 1);
    require(empty.scalarSize(V()) == .2, "empty reference fallback failed");
}

static void connectedRelaxation()
{
    // Nearby disconnected sheets must not attract a vertex off its original sheet.
    std::vector<V> p;
    std::vector<std::vector<size_t>> t;
    for (double z : { 0., .02 }) {
        const size_t offset = p.size();
        for (size_t x = 0; x < 5; ++x) {
            p.emplace_back(x, -1, z);
            p.emplace_back(x, 1, z);
        }
        for (size_t x = 0; x < 4; ++x) {
            size_t a = offset + 2 * x;
            t.push_back({ a, a + 2, a + 3 });
            t.push_back({ a, a + 3, a + 1 });
        }
    }
    SurfaceAnalysis analysis(SurfaceMesh(p, t), .2, 90, 1, 1);
    std::vector<V> output = { V(.25, 0, 0), V(3.5, 0, .08) };
    std::vector<std::unordered_set<size_t>> neighbors = { { 1 }, { 0 } };
    analysis.relaxSurface(output, neighbors, { false, true }, {}, 8);
    require(std::fabs(output[0].z()) < 1e-12, "relaxation jumped to a nearby sheet");
    require(output[0].x() > 3, "source traversal stalled at a triangle edge");
    require((output[1] - V(3.5, 0, .08)).length() < 1e-12, "relaxation moved a locked vertex");
    // A captured boundary cannot slide because of an off-curve neighbor.
    output = { V(.5, -.99, 0), V(2, 0, 0) };
    analysis.relaxSurface(output, neighbors, { false, true }, {}, 8);
    require(std::fabs(output[0].y() + 1) < 1e-12 && std::fabs(output[0].z()) < 1e-12,
        "relaxation lost its source curve binding");
    require(std::fabs(output[0].x() - .5) < 1e-12, "off-curve neighbor moved a bound vertex");
    output = { V(.5, -.99, 0), V(1, -1, 0), V(3, -1, 0) };
    analysis.relaxSurface(output, { { 1, 2 }, { 0 }, { 0 } }, { false, true, true }, {}, 8);
    require(output[0].x() > 1.9 && std::fabs(output[0].y() + 1) < 1e-12, "same-curve neighbors could not slide a bound vertex");
    // Snapping the triangle apex across its opposite edge would flip this face.
    output = { V(.5, -.995, 0), V(.7, -.999, 0), V(.3, -.999, 0) };
    const auto apex = output[0];
    neighbors = { { 1, 2 }, { 0, 2 }, { 0, 1 } };
    analysis.relaxSurface(output, neighbors, { false, true, true }, { { 0, 1, 2 } }, 8);
    require((output[0] - apex).length() < 1e-12, "curve capture folded an incident face");
    // Nearby interior rows and cross-curve neighbors must not pull a feature row.
    SurfaceAnalysis stencil(SurfaceMesh(p, t), .5, 90, 0, 0, true);
    output = { { .5, -.8, 0 }, { 0, -.8, 0 }, { 1, -.8, 0 } };
    neighbors = { { 1, 2 }, { 0 }, { 0 } };
    stencil.relaxSurface(output, neighbors, { false, true, true }, {}, 2);
    require((output[0] - V(.5, -.8, 0)).length() < 1e-10, "interior vertex captured a nearby curve");
    output = { { .5, -1, 0 }, { .25, -1, 0 }, { .75, -1, 0 }, { 2, 0, 0 } };
    neighbors = { { 1, 2, 3 }, { 0 }, { 0 }, { 0 } };
    stencil.relaxSurface(output, neighbors, { false, true, true, true }, {}, 2);
    require((output[0] - V(.5, -1, 0)).length() < 1e-10, "cross-curve neighbor disturbed a regular feature row");
    // Crossing a planar fan requires visiting faces tied at its center.
    p = { { 0, 0, 0 } };
    t.clear();
    for (size_t i = 0; i < 8; ++i)
        p.emplace_back(std::cos(i * M_PI / 4), std::sin(i * M_PI / 4), 0);
    for (size_t i = 0; i < 8; ++i)
        t.push_back({ 0, 1 + i, 1 + (i + 1) % 8 });
    SurfaceAnalysis fan(SurfaceMesh(p, t), 1, 90, 0, 0, false, false);
    output = { { .2, .05, 0 }, { -.7, -.2, 0 }, { -.7, -.1, 0 } };
    fan.relaxSurface(output, { { 1, 2 }, { 0, 2 }, { 0, 1 } }, { false, true, true }, { { 0, 1, 2 } }, 1);
    require((output[0] - V(-.25, -.05, 0)).length() < 1e-12, "source traversal pinned a vertex at a planar fan center");
}

static void curvedProtectedTube()
{
    std::vector<V> p;
    std::vector<std::vector<size_t>> t;
    for (size_t j = 0; j <= 12; ++j)
        for (size_t i = 0; i < 10; ++i) {
            const double a = j * M_PI / 24, b = i * M_PI / 5, r = 2 + .1 * std::cos(b);
            p.emplace_back(r * std::cos(a), r * std::sin(a), .1 * std::sin(b));
        }
    for (size_t j = 0; j < 12; ++j)
        for (size_t i = 0; i < 10; ++i) {
            const size_t a = 10 * j + i, b = 10 * j + (i + 1) % 10;
            t.push_back({ a, b, b + 10 });
            t.push_back({ a, b + 10, a + 10 });
        }
    SurfaceMesh mesh(p, t);
    SurfaceAnalysis analysis(mesh, .1, 90, 1, 1, true);
    const auto guidance = analysis.transfer(mesh, true);
    size_t protectedCurves = 0;
    for (const auto& chain : analysis.chains())
        if (chain.strength > 0 && !chain.directional) {
            ++protectedCurves;
            const size_t c = chain.corners[chain.corners.size() / 2];
            require(!guidance.featureCorners[c], "a shallow curved tube crease forced a grid axis");
            const auto binding = analysis.bindCurve(mesh.position(mesh.cornerVertex(c)), .001);
            require(binding.chain != SurfaceMesh::npos, "protected tube curve lost its finishing provenance");
        }
    require(protectedCurves > 0, "curved tube retained faceting as grid constraints");
}

static void protectedJunctions()
{
    // Five-way icosahedron junctions cannot all be quad-field directions.
    std::vector<V> p;
    std::vector<std::vector<size_t>> triangles;
    const double phi = (1 + std::sqrt(5.)) / 2;
    for (double a : { -1., 1. })
        for (double b : { -1., 1. }) {
            p.emplace_back(0, a, b * phi);
            p.emplace_back(a, b * phi, 0);
            p.emplace_back(b * phi, 0, a);
        }
    for (size_t a = 0; a < p.size(); ++a)
        for (size_t b = a + 1; b < p.size(); ++b)
            for (size_t c = b + 1; c < p.size(); ++c)
                if (std::fabs((p[a] - p[b]).lengthSquared() - 4) < 1e-10 && std::fabs((p[b] - p[c]).lengthSquared() - 4) < 1e-10 && std::fabs((p[c] - p[a]).lengthSquared() - 4) < 1e-10) {
                    if (V::dotProduct(V::crossProduct(p[b] - p[a], p[c] - p[a]), p[a]) > 0)
                        triangles.push_back({ a, b, c });
                    else
                        triangles.push_back({ a, c, b });
                }
    SurfaceMesh mesh(p, triangles);
    SurfaceAnalysis reference(mesh, .4, 90, 1, 1, true, false);
    require(triangles.size() == 20 && reference.chains().size() == 30, "invalid junction fixture");
    const auto guidance = reference.transfer(mesh, true);
    for (char corner : guidance.featureCorners)
        require(!corner, "five-way protected junction forced a quad axis");
    const V midpoint = (p[triangles[0][0]] + p[triangles[0][1]]) * .5;
    const auto binding = reference.bindCurve(midpoint, .05);
    require(binding.chain != SurfaceMesh::npos, "protected geometry lost its original curve");
    require((reference.projectCurve(binding, midpoint) - midpoint).length() < 1e-12, "protected curve finishing moved the source edge");
}

static void missingSurface()
{
    // Every point of a surviving half lies on the source, despite losing half
    // the sheet. Reverse distance must distinguish it from a complete quad.
    const std::vector<V> p = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } };
    SurfaceAnalysis reference(SurfaceMesh(p, { { 0, 1, 2 }, { 0, 2, 3 } }), 1, 90, 0, 0);
    require(reference.surfaceDistanceSquared(V(.75, .25, 0)) < 1e-20, "coplanar fragment left its source");
    require(reference.missingSurfaceError(p, { { 0, 1, 2 } }) > .01, "missing half was scored as a complete sheet");
    require(reference.missingSurfaceError(p, { { 0, 1, 2, 3 } }) < 1e-20, "complete quad has reverse fitting error");
    require(std::isinf(reference.missingSurfaceError({}, {})), "empty output has finite fitting error");
}

int main()
{
    try {
        missingSurface();
        curvedProtectedTube();
        protectedJunctions();
        cylinder();
        cylinder(true);
        creaseAndTransfer();
        connectedRelaxation();
        std::cout << "Surface analysis checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
