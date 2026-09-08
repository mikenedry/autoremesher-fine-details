// Geometry regressions and an optional headless OBJ driver.
#include <AutoRemesher/AutoRemesher>
#include <AutoRemesher/FrameField>
#include <AutoRemesher/IsotropicRemesher>
#include <AutoRemesher/MeshSeparator>
#include <AutoRemesher/MixedIntegerLeastSquares>
#include <AutoRemesher/QuadExtractor>
#include <AutoRemesher/QuadParameterizer>
#include <AutoRemesher/SurfaceAnalysis>
#include <isotropichalfedgemesh.h>
#include <isotropicremesher.h>
#include <map>
#include <set>
#define TINYOBJLOADER_IMPLEMENTATION
#include "../src/tiny_obj_loader.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <type_traits>

using Remesher = AutoRemesher::AutoRemesher;
using Vec3 = AutoRemesher::Vector3;
using Faces = std::vector<std::vector<size_t>>;
using Feature = Remesher::ReferenceSurface::EdgeFeature;
using UV = AutoRemesher::Vector2;

static void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

static void tube(size_t sides, size_t rings, std::vector<Vec3>& p, Faces& t,
    std::vector<std::vector<UV>>& uv)
{
    for (size_t j = 0; j < rings; ++j)
        for (size_t i = 0; i < sides; ++i)
            p.emplace_back(std::cos(2 * M_PI * i / sides), std::sin(2 * M_PI * i / sides), double(j));
    for (size_t j = 0; j + 1 < rings; ++j)
        for (size_t i = 0; i < sides; ++i) {
            const size_t a = j * sides + i, b = j * sides + (i + 1) % sides, c = b + sides, d = a + sides;
            t.push_back({ a, b, c });
            uv.push_back({ { double(i), double(j) }, { double(i + 1), double(j) }, { double(i + 1), double(j + 1) } });
            t.push_back({ a, c, d });
            uv.push_back({ { double(i), double(j) }, { double(i + 1), double(j + 1) }, { double(i), double(j + 1) } });
        }
}

static void checkReferences(Remesher& remesher, const std::vector<Vec3>& vertices,
    const Faces& faces)
{
    std::vector<bool> found(faces.size(), false);
    size_t vertexOffset = 0, triangleOffset = 0;
    for (const auto& prepared : remesher.preparedIslands()) {
        require(bool(prepared.reference), "missing reference");
        const auto& ref = *prepared.reference;
        require(bool(prepared.analysis), "missing persistent source analysis");
        require(prepared.analysis->faces().size() == ref.triangles.size(), "source metric domain");
        for (const auto& chain : prepared.analysis->chains()) {
            require(chain.first < ref.vertices.size() && chain.last < ref.vertices.size(), "chain endpoint provenance");
            for (size_t corner : chain.corners)
                require(corner < 3 * ref.triangles.size(), "chain corner provenance");
        }
        static_assert(std::is_const<typename std::remove_reference<decltype(ref)>::type>::value,
            "reference surface must be immutable");
        require(ref.vertices.size() == ref.sourceVertexIds.size(), "vertex provenance size");
        require(ref.triangles.size() == ref.sourceTriangleIds.size(), "triangle provenance size");
        require(ref.edgeFeatures.size() == 3 * ref.triangles.size(), "feature provenance size");
        for (size_t i = 0; i < ref.vertices.size(); ++i) {
            require(ref.sourceVertexIds[i] < vertices.size(), "source vertex domain");
            require((ref.vertices[i] - vertices[ref.sourceVertexIds[i]]).lengthSquared() == 0, "reference position changed");
        }
        for (size_t i = 0; i < ref.triangles.size(); ++i) {
            const size_t source = ref.sourceTriangleIds[i];
            require(source < faces.size() && !found[source], "lost or duplicated source triangle");
            found[source] = true;
            for (size_t corner = 0; corner < 3; ++corner) {
                const size_t local = ref.triangles[i][corner];
                require(local < ref.vertices.size(), "reference vertex domain");
                require(ref.sourceVertexIds[local] == faces[source][corner], "source corner order changed");
            }
        }
        require(prepared.vertexOffset == vertexOffset && prepared.triangleOffset == triangleOffset,
            "prepared island offsets");
        vertexOffset += prepared.vertexCount;
        triangleOffset += prepared.triangleCount;
        require(triangleOffset <= remesher.isotropicTriangles().size(), "prepared triangle range");
        for (size_t i = prepared.triangleOffset; i < triangleOffset; ++i)
            for (size_t v : remesher.isotropicTriangles()[i])
                require(v >= prepared.vertexOffset && v < vertexOffset, "working triangle crosses reference islands");
    }
    for (bool seen : found)
        require(seen, "source triangle missing after preparation");
    require(vertexOffset == remesher.isotropicVertices().size(), "unmapped working vertices");
    require(triangleOffset == remesher.isotropicTriangles().size(), "unmapped working triangles");
}

static void fixtures()
{
    // Two coincident, disconnected folded strips with different input indices.
    // Both contain a 90-degree crease and four open boundary edges.
    std::vector<Vec3> vertices { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 },
        { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
    const auto original = vertices;
    const Faces faces { { 0, 1, 2 }, { 1, 0, 3 }, { 4, 5, 6 }, { 5, 4, 7 } };
    Remesher remesher(vertices, faces);
    vertices[0] = { 99, 99, 99 }; // The caller's array is not the reference owner.
    remesher.setTargetTriangleCount(100);
    remesher.setSharpEdgeDegrees(30);
    remesher.remesh();
    checkReferences(remesher, original, faces);
    require(remesher.preparedIslands().size() == 2, "coincident islands aliased");
    for (const auto& prepared : remesher.preparedIslands()) {
        const auto& ref = *prepared.reference;
        size_t sharp = 0, boundary = 0;
        for (Feature edge : ref.edgeFeatures) {
            sharp += edge == Feature::Sharp;
            boundary += edge == Feature::Boundary;
        }
        require(sharp == 2 && boundary == 4, "original crease or boundary classification lost");
        require(ref.sharpEdgeDegrees == 30, "feature threshold provenance");
        require(prepared.vertexCount != ref.vertices.size(), "fixture did not resample");
    }
    const auto retained = remesher.preparedIslands().front().reference;
    remesher.setSharpEdgeDegrees(120);
    remesher.remesh();
    checkReferences(remesher, original, faces);
    require(retained->sharpEdgeDegrees == 30, "previous snapshot mutated on rerun");
    for (Feature edge : remesher.preparedIslands().front().reference->edgeFeatures)
        require(edge != Feature::Sharp, "new threshold not captured");
    remesher.setTargetTriangleCount(0);
    require(!remesher.remesh() && remesher.preparedIslands().empty(), "failed run exposes stale preparation");
    require(remesher.remeshedVertices().empty() && remesher.remeshedQuads().empty() && remesher.isotropicVertices().empty() && remesher.isotropicTriangles().empty() && remesher.isotropicTriangleUvs().empty(), "failed run exposes stale delivery");
    require(retained->vertices[0][0] == original[0][0], "retained reference expired");

    // Equal face values must retain distinct input row IDs.
    const Faces duplicates { { 0, 1, 2 }, { 0, 1, 2 } };
    std::vector<Faces> islands;
    std::vector<std::vector<size_t>> ids;
    AutoRemesher::MeshSeparator::splitToIslands(duplicates, islands, &ids);
    require(ids.size() == 2 && ids[0][0] == 0 && ids[1][0] == 1, "duplicate source face identity");
}

static void tinyIsland()
{
    const std::vector<Vec3> cube { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
        { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
    const Faces quads { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
    std::vector<Vec3> vertices = cube;
    Faces triangles;
    for (const auto& p : cube)
        vertices.push_back(p * .01 + Vec3(3, 0, 0));
    for (size_t offset : { size_t(0), size_t(8) })
        for (const auto& q : quads) {
            triangles.push_back({ q[0] + offset, q[1] + offset, q[2] + offset });
            triangles.push_back({ q[0] + offset, q[2] + offset, q[3] + offset });
        }
    Remesher remesher(vertices, triangles);
    remesher.setTargetTriangleCount(2000);
    remesher.setSharpEdgeDegrees(80);
    require(remesher.remesh(), "tiny island regression failed to remesh");
    bool large = false, small = false;
    for (const auto& q : remesher.remeshedQuads()) {
        const auto& p = remesher.remeshedVertices()[q[0]];
        large = large || p.x() < 2;
        small = small || p.x() > 2;
    }
    require(large && small, "global sizing removed the tiny disconnected part");
}

static void directionalQuadCover()
{
    // A very uneven triangulation must reproduce an affine rectangular grid,
    // including anisotropic spacing and a change of world units.
    const Faces triangles = { { 0, 1, 4 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 } };
    for (double worldScale : { 1., 7. })
        for (double uSpacing : { 1., 2. })
            for (bool adaptiveSizing : { false, true }) {
                std::vector<Vec3> positions = { { 0, 0, 0 }, { 4, 0, 0 }, { 4, 2, 0 }, { 0, 2, 0 }, { .07, .73, 0 } };
                for (auto& p : positions)
                    p = p * worldScale;
                const AutoRemesher::SurfaceMesh mesh(positions, triangles);
                const std::vector<Vec3> field(4, Vec3(1, 0, 0));
                const std::vector<double> u(4, uSpacing), v(4, 1);
                // Local detail may resolve below the component's nominal spacing floor.
                // Constant physical targets must survive the final transported sizing
                // stage, including unequal triangle areas and a change of world units.
                AutoRemesher::SurfaceGuidance sizing;
                sizing.faces.resize(4);
                sizing.spacingUpper = 10;
                sizing.spacingLower = 2.5;
                sizing.spacingAspect = 2.3;
                for (auto& face : sizing.faces)
                    face.adaptiveWeight = true;
                const std::vector<double> scalar(4, 1);
                AutoRemesher::QuadParameterizer::Result result;
                if (worldScale == 1 && uSpacing == 1 && !adaptiveSizing) {
                    const std::vector<char> shortCorners;
                    const std::vector<double> shortScale;
                    AutoRemesher::SurfaceGuidance shortGuidance;
                    std::vector<Vec3> outputField;
                    require(!AutoRemesher::QuadParameterizer::parameterize(positions, triangles, &field, 1, 90, &result,
                                nullptr, nullptr, nullptr, nullptr, &shortCorners, true),
                        "short feature array accepted");
                    require(!AutoRemesher::QuadParameterizer::parameterize(positions, triangles, &field, 1, 90, &result,
                                &shortScale, nullptr, nullptr, nullptr, nullptr, true),
                        "short sizing array accepted");
                    require(!AutoRemesher::FrameField::create(mesh, 90, &outputField, &shortGuidance), "short field guidance accepted");
                }
                require(AutoRemesher::QuadParameterizer::parameterize(positions, triangles, &field,
                            worldScale / mesh.averageEdgeLength(), 90, &result, adaptiveSizing ? &scalar : nullptr, &u, &v,
                            nullptr, nullptr, adaptiveSizing, adaptiveSizing ? &sizing : nullptr),
                    "directional quad cover failed");
                for (size_t f = 0; f < triangles.size(); ++f)
                    for (size_t k = 0; k < 3; ++k) {
                        const auto& p = positions[triangles[f][k]];
                        const auto& uv = result.triangleUvs[f][k];
                        require(std::fabs(uv.x() - p.x() / (worldScale * uSpacing)) < 1e-5 && std::fabs(uv.y() - p.y() / worldScale) < 1e-5, adaptiveSizing ? "final spacing changed a resolved affine grid" : "quad cover distorted an affine grid");
                    }
            }
}

static void oddTubePeriod()
{
    // Five cells around a tube require an odd chart translation. An even-only
    // lattice changes the circumference and can collapse narrow curved tubes.
    std::vector<Vec3> p;
    Faces t;
    std::vector<std::vector<UV>> uv;
    const size_t sides = 16;
    tube(sides, 7, p, t, uv);
    AutoRemesher::SurfaceMesh mesh(p, t);
    std::vector<Vec3> field(t.size(), Vec3(0, 0, 1));
    std::vector<double> u(t.size(), 1), v(t.size(), 2 * sides * std::sin(M_PI / sides) / 5);
    AutoRemesher::QuadParameterizer::Result result;
    require(AutoRemesher::QuadParameterizer::parameterize(p, t, &field, 1 / mesh.averageEdgeLength(),
                90, &result, nullptr, &u, &v, nullptr, nullptr, true),
        "odd tube cover failed");
    double circumference = 0;
    for (size_t i = 0; i < sides; ++i) {
        const auto d = result.triangleUvs[2 * i][1] - result.triangleUvs[2 * i][0];
        circumference += std::hypot(d.x(), d.y());
    }
    require(std::fabs(circumference - 5) < 1e-5, "tube circumference was forced onto an even lattice");
}

static void featureIntegerLattice()
{
    using System = AutoRemesher::MixedIntegerLeastSquares;
    // Independent near-half-integer targets admit only one variable per pass.
    // The cover's iteration allowance must scale beyond its old 100-pass cap.
    System slow(201);
    for (size_t i = 0; i < 201; ++i) {
        slow.setVariablePeriod(i, 1);
        slow.addEnergy({ { i, 1 } }, .49);
    }
    slow.finalizeConstraints();
    const size_t limit = slow.integerKernelVariableCount() + 1;
    for (size_t i = 0; i < limit && !slow.converged(); ++i)
        require(slow.solveIteration(true, true), "bounded integer rounding failed");
    require(slow.converged(), "integer rounding stopped before admitting every variable");
    for (bool star : { false, true }) {
        System system(6);
        for (size_t i = 0; i < 3; ++i) {
            system.setVariablePeriod(i, 1);
            system.setVariablePeriod(i + 3, 2);
            system.addConstraint(i, 1, i + 3, -1);
            system.addEnergy({ { i, 1 } }, .1 * i);
        }
        system.finalizeConstraints();
        require(system.solveIteration(false), "continuous lattice solve");
        require(system.separateIntegerCoordinates(0, 1, 1), "first feature separation rejected");
        require(system.separateIntegerCoordinates(star ? 0 : 1, 2, 1), "adjacent feature separation rejected after elimination");
        for (size_t i = 0; i < 100 && !system.converged(); ++i)
            require(system.solveIteration(true, true), "budgeted lattice solve failed");
        require(system.converged(), "budgeted rounding failed to converge");
        require(std::fabs(system.value(1) - system.value(0) - 2) < 1e-8 && std::fabs(system.value(2) - system.value(star ? 0 : 1) - 2) < 1e-8, "positive feature lattice spacing lost");
        for (size_t i = 0; i < 6; ++i)
            require(std::fabs(system.value(i) / 2 - std::round(system.value(i) / 2)) < 1e-8, "period-two cover lost integrality");
        require(!system.separateIntegerCoordinates(99, 0, 1), "out of range feature accepted");
    }
    // A sparse solve can land on either side of a half-integer by roundoff.
    // Original stored-float decisions must use the same tie for scoring and fixing.
    for (double target : { 1. - 1e-12, 1. + 1e-12 }) {
        System tie(1);
        tie.setVariablePeriod(0, 2);
        tie.addEnergy({ { 0, 1 } }, target);
        tie.finalizeConstraints();
        require(tie.solveIteration(false) && tie.solveIteration(true, true) && tie.converged(), "half-integer solve failed");
        require(std::fabs(tie.value(0) - 2) < 1e-8, "half-integer decision changed with solver roundoff");
    }
}

static void radialCover()
{
    std::vector<Vec3> p { { 0, 0, 0 } };
    Faces t;
    const size_t sides = 32;
    for (double r : { .35, .4, .45, 1. })
        for (size_t i = 0; i < sides; ++i)
            p.emplace_back(r * std::cos(2 * M_PI * i / sides), r * std::sin(2 * M_PI * i / sides), 0);
    for (size_t i = 0; i < sides; ++i)
        t.push_back({ 0, 1 + i, 1 + (i + 1) % sides });
    for (size_t j = 0; j < 3; ++j)
        for (size_t i = 0; i < sides; ++i) {
            const size_t a = 1 + j * sides + i, b = 1 + j * sides + (i + 1) % sides, c = b + sides, d = a + sides;
            t.push_back({ a, d, c });
            t.push_back({ a, c, b });
        }
    std::vector<Vec3> field;
    std::vector<char> features;
    for (const auto& f : t) {
        field.push_back((p[f[0]] + p[f[1]] + p[f[2]]).normalized());
        for (size_t k = 0; k < 3; ++k)
            features.push_back(std::fabs(p[f[k]].length() - p[f[(k + 1) % 3]].length()) < 1e-8);
    }
    AutoRemesher::SurfaceMesh mesh(p, t);
    AutoRemesher::QuadParameterizer::Result result;
    require(AutoRemesher::QuadParameterizer::parameterize(p, t, &field, .25 / mesh.averageEdgeLength(), 90, &result,
                nullptr, nullptr, nullptr, nullptr, &features, true),
        "radial cover failed");
    require(result.fullTurnVertices == std::vector<size_t> { 0 }, "radial full-turn center was lost modulo four");
    double area = 0;
    for (const auto& q : result.triangleUvs) {
        const double a = (q[1].x() - q[0].x()) * (q[2].y() - q[0].y()) - (q[1].y() - q[0].y()) * (q[2].x() - q[0].x());
        require(a >= -1e-6, "radial cover folded");
        area += a;
    }
    require(area > 10, "radial grid collapsed to a line");
}

int main(int argc, char** argv)
{
    try {
        if (argc == 1) {
            fixtures();
            tinyIsland();
            directionalQuadCover();
            oddTubePeriod();
            featureIntegerLattice();
            radialCover();
        } else {
            require(argc == 4, "usage: reference_surface_test [input.obj target_quads output.obj]");
            tinyobj::attrib_t attributes;
            std::vector<tinyobj::shape_t> shapes;
            std::vector<tinyobj::material_t> materials;
            std::string warning, error;
            require(tinyobj::LoadObj(&attributes, &shapes, &materials, &warning, &error, argv[1]), "OBJ load failed");
            std::vector<Vec3> vertices;
            for (size_t i = 0; i < attributes.vertices.size(); i += 3)
                vertices.push_back({ attributes.vertices[i], attributes.vertices[i + 1], attributes.vertices[i + 2] });
            Faces triangles;
            for (const auto& shape : shapes)
                for (size_t i = 0; i < shape.mesh.indices.size(); i += 3)
                    triangles.push_back({ size_t(shape.mesh.indices[i].vertex_index),
                        size_t(shape.mesh.indices[i + 1].vertex_index), size_t(shape.mesh.indices[i + 2].vertex_index) });
            Remesher remesher(vertices, triangles);
            remesher.setTargetTriangleCount(std::stoull(argv[2]) * 2);
            remesher.setScaling(1.0);
            remesher.setSharpEdgeDegrees(90.0);
            require(remesher.remesh(), "remesh failed");
            checkReferences(remesher, vertices, triangles);
            require(!remesher.remeshedQuads().empty(), "empty result");
            for (const auto& p : remesher.remeshedVertices())
                for (size_t axis = 0; axis < 3; ++axis)
                    require(std::isfinite(p[axis]), "non-finite output vertex");
            for (const auto& face : remesher.remeshedQuads()) {
                require(face.size() >= 3, "short output face");
                for (size_t v : face)
                    require(v < remesher.remeshedVertices().size(), "output face index out of range");
            }
            std::ofstream out(argv[3]);
            out << std::setprecision(17);
            for (const auto& p : remesher.remeshedVertices())
                out << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
            for (const auto& f : remesher.remeshedQuads()) {
                out << 'f';
                for (size_t v : f)
                    out << ' ' << v + 1;
                out << '\n';
            }
            require(out.good(), "OBJ write failed");
            size_t chains = 0, selected = 0;
            for (const auto& island : remesher.preparedIslands())
                for (const auto& chain : island.analysis->chains()) {
                    ++chains;
                    selected += chain.strength > 0;
                }
            std::cout << "Feature chains: " << selected << " selected of " << chains << '\n';
            std::cout << "Reference islands: " << remesher.preparedIslands().size()
                      << ", decimation used: " << remesher.decimated() << '\n';
        }
        std::cout << "Reference checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
