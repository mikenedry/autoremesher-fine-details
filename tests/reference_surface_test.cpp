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

int main(int argc, char** argv)
{
    try {
        if (argc == 1) {
            fixtures();
            tinyIsland();
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
