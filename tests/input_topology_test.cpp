// SPDX-License-Identifier: MIT
#include <AutoRemesher/AutoRemesher>
#include "../src/AutoRemesher/inputtopology.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <set>
#include <stdexcept>
using V = AutoRemesher::Vector3;
using Faces = std::vector<std::vector<size_t>>;
static void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
static void classification()
{
    const std::vector<V> base {{0,0,0}, {1,0,0}, {0,4,.08}, {1,-4,.08}, {.5,0,.004}};
    const Faces faces {{0,1,2}, {1,0,3}, {0,1,4}};
    using AutoRemesher::redundantInputFins;
    require(redundantInputFins(base, faces) == std::vector<size_t>{2}, "isolated negligible fin not recognized");
    require(redundantInputFins(base, {{0,1,2},{1,0,3}}).empty(), "open sheet changed");
    require(redundantInputFins(base, {{0,1,4}}).empty(), "isolated thin triangle changed");
    require(redundantInputFins(base, {{0,1,2},{1,0,4}}).empty(), "manifold thin detail changed");
    require(redundantInputFins(base, {{0,2,1},{0,1,3},{1,2,3},{2,0,3}}).empty(), "closed surface changed");
    require(redundantInputFins(base, {{0,1,2},{0,1,3},{0,1,4}}).empty(), "inconsistent winding repaired speculatively");
    require(redundantInputFins(base, {{0,1,2},{1,0,3},{0,1,4},{1,0,2}}).empty(), "four-face junction changed");
    for (double height : {0., .0049}) {
        auto v=base; v[4]=V(.5,0,height);
        require(redundantInputFins(v,faces)==std::vector<size_t>{2}, "negligible fin below bound not recognized");
    }
    for (double height : {.0051, .006, .1, 1.}) {
        auto v=base; v[4]=V(.5,0,height);
        require(redundantInputFins(v,faces).empty(), "meaningful raised fin removed");
    }
    for (double position : {-.1,0.,1.,1.1}) {
        auto v=base; v[4]=V(position,0,.004);
        require(redundantInputFins(v,faces).empty(), "fin beyond attachment interior removed");
    }
    // An outward branch is not an internal fold, even if it is equally small.
    auto v=base; v[4]=V(.5,0,-.004);
    require(redundantInputFins(v,faces).empty(), "outward negligible feature removed");
    v=base; v[2]=V(0,4,0); v[3]=V(1,-4,0);
    require(redundantInputFins(v,faces).empty(), "flat supporting sheet has no resolved inward side");
    v[4]=V(.5,0,0);
    require(redundantInputFins(v,faces)==std::vector<size_t>{2}, "exact collinear fin on flat sheet retained");
    v=base; v[2]=V(0,4,0); v[3]=V(1,-4,1e-15);
    require(redundantInputFins(v,faces).empty(), "unresolved supporting bend classified");
    v=base; v[4]=V(.5,.004,0);
    require(redundantInputFins(v,faces).empty(), "branch outside one supporting half-space removed");
    v=base; v[2]=V(0,1,0); v[3]=V(1,-1,0);
    require(redundantInputFins(v,faces).empty(), "fin large relative to supporting faces removed");
    v=base; v[3]=V(1,-4,1);
    require(redundantInputFins(v,faces).empty(), "crease treated as a continuous sheet");
    v=base; v[3]=V(1,4,0);
    require(redundantInputFins(v,faces).empty(), "overlapping supports accepted");
    v=base; v.push_back(V(.5,1,1)); v.push_back(V(.6,1,1));
    auto connected=faces; connected.push_back({4,5,6});
    require(redundantInputFins(v,connected).empty(), "connected branch removed");
    v=base; v[1]=v[0];
    require(redundantInputFins(v,faces).empty(), "zero-length attachment changed");
    v=base; v[4]=V(std::numeric_limits<double>::quiet_NaN(),0,0);
    require(redundantInputFins(v,faces).empty(), "nonfinite input classified");
    require(redundantInputFins(base, {{0,1,99}}).empty(), "invalid index classified");
    // A ring with an authored central opening must remain an open ring.
    const std::vector<V> ring {{-2,-2,0},{2,-2,0},{2,2,0},{-2,2,0},
                              {-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0}};
    const Faces ringFaces {{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}};
    require(redundantInputFins(ring,ringFaces).empty(), "authored opening changed");
    for (double scale : {1e-9,1.,1e9}) {
      for (int variant : {0,1,2,3}) {
        std::vector<V> transformed;
        for (size_t vertex=0; vertex<base.size(); ++vertex) {
            auto p=base[vertex];
            if (variant==1 && vertex==4) p=V(.5,0,-.004);
            if (variant==2 && vertex==2) p=V(0,4,0);
            if (variant==2 && vertex==3) p=V(1,-4,0);
            if (variant==3) {
                if (vertex==4) p=V(.5,0,0);
                // Axis permutation preserves exact collinearity in floating point.
                transformed.emplace_back((p.z()+1)*scale,(p.x()+2)*scale,(p.y()-3)*scale);
                continue;
            }
            const double x=.8*p.x()-.6*p.y(), y=.6*p.x()+.8*p.y();
            transformed.emplace_back((x+2)*scale,(.8*y-.6*p.z()-3)*scale,(.6*y+.8*p.z()+1)*scale);
        }
        std::vector<size_t> order {0,1,2};
        do {
            for (bool reverse : {false,true}) {
                Faces shuffled; size_t fin=0;
                for (size_t i=0;i<order.size();++i) {
                    auto f=faces[order[i]]; std::rotate(f.begin(),f.begin()+i,f.end());
                    if (reverse) std::reverse(f.begin(),f.end());
                    shuffled.push_back(f); if (order[i]==2) fin=i;
                }
                const auto expected=(variant==1 || variant==2) ? std::vector<size_t>{} : std::vector<size_t>{fin};
                require(redundantInputFins(transformed,shuffled)==expected,
                        "classification depends on scale, rotation, winding or face order");
            }
        } while (std::next_permutation(order.begin(),order.end()));
      }
    }
}
int main()
{
    try {
        classification();
        const std::vector<V> vertices {{0,0,0}, {1,0,0}, {0,4,.08}, {1,-4,.08}, {.5,0,.004}};
        const Faces faces {{0,1,2}, {0,1,4}, {1,0,3}};
        AutoRemesher::AutoRemesher engine(vertices, faces);
        engine.setTargetTriangleCount(64);
        engine.remesh();
        require(engine.removedInputFinTriangles()==std::vector<size_t>{1}, "repair IDs lost their original input numbering");
        std::set<size_t> retained;
        size_t boundaries=0;
        for (const auto& island : engine.preparedIslands()) {
            const auto& ref=*island.reference;
            for (auto edge : ref.edgeFeatures)
                boundaries += edge==AutoRemesher::AutoRemesher::ReferenceSurface::EdgeFeature::Boundary;
            for (size_t f=0;f<ref.triangles.size();++f)
                for (size_t c=0;c<3;++c) {
                    const auto original=faces[ref.sourceTriangleIds[f]][c];
                    require(ref.sourceVertexIds[ref.triangles[f][c]]==original, "retained face source mapping changed");
                    require((ref.vertices[ref.triangles[f][c]]-vertices[original]).lengthSquared()==0,
                            "input correction moved a retained vertex");
                }
        }
        require(boundaries==4, "repaired reference gained a boundary inside the continuous sheet");
        for (const auto& island : engine.preparedIslands())
            for (size_t id : island.reference->sourceTriangleIds) retained.insert(id);
        require(retained == std::set<size_t>({0,2}),
            "dangling near-collinear fin became an authored source boundary");
        require(faces.size() == 3 && vertices.size() == 5, "input was modified");
        std::cout << "Input topology regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
