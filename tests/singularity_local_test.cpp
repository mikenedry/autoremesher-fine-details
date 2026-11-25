// SPDX-License-Identifier: GPL-3.0-or-later
// Compare local charge updates and allocations with the original implementation.
#include <AutoRemesher/SingularitySimplifier>
#undef AUTO_REMESHER_SINGULARITY_SIMPLIFIER_H
#define SingularitySimplifier ReferenceSingularitySimplifier
#include <AutoRemesher/SingularitySimplifier>
#undef SingularitySimplifier
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <new>
#include <random>
#include <stdexcept>
#include <string>
// Observe full charge-vector-sized allocation requests, not elapsed time.
// The executable owns the static core, so both algorithms use this allocator.
namespace {
std::atomic<size_t> watchedBytes { 0 }, watchedRequests { 0 };
struct ChargeAllocations {
    explicit ChargeAllocations(size_t vertices)
    {
        watchedRequests = 0;
        watchedBytes = vertices * sizeof(int);
    }
    ~ChargeAllocations() { watchedBytes = 0; }
    size_t count() const { return watchedRequests.load(); }
};
}
void* operator new(size_t bytes)
{
    if (bytes == watchedBytes.load() && bytes != 0)
        ++watchedRequests;
    if (void* p = std::malloc(bytes ? bytes : 1))
        return p;
    throw std::bad_alloc();
}
void* operator new[](size_t bytes) { return ::operator new(bytes); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, size_t) noexcept { std::free(p); }
void operator delete[](void* p, size_t) noexcept { std::free(p); }
using namespace AutoRemesher;

int main()
{
    size_t cases = 0, cancelled = 0, originalAllocations = 0, localAllocations = 0;
    double originalSeconds = 0, localSeconds = 0;
    for (const size_t side : { 12u, 24u }) {
        for (bool closed : { false, true }) {
            std::vector<Vector3> points;
            std::vector<std::vector<size_t>> faces;
            for (size_t y = 0; y < side; ++y)
                for (size_t x = 0; x < side; ++x) {
                    const double u = 2 * M_PI * x / side, v = 2 * M_PI * y / side;
                    points.emplace_back(closed ? (2 + .7 * std::cos(v)) * std::cos(u) : x,
                        closed ? (2 + .7 * std::cos(v)) * std::sin(u) : y,
                        closed ? .7 * std::sin(v) : .1 * std::sin(u) * std::cos(v));
                }
            for (size_t y = 0; y < side - (closed ? 0 : 1); ++y)
                for (size_t x = 0; x < side - (closed ? 0 : 1); ++x) {
                    const size_t a = y * side + x, b = y * side + (x + 1) % side;
                    const size_t c = ((y + 1) % side) * side + x, d = ((y + 1) % side) * side + (x + 1) % side;
                    faces.push_back({ a, b, d });
                    faces.push_back({ a, d, c });
                }
            // Exercise vertices with no incident corners alongside larger regions.
            if (side == 24) {
                points.emplace_back(100, 100, 100);
                points.emplace_back(-100, -100, -100);
            }
            SurfaceMesh mesh(points, faces);
            std::vector<char> featureCorners(mesh.cornerCount(), 0);
            for (size_t c = 0; c < mesh.cornerCount(); ++c)
                if (mesh.cornerVertex(c) % side == side / 2 && mesh.cornerVertex(mesh.nextCorner(c)) % side == side / 2)
                    featureCorners[c] = 1;
            for (bool protectedEdges : { false, true }) {
                for (unsigned seed = 0; seed < 8; ++seed) {
                    const auto context = std::string(" side=") + std::to_string(side) + " closed=" + std::to_string(closed) + " protected=" + std::to_string(protectedEdges) + " seed=" + std::to_string(seed);
                    std::mt19937 rng(seed);
                    std::uniform_real_distribution<double> random(-1, 1);
                    std::vector<Vector3> reference, local;
                    for (const auto& face : faces) {
                        Vector3 normal = Vector3::crossProduct(points[face[1]] - points[face[0]], points[face[2]] - points[face[0]]).normalized();
                        Vector3 raw(random(rng), random(rng), random(rng));
                        Vector3 tangent = Vector3::crossProduct(normal, raw).normalized();
                        reference.push_back(tangent);
                    }
                    local = reference;
                    ReferenceSingularitySimplifier original(mesh, &reference);
                    SingularitySimplifier optimized(mesh, &local);
                    if (protectedEdges) {
                        original.setFeatureCorners(&featureCorners);
                        optimized.setFeatureCorners(&featureCorners);
                    }
                    {
                        ChargeAllocations probe(mesh.vertexCount());
                        const auto charges = original.vertexCharges();
                        if (charges.size() != mesh.vertexCount() || probe.count() != 1)
                            throw std::runtime_error("Charge allocation counter is not observing the core");
                    }
                    auto start = std::chrono::steady_clock::now();
                    {
                        ChargeAllocations count(mesh.vertexCount());
                        original.simplify();
                        originalAllocations += count.count();
                    }
                    originalSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
                    start = std::chrono::steady_clock::now();
                    {
                        ChargeAllocations count(mesh.vertexCount());
                        optimized.simplify();
                        localAllocations += count.count();
                    }
                    localSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
                    if (original.vertexCharges() != optimized.vertexCharges() || original.cancelledPairCount() != optimized.cancelledPairCount() || original.singularityCountBefore() != optimized.singularityCountBefore() || original.singularityCountAfter() != optimized.singularityCountAfter())
                        throw std::runtime_error("Local charge decisions differ from the original" + context);
                    for (size_t f = 0; f < local.size(); ++f)
                        if (local[f].x() != reference[f].x() ||
                            local[f].y() != reference[f].y() ||
                            local[f].z() != reference[f].z())
                            throw std::runtime_error("Optimized field differs from the original" + context);
                    cancelled += optimized.cancelledPairCount();
                    ++cases;
                }
            }
        }
    }
    std::cout << "Charge-sized requests: reference=" << originalAllocations
              << " candidate=" << localAllocations << "; seconds: reference="
              << originalSeconds << " candidate=" << localSeconds << std::endl;
    if (localAllocations >= originalAllocations)
        throw std::runtime_error("Local charge optimization did not reduce full-charge-sized allocations");
    if (!cancelled)
        throw std::runtime_error("No cancellation exercised");
    std::cout << "Exact full-field agreement across " << cases
              << " closed/open, small/large, protected/unprotected cases; cancelled pairs: "
              << cancelled << '\n';
}
