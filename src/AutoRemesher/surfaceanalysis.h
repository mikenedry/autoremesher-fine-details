#ifndef AUTO_REMESHER_SURFACE_ANALYSIS_H
#define AUTO_REMESHER_SURFACE_ANALYSIS_H
#include <AutoRemesher/SurfaceMesh>
#include <axisalignedboundingboxtree.h>
#include <memory>
#include <unordered_set>

namespace AutoRemesher {
struct SurfaceGuidance {
    struct Face {
        Vector3 direction;
        double confidence = 0, scale = 1, ratio = 1;
        double major = 0, minor = 0, radius = 0;
        bool adaptiveWeight = false;
    };
    std::vector<Face> faces;
    std::vector<char> featureCorners;
    double spacingLower = 0, spacingUpper = 0, spacingAspect = 1;
};

// One immutable analysis of an original island, sampled after mesh mutations.
class SurfaceAnalysis {
public:
    struct Chain {
        std::vector<size_t> corners; // Original island corners: provenance, not working edge IDs.
        size_t first = 0, last = 0;
        double length = 0, strength = 0;
        bool closed = false, directional = true;
    };
    SurfaceAnalysis(const SurfaceMesh& mesh, double length, double sharpDegrees,
        double adaptivity, double anisotropy, bool featureLayout = false, bool measureCurvature = true, bool deferSpacing = false);
    SurfaceAnalysis(const SurfaceAnalysis&) = delete; // The tree refers to this object's box array.
    // sameTopology requires the exact mesh used to construct this analysis.
    SurfaceGuidance transfer(const SurfaceMesh& mesh, bool sameTopology = false) const;
    double scalarSize(const Vector3& position) const;
    struct CurveBinding {
        size_t chain = SurfaceMesh::npos, vertex = SurfaceMesh::npos;
    };
    CurveBinding bindCurve(const Vector3& position, double radius, bool boundaryOnly = false) const;
    Vector3 projectCurve(const CurveBinding& binding, const Vector3& position) const;
    size_t finishCurves(std::vector<Vector3>& vertices,
        const std::vector<std::vector<size_t>>& faces, size_t iterations = 2) const;
    void relaxSurface(std::vector<Vector3>& vertices,
        const std::vector<std::unordered_set<size_t>>& neighbors,
        const std::vector<bool>& locked, const std::vector<std::vector<size_t>>& polygons, size_t iterations) const;
    const std::vector<Chain>& chains() const { return m_chains; }
    const std::vector<SurfaceGuidance::Face>& faces() const { return m_faces; }
    double length() const { return m_length; }
    bool featureLayout() const { return m_featureLayout; }
    bool onSourceBoundary(const Vector3& position) const;
    bool onSourceBoundary(const Vector3& first, const Vector3& second) const;
    double surfaceDistanceSquared(const Vector3& position, Vector3* normal = nullptr) const;
    double missingSurfaceError(const std::vector<Vector3>& vertices, const std::vector<std::vector<size_t>>& faces) const;

private:
    void traceFeatureChains(double sharpDegrees);
    void resolveFeatureJunctions();
    void measureFaceCurvature();
    void measureVertexCurvature();
    double initializeDirectionalSizes(double adaptivity, double anisotropy);
    void regularizeDirectionalSizes(double aspect, bool deferSpacing);
    void buildSourceIndex();

    SurfaceMesh m_mesh;
    double m_length;
    std::vector<Chain> m_chains;
    std::vector<double> m_features;
    std::vector<size_t> m_cornerChain;
    std::vector<SurfaceGuidance::Face> m_faces;
    bool m_featureLayout = false;
    std::vector<AxisAlignedBoudingBox> m_boxes;
    std::unique_ptr<AxisAlignedBoudingBoxTree> m_tree;
    size_t nearestFace(const Vector3& p) const;
};
}
#endif
