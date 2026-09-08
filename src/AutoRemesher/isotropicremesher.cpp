/*
 *  Copyright (c) 2020 Jeremy HU <jeremy-at-dust3d dot org>. All rights reserved.
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:

 *  The above copyright notice and this permission notice shall be included in all
 *  copies or substantial portions of the Software.

 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *  SOFTWARE.
 */
#include <AutoRemesher/IsotropicRemesher>
#include <AutoRemesher/Vector3>
#include <cstdio>

#include <isotropichalfedgemesh.h>
#include <isotropicremesher.h>

namespace AutoRemesher {

bool IsotropicRemesher::remesh()
{
    std::vector<::Vector3> inputVertices;
    inputVertices.reserve(m_vertices.size());
    for (const auto& position : m_vertices)
        inputVertices.push_back(::Vector3(position.x(), position.y(), position.z()));

    ::IsotropicRemesher remesher(&inputVertices, &m_triangles);
    if (m_targetEdgeLength > 0)
        remesher.setTargetEdgeLength(m_targetEdgeLength);
    if (m_vertexTargetEdgeLengths != nullptr)
        remesher.setVertexTargetEdgeLengths(m_vertexTargetEdgeLengths);
    remesher.setSharpEdgeIncludedAngle(180.0 - m_sharpEdgeDegrees);
    remesher.setSmoothNormalDegrees(m_smoothNormalDegrees);
    if (m_progressHandler)
        remesher.setProgressHandler(m_progressHandler);
    // Bound the feature trial to half Remesher's 50,000-point working budget.
    if (m_refineOnly)
        remesher.setRefinementVertexLimit(25000);
    remesher.remesh(m_refineOnly ? 0 : m_remeshIterations);

    IsotropicHalfedgeMesh* halfedgeMesh = remesher.remeshedHalfedgeMesh();
    if (nullptr == halfedgeMesh)
        return false;

    // Remesher prepare/refinement_coordinator.cpp balances diagonals before
    // field measurement. Reuse our flip guard only on coplanar source patches.
    if (m_refineOnly && m_balanceDiagonals) {
        halfedgeMesh->updateTriangleNormals();
        for (auto* face = halfedgeMesh->moveToNextFace(nullptr); face;) {
            auto* start = face->halfedge;
            face = halfedgeMesh->moveToNextFace(face);
            auto* edge = start;
            do {
                auto* next = edge->nextHalfedge;
                auto* other = edge->oppositeHalfedge;
                if (other && ::Vector3::dotProduct(edge->leftFace->_normal, other->leftFace->_normal) > 1 - 1e-10 && halfedgeMesh->flipEdge(edge))
                    break;
                edge = next;
            } while (edge != start);
        }
    }
    size_t outputIndex = 0;
    for (IsotropicHalfedgeMesh::Vertex* vertex = halfedgeMesh->moveToNextVertex(nullptr);
        nullptr != vertex;
        vertex = halfedgeMesh->moveToNextVertex(vertex)) {
        vertex->outputIndex = outputIndex++;
        m_remeshedVertices.push_back(Vector3 {
            vertex->position.x(),
            vertex->position.y(),
            vertex->position.z() });
    }
    for (IsotropicHalfedgeMesh::Face* face = halfedgeMesh->moveToNextFace(nullptr);
        nullptr != face;
        face = halfedgeMesh->moveToNextFace(face)) {
        m_remeshedTriangles.push_back(std::vector<size_t> {
            face->halfedge->previousHalfedge->startVertex->outputIndex,
            face->halfedge->startVertex->outputIndex,
            face->halfedge->nextHalfedge->startVertex->outputIndex });
    }

    return true;
}

void IsotropicRemesher::debugExportObj(const char* filename)
{
    FILE* fp = fopen(filename, "wb");
    for (const auto& it : m_remeshedVertices) {
        fprintf(fp, "v %f %f %f\n",
            it[0], it[1], it[2]);
    }
    for (const auto& it : m_remeshedTriangles) {
        fprintf(fp, "f %zu %zu %zu\n",
            it[0] + 1, it[1] + 1, it[2] + 1);
    }
    fclose(fp);
}

}
