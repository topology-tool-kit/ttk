/// \ingroup base
/// \class ttk::ftr::Mesh
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-25
///
/// \brief TTK %FTRGraph mesh related operations
///
/// This class manage mesh request and is based on Triangulation.
/// It add the ability to pre-sort some requests.
///
/// \sa ttk::FTRGraph

#pragma once

// core includes
#include <Triangulation.h>

// local includes
#include "FTRCommon.h"
#include "FTRDataTypes.h"

// c++ includes
#include <vector>

namespace ttk {
  namespace ftr {
    // bit field havind 2 value
    // struct edgeScheme {
    //    unsigned int scheme:1;
    // };
    // keep the internal order of boundary vertices of an edge
    // 0 : 0 1
    // 1 : 1 0
    using edgeScheme = char;

    // bit field having 2^3 = 8 value.
    // only 6 are used
    //
    // keep the internal order of edges composing a triangle
    // 0 : 1 2 3
    // 1 : 1 3 2
    // 2 : 2 1 3
    // 3 : 2 3 1
    // 4 : 3 1 2
    // 5 : 3 2 1
    struct triScheme {
      unsigned int scheme : 3;

      triScheme() {
      }

      triScheme(const char u) : scheme(u) {
      }

      // direct acess to data
      operator char() const {
        return scheme;
      }
    };

    class Mesh : public Allocable {
    private:
      ttk::Triangulation *tri_;
      idVertex nbVerts_;
      idEdge nbEdges_;
      idCell nbTriangles_;
      std::vector<edgeScheme> edgesSortId_;
      std::vector<triScheme> trianglesSortId_;

    public:
      // Init
      explicit Mesh(Triangulation *tri) : tri_{tri} {
      }

      Mesh() : Mesh{nullptr} {
      }

      void setTriangulation(Triangulation *tri) {
        tri_ = tri;
      }

      Triangulation *getTriangulation() {
        return tri_;
      }

      void alloc(void) override {
        edgesSortId_.resize(nbEdges_);
        trianglesSortId_.resize(nbTriangles_);
      }

      void init(void) override {
      }

      void preprocess(void) {
        tri_->preprocessVertexNeighbors();
        tri_->preprocessVertexEdges();
        tri_->preprocessVertexTriangles();
        tri_->preprocessVertexStars();
        tri_->preprocessTriangleEdges();

        nbVerts_ = tri_->getNumberOfVertices();
        nbEdges_ = tri_->getNumberOfEdges();
        nbTriangles_ = tri_->getNumberOfTriangles();
      }

      // Facade Triangulation

      inline SimplexId getDimensionality(void) const {
        return tri_->getDimensionality();
      }

      inline idVertex getNumberOfVertices(void) const {
        return nbVerts_;
      }

      inline idEdge getNumberOfEdges(void) const {
        return nbEdges_;
      }

      inline idCell getNumberOfTriangles(void) const {
        return nbTriangles_;
      }

      inline idCell getNumberOfCells(void) const {
        return tri_->getNumberOfCells();
      }

      inline idVertex getVertexNeighborNumber(const idVertex v) const {
        return tri_->getVertexNeighborNumber(v);
      }

      inline void getVertexNeighbor(const idVertex v,
                                    const idVertex local_v,
                                    idVertex &res) {
        tri_->getVertexNeighbor(v, local_v, res);
      }

      inline idEdge getVertexEdgeNumber(const idEdge e) const {
        return tri_->getVertexEdgeNumber(e);
      }

      inline void getVertexEdge(const idVertex v,
                                const idEdge local_e,
                                idEdge &res) const {
        tri_->getVertexEdge(v, local_e, res);
      }

      inline idCell getVertexTriangleNumber(const idVertex v) const {
        return tri_->getVertexTriangleNumber(v);
      }

      inline void getVertexTriangle(const idVertex v,
                                    const idCell local_c,
                                    idCell &res) const {
        tri_->getVertexTriangle(v, local_c, res);
      }

      inline void getEdgeVertex(const idEdge e,
                                const char local_id,
                                idVertex &res) const {
        tri_->getEdgeVertex(e, local_id, res);
      }

      inline idCell getEdgeTriangleNumber(const idEdge e) const {
        return tri_->getEdgeTriangleNumber(e);
      }

      inline void getEdgeTriangle(const idEdge e,
                                  const idCell local_t,
                                  idCell &res) const {
        tri_->getEdgeTriangle(e, local_t, res);
      }

      inline idEdge getTriangleEdgeNumber(const idCell t) const {
        return tri_->getTriangleEdgeNumber(t);
      }

      inline void getTriangleEdge(const idCell t,
                                  const char local_id,
                                  idEdge &res) const {
        tri_->getTriangleEdge(t, local_id, res);
      }

      // Extend Triangulation by pre-sorting scalars
      // This is done by keeping id of the internal scheme of the simplex to
      // reorder it instantly

      // precompute

      void preSortEdges(VertCompFN lowerThan);

      void preSortTriangles(VertCompFN lowerThan);

      // get simplex

      orderedEdge getOrderedEdge(const idEdge e,
                                 const bool increasingOrder) const;

      void getOrderedEdge(const idEdge e,
                          const bool increasingOrder,
                          orderedEdge &oEdge) const;

      orderedTriangle getOrderedTriangle(const idCell t,
                                         const bool increasingOrder) const;

      void getOrderedTriangle(const idCell t,
                              const bool increasingOrder,
                              orderedTriangle &oTri) const;

      bool compareEdges(const idEdge e0,
                        const idEdge e1,
                        VertCompFN lowerThan) const;

      bool compareLinks(const linkEdge &l0,
                        const linkEdge &l1,
                        VertCompFN lowerTan) const;

      // tools

      std::string printEdges(void) const;

      std::string printEdge(const idEdge e) const;
    };
  } // namespace ftr
} // namespace ttk
