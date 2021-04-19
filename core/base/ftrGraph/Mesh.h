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

    template <typename triangulationType>
    class Mesh : public Allocable {
    private:
      triangulationType *tri_{};
      idVertex nbVerts_{};
      idEdge nbEdges_{};
      idCell nbTriangles_{};
      std::vector<edgeScheme> edgesSortId_{};
      std::vector<triScheme> trianglesSortId_{};

    public:
      // Init
      explicit Mesh(triangulationType *tri) : tri_{tri} {
      }

      Mesh() {
      }

      void setTriangulation(triangulationType *tri) {
        tri_ = tri;
      }

      triangulationType *getTriangulation() {
        return tri_;
      }

      void alloc(void) override {
        edgesSortId_.resize(nbEdges_);
        trianglesSortId_.resize(nbTriangles_);
      }

      void init(void) override {
      }

      void preprocess(void) {
        tri_->preconditionVertexNeighbors();
        tri_->preconditionVertexEdges();
        tri_->preconditionVertexTriangles();
        tri_->preconditionVertexStars();
        tri_->preconditionTriangleEdges();

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

      void preSortEdges(const VertCompFN &lowerThan);

      void preSortTriangles(const VertCompFN &lowerThan);

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
                        const VertCompFN &lowerThan) const;

      bool compareLinks(const linkEdge &l0,
                        const linkEdge &l1,
                        const VertCompFN &lowerTan) const;

      // tools

      std::string printEdges(void) const;

      std::string printEdge(const idEdge e) const;
    };

    template <typename triangulationType>
    void Mesh<triangulationType>::preSortEdges(const VertCompFN &lowerThan) {
// Can't be parallel on a vector of bool !!
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(static)
#endif
      for(idEdge ei = 0; ei < nbEdges_; ++ei) {
        idVertex v0, v1;
        tri_->getEdgeVertex(ei, 0, v0);
        tri_->getEdgeVertex(ei, 1, v1);

        edgesSortId_[ei] = lowerThan(v0, v1);
      }
    }

    template <typename triangulationType>
    void
      Mesh<triangulationType>::preSortTriangles(const VertCompFN &lowerThan) {
// require edges to be already sorted
//
// Not sure we can parallelie on a vector of bitfields
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(static)
#endif
      for(idCell ti = 0; ti < nbTriangles_; ++ti) {
        idEdge e0, e1, e2;
        getTriangleEdge(ti, 0, e0);
        getTriangleEdge(ti, 1, e1);
        getTriangleEdge(ti, 2, e2);

        if(compareEdges(e0, e1, lowerThan)) {
          // 1 2 3
          // 1 3 2
          // 2 3 1
          if(compareEdges(e1, e2, lowerThan)) {
            // 0 : 1 2 3
            trianglesSortId_[ti] = 0;
          } else if(compareEdges(e0, e2, lowerThan)) {
            // 1 : 1 3 2
            trianglesSortId_[ti] = 1;
          } else {
            // 3 : 2 3 1
            trianglesSortId_[ti] = 3;
          }
        } else {
          // 2 1 3
          // 3 1 2
          // 3 2 1
          if(compareEdges(e0, e2, lowerThan)) {
            // 2 : 2 1 3
            trianglesSortId_[ti] = 2;
          } else if(compareEdges(e1, e2, lowerThan)) {
            // 4 : 3 1 2
            trianglesSortId_[ti] = 4;
          } else {
            // 5 : 3 2 1
            trianglesSortId_[ti] = 5;
          }
        }
      }
    }

    template <typename triangulationType>
    orderedEdge Mesh<triangulationType>::getOrderedEdge(
      const idEdge e, const bool increasingOrder) const {
      idVertex v0, v1;
      tri_->getEdgeVertex(e, 0, v0);
      tri_->getEdgeVertex(e, 1, v1);
      if(edgesSortId_[e] == static_cast<edgeScheme>(increasingOrder)) {
        return std::make_tuple(v0, v1);
      } else {
        return std::make_tuple(v1, v0);
      }
    }

    template <typename triangulationType>
    void Mesh<triangulationType>::getOrderedEdge(const idEdge e,
                                                 const bool increasingOrder,
                                                 orderedEdge &oEdge) const {
      idVertex v0, v1;
      tri_->getEdgeVertex(e, 0, v0);
      tri_->getEdgeVertex(e, 1, v1);
      if(edgesSortId_[e] == static_cast<edgeScheme>(increasingOrder)) {
        std::get<0>(oEdge) = v0;
        std::get<1>(oEdge) = v1;
      } else {
        std::get<0>(oEdge) = v1;
        std::get<1>(oEdge) = v0;
      }
    }

    template <typename triangulationType>
    orderedTriangle Mesh<triangulationType>::getOrderedTriangle(
      const idCell ti, const bool increasingOrder) const {
      idEdge e0, e1, e2;
      getTriangleEdge(ti, 0, e0);
      getTriangleEdge(ti, 1, e1);
      getTriangleEdge(ti, 2, e2);

      switch(trianglesSortId_[ti]) {
        case 0:
          if(increasingOrder)
            return std::make_tuple(e0, e1, e2);
          else
            return std::make_tuple(e2, e1, e0);

        case 1:
          if(increasingOrder)
            return std::make_tuple(e0, e2, e1);
          else
            return std::make_tuple(e1, e2, e0);

        case 2:
          if(increasingOrder)
            return std::make_tuple(e1, e0, e2);
          else
            return std::make_tuple(e2, e0, e1);

        case 3:
          if(increasingOrder)
            return std::make_tuple(e2, e0, e1);
          else
            return std::make_tuple(e1, e0, e2);

        case 4:
          if(increasingOrder)
            return std::make_tuple(e1, e2, e0);
          else
            return std::make_tuple(e0, e2, e1);

        case 5:
          if(increasingOrder)
            return std::make_tuple(e2, e1, e0);
          else
            return std::make_tuple(e0, e1, e2);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      std::cerr << "[Mesh]: error, unknownk case in triangleSortId_ : "
                << static_cast<unsigned>(trianglesSortId_[ti]) << " at " << ti
                << std::endl;
#endif

      return std::make_tuple(nullEdge, nullEdge, nullEdge);
    }

    template <typename triangulationType>
    void Mesh<triangulationType>::getOrderedTriangle(
      const idCell ti,
      const bool increasingOrder,
      orderedTriangle &oTriangle) const {
      idEdge e0, e1, e2;
      getTriangleEdge(ti, 0, e0);
      getTriangleEdge(ti, 1, e1);
      getTriangleEdge(ti, 2, e2);

      switch(trianglesSortId_[ti]) {
        case 0:
          if(increasingOrder) {
            std::get<0>(oTriangle) = e0;
            std::get<1>(oTriangle) = e1;
            std::get<2>(oTriangle) = e2;
          } else {
            std::get<0>(oTriangle) = e2;
            std::get<1>(oTriangle) = e1;
            std::get<2>(oTriangle) = e0;
          }
          break;

        case 1:
          if(increasingOrder) {
            std::get<0>(oTriangle) = e0;
            std::get<1>(oTriangle) = e2;
            std::get<2>(oTriangle) = e1;
          } else {
            std::get<0>(oTriangle) = e1;
            std::get<1>(oTriangle) = e2;
            std::get<2>(oTriangle) = e0;
          }
          break;

        case 2:
          if(increasingOrder) {
            std::get<0>(oTriangle) = e1;
            std::get<1>(oTriangle) = e0;
            std::get<2>(oTriangle) = e2;
          } else {
            std::get<0>(oTriangle) = e2;
            std::get<1>(oTriangle) = e0;
            std::get<2>(oTriangle) = e1;
          }
          break;

        case 3:
          if(increasingOrder) {
            std::get<0>(oTriangle) = e2;
            std::get<1>(oTriangle) = e0;
            std::get<2>(oTriangle) = e1;
          } else {
            std::get<0>(oTriangle) = e1;
            std::get<1>(oTriangle) = e0;
            std::get<2>(oTriangle) = e2;
          }
          break;

        case 4:
          if(increasingOrder) {
            std::get<0>(oTriangle) = e1;
            std::get<1>(oTriangle) = e2;
            std::get<2>(oTriangle) = e0;
          } else {
            std::get<0>(oTriangle) = e0;
            std::get<1>(oTriangle) = e2;
            std::get<2>(oTriangle) = e1;
          }
          break;

        case 5:
          if(increasingOrder) {
            std::get<0>(oTriangle) = e2;
            std::get<1>(oTriangle) = e1;
            std::get<2>(oTriangle) = e0;
          } else {
            std::get<0>(oTriangle) = e0;
            std::get<1>(oTriangle) = e1;
            std::get<2>(oTriangle) = e2;
          }
          break;
        default:
#ifndef TTK_ENABLE_KAMIKAZE
          std::cerr << "[Mesh]: error, unknownk case in triangleSortId_ : "
                    << static_cast<unsigned>(trianglesSortId_[ti]) << " at "
                    << ti << std::endl;
#endif
          break;
      }
    }

    // tools

    template <typename triangulationType>
    bool Mesh<triangulationType>::compareEdges(
      const idEdge e0, const idEdge e1, const VertCompFN &lowerThan) const {
      idVertex e0v0, e0v1;
      if(edgesSortId_[e0] == 0) {
        getEdgeVertex(e0, 0, e0v0);
        getEdgeVertex(e0, 1, e0v1);
      } else {
        getEdgeVertex(e0, 1, e0v0);
        getEdgeVertex(e0, 0, e0v1);
      }

      idVertex e1v0, e1v1;
      if(edgesSortId_[e1] == 0) {
        getEdgeVertex(e1, 0, e1v0);
        getEdgeVertex(e1, 1, e1v1);
      } else {
        getEdgeVertex(e1, 1, e1v0);
        getEdgeVertex(e1, 0, e1v1);
      }

      if(e0v0 == e1v0)
        return lowerThan(e0v1, e1v1);
      return lowerThan(e0v0, e1v0);
    }

    template <typename triangulationType>
    bool
      Mesh<triangulationType>::compareLinks(const linkEdge &l0,
                                            const linkEdge &l1,
                                            const VertCompFN &lowerThan) const {
      if(std::get<0>(l0) == std::get<0>(l1)) {
        return compareEdges(std::get<1>(l0), std::get<1>(l1), lowerThan);
      }
      return compareEdges(std::get<0>(l0), std::get<0>(l1), lowerThan);
    }

    template <typename triangulationType>
    std::string Mesh<triangulationType>::printEdges(void) const {
      std::stringstream res;
      res << " edges: " << std::endl;
      for(idEdge i = 0; i < nbEdges_; ++i) {
        idVertex v0, v1;
        getEdgeVertex(i, 0, v0);
        getEdgeVertex(i, 1, v1);
        res << i << " : " << v0 << " " << v1 << std::endl;
      }

      return res.str();
    }

    template <typename triangulationType>
    std::string Mesh<triangulationType>::printEdge(const idEdge e) const {
      std::stringstream res;
      idVertex v0, v1;
      getEdgeVertex(e, 0, v0);
      getEdgeVertex(e, 1, v1);
      res << v0 << " " << v1;
      return res.str();
    }

  } // namespace ftr
} // namespace ttk
