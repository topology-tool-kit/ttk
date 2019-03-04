#ifndef FTRGRAPHPRINT_TEMPLATE_H
#define FTRGRAPHPRINT_TEMPLATE_H

#include "FTRGraph.h"

namespace ttk {
  namespace ftr {

    template <typename ScalarType>
    std::string FTRGraph<ScalarType>::printMesh(void) const {
      std::stringstream res;

      res << "cells     : " << mesh_.getNumberOfCells() << std::endl;
      res << "triangles : " << mesh_.getNumberOfTriangles() << std::endl;
      res << "edges     : " << mesh_.getNumberOfEdges() << std::endl;
      res << "vertices  : " << mesh_.getNumberOfVertices() << std::endl;

      return res.str();
    }

    template <typename ScalarType>
    std::string FTRGraph<ScalarType>::printEdge(
      const idEdge edgeId, const Propagation *const localPropagation) const {
      const orderedEdge oEdge
        = mesh_.getOrderedEdge(edgeId, localPropagation->goUp());
      std::stringstream res;
      res << "e" << edgeId << ":";
      res << "(";
      res << std::get<0>(oEdge) << " - " << std::get<1>(oEdge);
      res << ")";
      return res.str();
    }

    template <typename ScalarType>
    std::string FTRGraph<ScalarType>::printTriangle(
      const idCell cellId, const Propagation *const localPropagation) const {
      std::stringstream res;
      const orderedTriangle oTriangle
        = mesh_.getOrderedTriangle(cellId, localPropagation->goUp());
      orderedEdge e0, e1, e2;

      e0 = mesh_.getOrderedEdge(
        std::get<0>(oTriangle), localPropagation->goUp());
      e1 = mesh_.getOrderedEdge(
        std::get<1>(oTriangle), localPropagation->goUp());
      e2 = mesh_.getOrderedEdge(
        std::get<2>(oTriangle), localPropagation->goUp());

      res << "t" << cellId << ":";
      res << "{";
      res << " " << printEdge(std::get<0>(oTriangle), localPropagation);
      res << " " << printEdge(std::get<1>(oTriangle), localPropagation);
      res << " " << printEdge(std::get<2>(oTriangle), localPropagation);
      res << " }";
      return res.str();
    }

    template <typename ScalarType>
    void FTRGraph<ScalarType>::printGraph(const int verbosity) const {
      std::cout << graph_.print(verbosity) << std::endl;
    }

    template <typename ScalarType>
    void FTRGraph<ScalarType>::printTime(DebugTimer &timer,
                                         const std::string &msg,
                                         const int lvl) const {
      std::ostringstream outString(std::string(lvl, ' '));
      outString << msg << timer.getElapsedTime() << std::endl;
      dMsg(std::cout, outString.str(), lvl);
    }

  } // namespace ftr
} // namespace ttk

#endif /* end of include guard: FTRGRAPHPRINT_TEMPLATE_H */
