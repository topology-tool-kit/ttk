#pragma once

#include "FTRGraph.h"

namespace ttk {
  namespace ftr {

    template <typename ScalarType, typename triangulationType>
    std::string FTRGraph<ScalarType, triangulationType>::printMesh(void) const {
      std::stringstream res;

      res << "cells     : " << mesh_.getNumberOfCells() << std::endl;
      res << "triangles : " << mesh_.getNumberOfTriangles() << std::endl;
      res << "edges     : " << mesh_.getNumberOfEdges() << std::endl;
      res << "vertices  : " << mesh_.getNumberOfVertices() << std::endl;

      return res.str();
    }

    template <typename ScalarType, typename triangulationType>
    std::string FTRGraph<ScalarType, triangulationType>::printEdge(
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

    template <typename ScalarType, typename triangulationType>
    std::string FTRGraph<ScalarType, triangulationType>::printTriangle(
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

    template <typename ScalarType, typename triangulationType>
    void FTRGraph<ScalarType, triangulationType>::printGraph(
      const int verbosity) const {
      std::cout << graph_.print(verbosity) << std::endl;
    }

    template <typename ScalarType, typename triangulationType>
    void FTRGraph<ScalarType, triangulationType>::printTime(
      Timer &timer, const std::string &msg) const {
      this->printMsg(msg, 1.0, timer.getElapsedTime(), this->threadNumber_);
    }

  } // namespace ftr
} // namespace ttk
