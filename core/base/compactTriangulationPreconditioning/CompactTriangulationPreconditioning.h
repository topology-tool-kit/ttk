/// \ingroup base
/// \class ttk::CompactTriangulationPreconditioning
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date May 2021.
///
/// \brief TTK processing package for mesh preprocessing before using
/// TopoCluster.
///
/// Given a simplicial mesh, this class uses PR star octree to divide it into
/// multiple clusters, and the cluster of each vertex is written as a new scalar
/// field.
///
/// \b Related \b publications \n
/// "The PR-star octree: A spatio-topological data structure for tetrahedral
/// meshes." Kenneth Weiss, Leila Floriani, Riccardo Fellegara, and Marcelo
/// Velloso In Proceedings of the 19th ACM SIGSPATIAL International Conference
/// on Advances in Geographic Information Systems, 2011.
///
/// "TopoCluster: A Localized Data Structure for Topology-based Visualization"
/// Guoxi Liu, Federico Iuricich, Riccardo Fellegara, and Leila De Floriani
/// IEEE Transactions on Visualization and Computer Graphics, 2021.
///
/// \sa ttk::TopoCluster
/// \sa ttkCompactTriangulationPreconditioning.cpp %for a usage example.

#pragma once

// ttk common includes
#include <Debug.h>
#include <Octree.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The CompactTriangulationPreconditioning class provides methods to compute
   * the indexing for each vertex and then used for TopoCluster data structure.
   */
  class CompactTriangulationPreconditioning : virtual public Debug {

  public:
    CompactTriangulationPreconditioning() {
      this->setDebugMsgPrefix(
        "CompactTriangulationPreconditioning"); // inherited from Debug: prefix
                                                // will be printed at the
      // beginning of every msg
    };
    ~CompactTriangulationPreconditioning(){};

    template <class triangulationType = ttk::AbstractTriangulation>
    int execute(const triangulationType *triangulation,
                const int &argument) const {
      // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(triangulation->getNumberOfVertices())},
      });
      this->printMsg(ttk::debug::Separator::L1);

      // -----------------------------------------------------------------------
      // Compute Octree
      // -----------------------------------------------------------------------
      {
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

#ifndef TTK_ENABLE_KAMIKAZE
        if(!triangulation)
          return -1;
        if(!vertices)
          return -2;
        if(!nodes)
          return -3;
        if(!cells)
          return -4;
#endif

        SimplexId vertexNumber = triangulation->getNumberOfVertices();
        SimplexId cellNumber = triangulation->getNumberOfCells();

        // create the octree
        Octree preOctree(triangulation, argument);
        for(SimplexId i = 0; i < vertexNumber; i++) {
          preOctree.insertVertex(i);
        }

        for(SimplexId i = 0; i < cellNumber; i++) {
          preOctree.insertCell(i);
        }

        if(preOctree.verifyTree(vertexNumber)) {
          cerr << "[PreprocessStellar] The construction of the tree failed!\n";
          return -1;
        }

        preOctree.reindex(vertices, nodes, cells);
        this->printMsg({
          {"Size of vertex vector", std::to_string(vertices->size())},
          {"Size of cell vector", std::to_string(cells->size())},
        });
      }

      // ---------------------------------------------------------------------
      // print global performance
      // ---------------------------------------------------------------------
      {
        this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
        this->printMsg(
          "Complete", 1, globalTimer.getElapsedTime() // global progress, time
        );
        this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
      }

      return 1; // return success
    }

  protected:
    std::vector<SimplexId> *vertices, *nodes, *cells;

  }; // CompactTriangulationPreconditioning class

} // namespace ttk