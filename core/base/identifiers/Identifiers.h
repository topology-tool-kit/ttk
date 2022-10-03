/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::Identifiers
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %Identifiers class that computes for each vertex of
/// a triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'Identifiers'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Geometry.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The Identifiers class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class Identifiers : virtual public Debug {

  public:
    Identifiers();

    int setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }

    inline void setDomainDimension(const int &dimension) {
      dimension_ = dimension;
    }

    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    int preconditionTriangulation(ttk::AbstractTriangulation *triangulation) {
      setDomainDimension(triangulation->getDimensionality());
      setVertexNumber(triangulation->getNumberOfVertices());
      return 0;
    }

    /**
     * TODO 3: Implmentation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */
    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int execute(dataType *outputData,
                const dataType *inputData,
                const triangulationType *triangulation) const {
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
      // Compute Vertex Averages
      // -----------------------------------------------------------------------
      {
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

        // print the progress of the current subprocedure (currently 0%)
        this->printMsg("Computing Averages",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_, ttk::debug::LineMode::REPLACE);

        // print the progress of the current subprocedure with elapsed time
        this->printMsg("Computing Averages",
                       1, // progress
                       localTimer.getElapsedTime(), this->threadNumber_);
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
    int nbPoint_{0};
    std::map<ttk::SimplexId, ttk::SimplexId> vertGtoL_;
    std::vector<int> neighbors_;
    double bounds_[6];
    MPI_Datatype mpiResponseType_;
    ttk::SimplexId vertexNumber_{};
    int dimension_{};

    template <class triangulationType>
    void inline findPoint(ttk::SimplexId &id,
                          float x,
                          float y,
                          float z,
                          triangulationType *triangulation) {
      float pointToFind[3] = {x, y, z};
      float dist{0};
      float p[3];
      triangulation->getVertexPoint(0, p[0], p[1], p[2]);
      float minDist = Geometry::distance(pointToFind, p);
      ttk::SimplexId indexMin = 0;
      for(int i = 1; i < vertexNumber_; i++) {
        triangulation->getVertexPoint(i, p[0], p[1], p[2]);
        dist = Geometry::distance(pointToFind, p, dimension_);
        if(dist < minDist) {
          minDist = dist;
          indexMin = i;
        }
      }
      id = indexMin;
    }

  }; // Identifiers class

} // namespace ttk
