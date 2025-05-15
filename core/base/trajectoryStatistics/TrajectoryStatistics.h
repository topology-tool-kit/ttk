/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::TrajectoryStatistics
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %TrajectoryStatistics class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'TrajectoryStatistics'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The TrajectoryStatistics class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class TrajectoryStatistics : virtual public Debug {

  public:
    TrajectoryStatistics();

    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      triangulation->preconditionVertexNeighbors();
      return triangulation->preconditionVertexNeighbors();
    }

     int findSurface(
                    ttk::SimplexId                          vertexId,
                    std::vector<ttk::SimplexId>            &surfVertex,
                    const std::vector<std::vector<double>> &vertexScalars,
                    std::vector<char>                      &visited,
                    const double                            threshold,
                    int                                     frame,
                    const ttk::AbstractTriangulation       *triangulation
    ); 

    /**
     * TODO 3: Implementation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */
    int execute(std::vector<std::vector<int>> &trajTime,         //input
                std::vector<std::vector<double>> &trajX,         //input 
                std::vector<std::vector<double>> &trajY,         //input
                std::vector<std::vector<double>> &trajZ,         //input 
                std::vector<std::vector<int>>    &trajVertexId,   //input
                std::vector<std::vector<double>> &vertexScalars,
                std::vector<int> &startFrames,          //output
                std::vector<int> &endFrames,            //output
                std::vector<int> &durations,            //output
                std::vector<double> &VX,
                std::vector<double> &VY,
                std::vector<double> &surfMin,
                std::vector<double> &surfMax,
                std::vector<double> &surfMoy,
                std::vector<std::vector<ttk::SimplexId>> &allVertexDebris,
                int frameSurf,
                ttk::AbstractTriangulation*triangulation);

  }; // TrajectoryStatistics class

} // namespace ttk
