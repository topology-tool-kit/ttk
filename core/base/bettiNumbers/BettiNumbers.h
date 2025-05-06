/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::BettiNumbers
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %BettiNumbers class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'BettiNumbers'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <UnionFind.h>

namespace ttk {

  /**
   * The BettiNumbers class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class BettiNumbers : virtual public Debug {

  public:
    BettiNumbers();

    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation)  {
      triangulation->preconditionVertexNeighbors();
      triangulation->preconditionEdges();
      
      triangulation_ = triangulation;
      
      return 0;
    }

    int execute();

    int getB0() { return B0; };

  private:
    int B0 = -1;
    ttk::AbstractTriangulation *triangulation_;

  }; // BettiNumbers class

} // namespace ttk
