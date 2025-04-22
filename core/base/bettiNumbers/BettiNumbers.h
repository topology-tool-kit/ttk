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
      triangulation_ = triangulation;
      return triangulation->preconditionEdges();
    }

    /**
     * TODO 3: Implementation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */
    int execute();

    int getB0() const { return B0_; }
    
  private:
    ttk::AbstractTriangulation *triangulation_ = nullptr;
    int B0_ = -1;
    int B1_ = -1;
    int B2_ = -1;

  }; // BettiNumbers class

  } // namespace ttk
