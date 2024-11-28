/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::MorseSmallComplexStability
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %MorseSmallComplexStability class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'MorseSmallComplexStability'
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
   * The MorseSmallComplexStability class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class MorseSmallComplexStability : virtual public Debug {

  public:
    MorseSmallComplexStability();

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

  }; // MorseSmallComplexStability class

} // namespace ttk
