/// \ingroup base
/// \class ttk::RipsPersistenceDiagram
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.
///
/// \brief TTK base class that computes the persistence diagram of a Rips
/// complex.
///
/// This module defines the %RipsPersistenceDiagram class that takes a point
/// cloud or a distance matrix and computes the persistence diagram of its Rips
/// complex.
///
/// \sa ttkRipsPersistenceDiagram.cpp %for a usage example.

#pragma once

// ttk common includes
#include <Debug.h>

#include <FastRipsPersistenceDiagram2.h>
#include <ripser.h>

namespace ttk {

  /**
   * The RipsPersistenceDiagram class provides a method to call the code Ripser
   * in order to compute the persistence diagram of the Rips complex of the
   * input.
   */
  class RipsPersistenceDiagram : virtual public Debug {

  public:
    enum class BACKEND : std::uint8_t {
      RIPSER = 0,
      GEOMETRY = 1,
    };

    RipsPersistenceDiagram();

    /**
     * @brief Main entry point
     *
     * @param[in] points Input point cloud or input distance matrix
     * @param[out] ph Persistence diagram
     */
    int execute(const rpd::PointCloud &points,
                rpd::MultidimensionalDiagram &ph) const;

  protected:
    /** BackEnd */
    BACKEND BackEnd{BACKEND::RIPSER};
    /** Max dimension of computed persistence diagram */
    int HomologyMaximumDimension{1};
    /** Rips diameter threshold */
    double SimplexMaximumDiameter{rpd::inf};
    /** Field of coefficients */
    int FieldOfCoefficients{2};
    /** is input a distance matrix */
    bool InputIsDistanceMatrix{false};

  }; // RipsPersistenceDiagram class

} // namespace ttk
