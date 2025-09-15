/// \ingroup base
/// \class ttk::DelaunayRipsPersistenceDiagram
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date July 2025.
///
/// \brief TTK base class that computes the persistence diagram of a
/// Delaunay-Rips filtration.
///
/// This module defines the %DelaunayRipsPersistenceDiagram class that takes a
/// point cloud and computes the persistence diagram of its Delaunay-Rips
/// filtration. It can also compute 1-dimensional and 2-dimensional persistence
/// generators for point clouds in R2 or R3.
///
/// \sa ttkDelaunayRipsPersistenceDiagram.cpp %for a usage example.

#pragma once

// ttk common includes
#include <Debug.h>

#include <FastRipsPersistenceDiagram2.h>
#include <geoPH3.h>
#include <geoPHd.h>

namespace ttk {

  /**
   * The DelaunayRipsPersistenceDiagram class provides a method to call the
   * relevant code in order to compute the persistence diagram of the Delaunay-
   * Rips filtration of the input point cloud.
   */
  class DelaunayRipsPersistenceDiagram : virtual public Debug {

  public:

    DelaunayRipsPersistenceDiagram();

    /**
    * @brief Main entry point (no generators)
    *
    * @param[in] points Input point cloud
    * @param[out] ph Persistence diagram
    */
    int execute(const PointCloud &points,
                MultidimensionalDiagram &ph) const;

    /**
    * @brief Main entry point
    *
    * @param[in] points Input point cloud
    * @param[out] ph Persistence diagram
    * @param[out] generators1 1-dimensional persistent generators, if required
    * @param[out] generators2 2-dimensional persistent generators, if required
    */
    int execute(const PointCloud &points,
                MultidimensionalDiagram &ph,
                std::vector<Generator1> &generators1,
                std::vector<Generator2> &generators2) const;

  }; // DelaunayRipsPersistenceDiagram class

} // namespace ttk
