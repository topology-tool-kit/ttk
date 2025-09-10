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
/// filtration.
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
    * @brief Main entry point
    *
    * @param[in] points Input point cloud
    * @param[out] ph Persistence diagram
    * @param[out] generators Persistent generators, if required
    */
    int execute(const rpd::PointCloud &points,
                rpd::MultidimensionalDiagram &ph,
                std::vector<rpd::Generator> &generators) const;

  protected:
    /** output generators */
    bool OutputGenerators{false};

  }; // DelaunayRipsPersistenceDiagram class

} // namespace ttk
