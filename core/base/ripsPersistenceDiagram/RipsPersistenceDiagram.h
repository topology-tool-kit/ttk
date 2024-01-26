/// \ingroup base
/// \class ttk::RipsPersistenceDiagram
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.
///
/// This module defines the %RipsPersistenceDiagram class that computes the persistence diagram
/// of the Rips complex of the input point cloud.

/// \param Input Input table (vtkTable)
/// \param Output PersistenceDiagram (vtkUnstructuredGrid)
///

#pragma once

// ttk common includes
#include <Debug.h>

#include <ripser.h>

namespace ttk {

  /**
   * The RipsPersistenceDiagram class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class RipsPersistenceDiagram : virtual public Debug {

  public:
    RipsPersistenceDiagram();

    /**
     * @brief Main entry point
     *
     * @param[in] points Input point cloud in any dimension
     * @param[out] ph Computed Rips persistence diagram
     */
    int execute(const std::vector<std::vector<double>> &points, std::vector<std::vector<pers_pair_t> >& ph) const;

  protected:
    /** Max dimension of computed persistence diagram */
    int MaxDim{1};
    /** Rips threshold */
    double Threshold{1.0};
    /** is input a distance matrix */
    int InputIsDistanceMatrix{0};

  }; // RipsPersistenceDiagram class

} // namespace ttk
