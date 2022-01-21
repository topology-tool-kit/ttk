/// \ingroup base
/// \class ttk::RipsComplex
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date January 2022
///
/// \brief TTK VTK-filter that computes a Rips complex.
///
/// \param Input Input table (vtkTable)
/// \param Output Triangulation (vtkUnstructuredGrid)
///
/// \brief TTK VTK-filter that takes a matrix (vtkTable) as input and
/// computes a Rips complex from it to generate an explicit
/// triangulation.
///
/// \sa ttk::Triangulation
/// \sa ttkRipsComplex.cpp %for a usage example.

#pragma once

#include <Debug.h>

namespace ttk {

  class RipsComplex : virtual public Debug {
  public:
    RipsComplex();

    int execute(std::vector<SimplexId> &connectivity,
                std::vector<double> &diameters,
                const std::vector<std::vector<double>> &inputMatrix) const;

  protected:
    int computeDistanceMatrix(
      std::vector<std::vector<double>> &distanceMatrix,
      const std::vector<std::vector<double>> &inputMatrix) const;

    /** Dimension of the generated complex */
    int OutputDimension{2};
    /** Distance threshold */
    double Epsilon{1.0};
    /** If input matrix is a distance matrix (compute it otherwise) */
    bool InputIsADistanceMatrix{false};
  };

} // namespace ttk
