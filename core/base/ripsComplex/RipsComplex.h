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

    /**
     * @brief Main entry point
     *
     * @param[out] connectivity Cell connectivity array (VTK format)
     * @param[out] diameters Cell diameters
     * @param[out] diamStats Min, mean and max cell diameters around point
     * @param[in] inputMatrix Either coordinates array or distance matrix
     * @param[out] density Gaussian density array on points
     */
    int execute(std::vector<SimplexId> &connectivity,
                std::vector<double> &diameters,
                std::array<double *const, 3> diamStats,
                const std::vector<std::vector<double>> &inputMatrix,
                double *const density = nullptr) const;

  protected:
    /**
     * @brief Compute distance matrix from coordinates matrix
     *
     * @param[out] distanceMatrix Output distance matrix
     * @param[in] inputMatrix Input coordinates matrix
     */
    int computeDistanceMatrix(
      std::vector<std::vector<double>> &distanceMatrix,
      const std::vector<std::vector<double>> &inputMatrix) const;

    /**
     * @brief Compute diameter statistics on points
     *
     * @param[in] nPoints Number of input points
     * @param[out] diamStats Min, mean and max cell diameters around point
     * @param[in] connectivity Cell connectivity array (pre-filled)
     * @param[in] cellDiameters Cell diameters
     */
    int computeDiameterStats(const SimplexId nPoints,
                             std::array<double *const, 3> diamStats,
                             const std::vector<SimplexId> &connectivity,
                             const std::vector<double> &cellDiameters) const;

    /**
     * @brief Compute Gaussian density on points
     *
     * @param[out] density Gaussian density array on points
     * @param[in] distanceMatrix Distance matrix between points
     */
    int computeGaussianDensity(
      double *const density,
      const std::vector<std::vector<double>> &distanceMatrix) const;

    /** Dimension of the generated complex */
    int OutputDimension{2};
    /** Distance threshold */
    double Epsilon{1.0};
    /** Standard Deviation for Gaussian density */
    double StdDev{1.0};
    /** Compute the Gaussian density from the distance matrix */
    bool ComputeGaussianDensity{false};
    /** If input matrix is a distance matrix (compute it otherwise) */
    bool InputIsADistanceMatrix{false};
  };

} // namespace ttk
