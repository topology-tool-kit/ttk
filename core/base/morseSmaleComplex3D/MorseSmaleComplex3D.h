/// \ingroup base
/// \class ttk::MorseSmaleComplex3D::
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK %morseSmaleComplex3D processing package.
///
/// %MorseSmaleComplex3D is a TTK processing package that takes a scalar field
/// on the input and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkMorseSmaleComplex3D.cpp %for a usage example.

#pragma once

// base code includes
#include <AbstractMorseSmaleComplex.h>

#include <set>

namespace ttk {

  /**
   * Class specialized in building the Morse-Smale complex
   * of 3D triangulation.
   */
  class MorseSmaleComplex3D : public AbstractMorseSmaleComplex {

  public:
    MorseSmaleComplex3D();
    ~MorseSmaleComplex3D();

    /**
     * Main function for computing the whole Morse-Smale complex.
     */
    template <typename dataType, typename idtype, typename triangulationType>
    int execute(const triangulationType &triangulation);

    /**
     * Compute the (saddle1, saddle2) pairs not detected by the
     * contour tree.
     */
    template <typename dataType, typename idType, typename triangulationType>
    int computePersistencePairs(
      std::vector<std::tuple<SimplexId, SimplexId, dataType>>
        &pl_saddleSaddlePairs,
      const triangulationType &triangulation);

    template <typename dataType>
    int setAugmentedCriticalPoints(const std::vector<dcg::Cell> &criticalPoints,
                                   SimplexId *ascendingManifold,
                                   SimplexId *descendingManifold) const;

    /**
     * Compute the descending 1-separatrices by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getAscendingSeparatrices1(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;

    /**
     * Compute the saddle-connectors by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getSaddleConnectors(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;

    /**
     * Compute the 2-separatrices by reading into the discrete
     * gradient from the maxima.
     */
    template <typename triangulationType>
    int getDescendingSeparatrices2(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * Compute the geometrical embedding of the descending
     * 2-separatrices. This function needs the following
     * internal pointers to be set:
     * outputSeparatrices2_numberOfPoints_
     * outputSeparatrices2_points_
     * outputSeparatrices2_numberOfCells_
     * outputSeparatrices2_cells_
     * inputScalarField_
     */
    template <typename dataType, typename triangulationType>
    int setDescendingSeparatrices2(
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * Find all tetras in the star of edgeId
     *
     * (primal: star of edgeId -> dual: vertices of polygon)
     */
    template <typename triangulationType>
    int getDualPolygon(const SimplexId edgeId,
                       std::vector<SimplexId> &polygon,
                       const triangulationType &triangulation) const;

    /**
     * Sort the polygon vertices to be clockwise
     */
    template <typename triangulationType>
    int sortDualPolygonVertices(std::vector<SimplexId> &polygon,
                                const triangulationType &triangulation) const;

    /**
     * Compute the 2-separatrices by reading into the discrete
     * gradient from the minima.
     */
    template <typename triangulationType>
    int getAscendingSeparatrices2(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * Compute the geometrical embedding of the ascending
     * 2-separatrices. This function needs the following
     * internal pointers to be set:
     * outputSeparatrices2_numberOfPoints_
     * outputSeparatrices2_points_
     * outputSeparatrices2_numberOfCells_
     * outputSeparatrices2_cells_
     * inputScalarField_
     */
    template <typename dataType, typename triangulationType>
    int setAscendingSeparatrices2(
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * @brief Flatten the vectors of vectors into their first component
     */
    void flattenSeparatricesVectors(
      std::vector<std::vector<ttk::Separatrix>> &separatrices,
      std::vector<std::vector<std::vector<ttk::dcg::Cell>>>
        &separatricesGeometry) const;
  };
} // namespace ttk

template <typename dataType, typename idType, typename triangulationType>
int ttk::MorseSmaleComplex3D::computePersistencePairs(
  std::vector<std::tuple<SimplexId, SimplexId, dataType>> &pl_saddleSaddlePairs,
  const triangulationType &triangulation) {

  const dataType *scalars = static_cast<const dataType *>(inputScalarField_);

  std::vector<std::array<dcg::Cell, 2>> dmt_pairs;
  {
    // simplify to be PL-conformant
    discreteGradient_.setDebugLevel(debugLevel_);
    discreteGradient_.setThreadNumber(threadNumber_);
    discreteGradient_.setCollectPersistencePairs(false);
    discreteGradient_.buildGradient<dataType, idType>(triangulation);
    discreteGradient_.reverseGradient<dataType, idType>(triangulation);

    // collect saddle-saddle connections
    discreteGradient_.setCollectPersistencePairs(true);
    discreteGradient_.setOutputPersistencePairs(&dmt_pairs);
    discreteGradient_.reverseGradient<dataType, idType>(triangulation, false);
  }

  // transform DMT pairs into PL pairs
  for(const auto &pair : dmt_pairs) {
    const SimplexId v0
      = discreteGradient_.getCellGreaterVertex(pair[0], triangulation);
    const SimplexId v1
      = discreteGradient_.getCellGreaterVertex(pair[1], triangulation);
    const dataType persistence = scalars[v1] - scalars[v0];

    if(v0 != -1 and v1 != -1 and persistence >= 0) {
      if(!triangulation.isVertexOnBoundary(v0)
         or !triangulation.isVertexOnBoundary(v1)) {
        pl_saddleSaddlePairs.emplace_back(v0, v1, persistence);
      }
    }
  }
  return 0;
}
