/// \ingroup base
/// \class ttk::NumericalIntegralLines
/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
/// \date May 2026
/// \date NumericalIntegralLines extractor wrapping the DiscreteGradient class.
///
/// \brief TTK convenience class wrapping the DiscreteGradient class for
/// the easy extraction of vpaths.
///
/// Given a simplexId and dimension, this class returns a descending (or
/// ascending) vpath started in the given input simplex.
///
/// \sa NumericalIntegralLines.cpp %for an alternative integral line backend.
/// \sa DiscreteGradient.cpp %for the core mechanisms.
/// \sa ttkNumericalIntegralLines.cpp %for a usage example.
///

#pragma once

// base code includes
#include <Triangulation.h>
// std includes

namespace ttk {
  namespace nil {

    struct PathPoint{
      SimplexId             simplexId_;
      int                   simplexDimension_;
      std::vector<float>    barycentricWeights_;
    };

    class NumericalIntegralLines : virtual public Debug {

    public:
      NumericalIntegralLines();
      ~NumericalIntegralLines() override;

      template <class dataType, class triangulationType>
        int computeEndPoint(const triangulationType *triangulation,
          const PathPoint &start, PathPoint &end) const;

      /**
       * @brief Compute a single numerical integral line.
       *
       * @param seed (SimplexId, dimension)
       * @param barycentricWeights Weights for the input seed.
       * @param output Output integral line (vector of 3D points).
       * @param isForawrd Forward or backward line (default: forward).
       */
      template <class dataType, class triangulationType>
        int computeIntegralLine(const triangulationType *triangulation,
          const std::pair<SimplexId, int> &seed,
          const std::vector<float> &barycentricWeights,
          std::vector<PathPoint> &output,
          const bool &isForward = false) const;

      template <class dataType, class triangulationType>
        int computeNumericalGradient(const triangulationType *triangulation,
          const int &simplexDimension, const int &simplexId,
          std::vector<float> &gradient) const;

      /**
       * @brief Compute numerical integral lines.
       *
       * @param output Vector storing the output vpaths (1 entry per seed,
       * with possibly multiple v-path per seed).
       * @param isForward Forward or backward vpath (default: forward).
       */
      template <class dataType, class triangulationType>
      int execute(const triangulationType *triangulation,
        const std::vector<std::pair<SimplexId, int>> &seeds,
        std::vector<std::vector<PathPoint>> &output,
        const bool &isForward = false) const;

      template <class triangulationType>
        int getVertexIdentifiers(const triangulationType *triangulation,
          const int &simplexDimension, const int &simplexId,
          std::vector<SimplexId> vertexIdentifiers) const;

      /**
       * @brief Triangulation preconditioning.
       */
      inline void preconditionTriangulation(AbstractTriangulation *triangulation){
        // precondition simplex2face
        // precondition face2cofacets
      }

      inline void setInputScalarField(const void *const scalars){
        scalars_ = scalars;
      }

    protected:
      int maximumIterationNumber_{1000000000};
      const void *scalars_;
    };
  } // namespace nil
} // namespace ttk

template <class dataType, class triangulationType>
  int ttk::nil::NumericalIntegralLines::computeEndPoint(
    const triangulationType *triangulation,
    const ttk::nil::PathPoint &start, ttk::nil::PathPoint &end) const{

  std::vector<float> gradient(3);

  computeNumericalGradient<dataType, triangulationType>(
    triangulation, start.simplexDimension_, start.simplexId_, gradient);

  return 0;
}

template <class dataType, class triangulationType>
  int ttk::nil::NumericalIntegralLines::computeIntegralLine(
    const triangulationType *triangulation,
    const std::pair<SimplexId, int> &seed,
    const std::vector<float> &startBarycentricWeights,
    std::vector<PathPoint> &output,
    const bool &isForward) const{

  output.clear();

  PathPoint startPoint, endPoint;

  startPoint.simplexId_ = seed.first;
  startPoint.simplexDimension_ = seed.second;
  startPoint.barycentricWeights_ = startBarycentricWeights;

  for(int i = 0; i < (int) maximumIterationNumber_; i++){

    computeEndPoint<dataType, triangulationType>(triangulation, startPoint, endPoint);
  }

  return 0;
}

// TODO
// move that function to the geometry class

template <class dataType, class triangulationType>
  int ttk::nil::NumericalIntegralLines::computeNumericalGradient(
    const triangulationType *triangulation,
    const int &simplexDimension, const int &simplexId,
    std::vector<float> &gradient) const{

  gradient = {0, 0, 0};

  if(!simplexDimension)
    return -1;

  std::vector<SimplexId> vertexIdentifiers;

  getVertexIdentifiers(triangulation, simplexDimension, simplexId, vertexIdentifiers);

  const int vertexNumber = vertexIdentifiers.size();

  std::vector<std::array<float, 3>> vertexPoints(vertexNumber);
  std::vector<dataType> vertexScalars(vertexNumber);

  for(int i = 0; i < (int) vertexNumber; i++){
    triangulation->getVertexPoint(vertexIdentifiers[i],
      vertexPoints[i][0], vertexPoints[i][1], vertexPoints[i][2]);
    vertexScalars[i] = ((dataType *) scalars_)[vertexIdentifiers[i]];
  }

  // build edge vectors and corresponding differences, wrt v0
  std::vector<std::array<float, 3>> edgeVectors(simplexDimension);
  std::vector<dataType> edgeDifferences(simplexDimension);

  for(int i = 0; i < simplexDimension; i++) {
    for(int c = 0; c < 3; c++)
      edgeVectors[i][c] = vertexPoints[i + 1][c] - vertexPoints[0][c];
    edgeDifferences[i] = vertexScalars[i + 1] - vertexScalars[0];
  }

  // Gram matrix  gramMatrix[i][j] = edgeVectors[i] . edgeVectors[j]
  // (simplexDimension x simplexDimension)
  std::vector<std::vector<float>>
    gramMatrix(simplexDimension, std::vector<float>(simplexDimension, 0));

  for(int i = 0; i < simplexDimension; i++)
    for(int j = 0; j < simplexDimension; j++)
      gramMatrix[i][j] = ttk::Geometry::dotProduct(
        edgeVectors[i].data(), edgeVectors[j].data());

  // Gaussian elimintation
  std::vector<std::vector<float>>
    augmentedMatrix(simplexDimension, std::vector<float>(simplexDimension + 1));
  for(int i = 0; i < simplexDimension; ++i) {
    for(int j = 0; j < simplexDimension; ++j)
      augmentedMatrix[i][j] = gramMatrix[i][j];
    augmentedMatrix[i][simplexDimension] = edgeDifferences[i];
  }

  for(int col = 0; col < simplexDimension; col++) {
    // Partial pivot
    int pivot = col;
    for(int row = col + 1; row < simplexDimension; row++)
      if(std::abs(augmentedMatrix[row][col]) >
        std::abs(augmentedMatrix[pivot][col]))
        pivot = row;
    std::swap(augmentedMatrix[col], augmentedMatrix[pivot]);

    const float diagVal = augmentedMatrix[col][col];
    if(std::abs(diagVal) < powf(10, -FLT_DIG))
      return -2;

    for(int row = 0; row < simplexDimension; row++) {
      if(row == col)
        continue;
      const float factor = augmentedMatrix[row][col] / diagVal;
      for(int j = col; j <= simplexDimension; ++j)
        augmentedMatrix[row][j] -= factor * augmentedMatrix[col][j];
    }
  }

  std::vector<float> alpha(simplexDimension);
  for(int i = 0; i < simplexDimension; i++)
    alpha[i] = augmentedMatrix[i][simplexDimension] / augmentedMatrix[i][i];

  // reconstruct the 3D gradient
  for(int i = 0; i < simplexDimension; i++)
    for(int c = 0; c < 3; c++)
      gradient[c] += alpha[i] * edgeVectors[i][c];

  return 0;
}

template <class dataType, class triangulationType>
int ttk::nil::NumericalIntegralLines::execute(const triangulationType *triangulation,
    const std::vector<std::pair<SimplexId, int>> &seeds,
    std::vector<std::vector<PathPoint>> &output,
    const bool &isForward) const{

  Timer t;

  output.resize(seeds.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif
  for(int i = 0; i < (int) seeds.size(); i++){
    std::vector<float>
      barycentricWeights(seeds[i].second + 1, 1/(seeds[i].second + 1));
    computeIntegralLine<dataType, triangulationType>(
      triangulation, seeds[i], barycentricWeights, output[i], isForward);

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
    printMsg("  - Seed-#"
      + std::to_string(seeds[i].first)
      + " (dim: "
      + std::to_string(seeds[i].second)
      + ", f: "
      + std::to_string(isForward)
      + "): "
      + std::to_string(output[i].size()) + " point(s).",
        debug::Priority::DETAIL);
  }

  printMsg("Computed numerical integral line(s) from "
    + std::to_string(output.size())
    + " seed(s)"
    , 1,
    t.getElapsedTime(), threadNumber_);

  return 0;
}

template <class triangulationType>
  int ttk::nil::NumericalIntegralLines::getVertexIdentifiers(
    const triangulationType *triangulation,
    const int &simplexDimension, const int &simplexId,
    std::vector<SimplexId> vertexIdentifiers) const{

  switch(simplexDimension){
    case 0:
      vertexIdentifiers = {simplexId};
      break;
    case 1:
      vertexIdentifiers.resize(2);
      triangulation->getEdgeVertex(simplexId, 0, vertexIdentifiers[0]);
      triangulation->getEdgeVertex(simplexId, 1, vertexIdentifiers[1]);
      break;
    case 2:
      vertexIdentifiers.resize(3);
      triangulation->getTriangleVertex(simplexId, 0, vertexIdentifiers[0]);
      triangulation->getTriangleVertex(simplexId, 1, vertexIdentifiers[1]);
      triangulation->getTriangleVertex(simplexId, 2, vertexIdentifiers[2]);
      break;
    case 3:
      vertexIdentifiers.resize(4);
      triangulation->getCellVertex(simplexId, 0, vertexIdentifiers[0]);
      triangulation->getCellVertex(simplexId, 1, vertexIdentifiers[1]);
      triangulation->getCellVertex(simplexId, 2, vertexIdentifiers[2]);
      triangulation->getCellVertex(simplexId, 3, vertexIdentifiers[3]);
      break;
    default:
      return -1;
      break;
  }

  return 0;
}
