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
#include <Geometry.h>
#include <Triangulation.h>
// std includes
#include <array>
#include <cfloat>
#include <cmath>
#include <limits>

namespace ttk {
  namespace nil {

    struct PathPoint {
      SimplexId simplexId_;
      int simplexDimension_;
      std::vector<float> barycentricWeights_;
    };

    /// Status of an elementary advection step (see doGradientStep()).
    enum StepStatus {
      /// The advection carries on (possibly in another simplex).
      REGULAR_STEP = 0,
      /// The advection reached the boundary of the domain.
      BOUNDARY_REACHED = 1,
      /// The advection reached a maximum (a minimum if backward).
      EXTREMUM_REACHED = 2
    };

    class NumericalIntegralLines : virtual public Debug {

    public:
      NumericalIntegralLines();
      ~NumericalIntegralLines() override;

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

      /**
       * @brief Compute the gradient of the piecewise linear scalar field,
       * restricted to the affine hull of the input simplex.
       *
       * @param simplexDimension Dimension of the input simplex.
       * @param simplexId Identifier of the input simplex.
       * @param gradient Output 3D gradient vector.
       * @param barycentricGradient Optional output expression of the gradient
       * in the edge basis (v1 - v0, ... vd - v0) of the simplex. This is also
       * the variation of the barycentric weights (but for the first one)
       * induced by a displacement along the gradient.
       */
      template <class dataType, class triangulationType>
      int computeNumericalGradient(const triangulationType *triangulation,
                                   const int &simplexDimension,
                                   const SimplexId &simplexId,
                                   std::vector<float> &gradient,
                                   std::vector<float> *barycentricGradient
                                   = nullptr) const;

      /**
       * @brief Elementary step of advection.
       *
       * Since the gradient of a piecewise linear scalar field is constant
       * within a simplex, the integral line is a straight segment there.
       * Hence, this step is integrated exactly: the current point is advected
       * within its simplex until it reaches its boundary. Then, the simplex in
       * which the advection carries on is identified (along with the
       * barycentric coordinates of the advected point within it).
       *
       * @param current Input point (simplex, dimension, barycentric weights).
       * @param isForward Forward or backward advection.
       * @param next Output point (simplex, dimension, barycentric weights).
       * @return StepStatus upon success (negative values otherwise). When the
       * advection cannot carry on (boundary of the domain or extremum), the
       * output simplex is the input one (with updated barycentric weights).
       */
      template <class dataType, class triangulationType>
      int doGradientStep(const triangulationType *triangulation,
                         const PathPoint &current,
                         const bool &isForward,
                         PathPoint &next) const;

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

      /**
       * @brief Compute the variation of the barycentric weights of a point
       * advected within the input simplex, per unit of arc length.
       *
       * @param slope Optional output slope of the scalar field along the
       * (unit) advection direction.
       * @return 0 if the advection can be carried on within the simplex
       * (negative values otherwise, in particular if the restriction of the
       * scalar field to the simplex is uniform).
       */
      template <class dataType, class triangulationType>
      int getBarycentricVelocity(const triangulationType *triangulation,
                                 const int &simplexDimension,
                                 const SimplexId &simplexId,
                                 const bool &isForward,
                                 std::vector<float> &velocity,
                                 float *slope = nullptr) const;

      /**
       * @brief Identify the simplex which takes the flow over at the input
       * point.
       *
       * The flow is taken over by the coface of the input simplex which
       * maximizes the slope of the scalar field (the steepest one) among the
       * cofaces which admit the advection. Note that the gradient restricted
       * to a face of a simplex is the projection of the gradient of the
       * simplex: the slope of a coface is therefore always larger than (or
       * equal to) that of its faces.
       *
       * @param point Input point (simplex, dimension, barycentric weights).
       * @param isForward Forward or backward advection.
       * @param next Output point (same location, expressed in the coface).
       * @return 0 if a coface takes the flow over (negative values otherwise,
       * in which case \p next is left untouched).
       */
      template <class dataType, class triangulationType>
      int getFlowSimplex(const triangulationType *triangulation,
                         const PathPoint &point,
                         const bool &isForward,
                         PathPoint &next) const;

      /**
       * @brief Retrieve the identifier of the face of the input simplex which
       * is spanned by the input vertices.
       */
      template <class triangulationType>
      int getFaceIdentifier(const triangulationType *triangulation,
                            const int &simplexDimension,
                            const SimplexId &simplexId,
                            const std::vector<SimplexId> &faceVertices,
                            const int &faceDimension,
                            SimplexId &faceIdentifier) const;

      /**
       * @brief Compute the 3D coordinates of a path point.
       */
      template <class triangulationType>
      int getPointCoordinates(const triangulationType *triangulation,
                              const PathPoint &point,
                              std::array<float, 3> &coordinates) const;

      /**
       * @brief Retrieve the cofaces of the input simplex (i.e. the simplices
       * of the star of the input simplex, of higher dimension), as a list of
       * (identifier, dimension) pairs.
       */
      template <class triangulationType>
      int getCofaces(const triangulationType *triangulation,
                     const int &simplexDimension,
                     const SimplexId &simplexId,
                     std::vector<std::pair<SimplexId, int>> &cofaces) const;

      /**
       * @brief Retrieve the face of the input simplex which supports the point
       * of input barycentric weights (i.e. the face spanned by the vertices of
       * non-zero weight).
       */
      template <class triangulationType>
      int getSubSimplex(const triangulationType *triangulation,
                        const int &simplexDimension,
                        const SimplexId &simplexId,
                        const std::vector<float> &barycentricWeights,
                        PathPoint &subSimplex) const;

      template <class triangulationType>
      int getVertexIdentifiers(const triangulationType *triangulation,
                               const int &simplexDimension,
                               const SimplexId &simplexId,
                               std::vector<SimplexId> &vertexIdentifiers) const;

      /**
       * @brief Check if the input simplex is on the boundary of the domain.
       */
      template <class triangulationType>
      bool isOnDomainBoundary(const triangulationType *triangulation,
                              const int &simplexDimension,
                              const SimplexId &simplexId) const;

      /**
       * @brief Check if an advection of input velocity can be carried on
       * within a simplex, from a point of input barycentric weights (i.e. the
       * advection does not immediately leave the simplex).
       */
      static inline bool
        isMotionAdmissible(const std::vector<float> &barycentricWeights,
                           const std::vector<float> &velocity) {

        float maximumVelocity = 0;
        for(int i = 0; i < (int)velocity.size(); i++)
          if(std::abs(velocity[i]) > maximumVelocity)
            maximumVelocity = std::abs(velocity[i]);

        if(!(maximumVelocity > 0))
          return false;

        for(int i = 0; i < (int)barycentricWeights.size(); i++)
          if((barycentricWeights[i] <= barycentricEpsilon_)
             && (velocity[i] < -relativeEpsilon_ * maximumVelocity))
            // the advection immediately exits through the i-th face
            return false;

        return true;
      }

      /**
       * @brief Express the barycentric weights of a point, given for a
       * simplex, in the basis of one of its cofaces.
       */
      static inline int
        mapBarycentricWeights(const std::vector<SimplexId> &sourceVertices,
                              const std::vector<float> &sourceWeights,
                              const std::vector<SimplexId> &targetVertices,
                              std::vector<float> &targetWeights) {

        targetWeights.clear();
        targetWeights.resize(targetVertices.size(), 0);

        for(int i = 0; i < (int)sourceVertices.size(); i++) {
          bool isFound = false;
          for(int j = 0; j < (int)targetVertices.size(); j++) {
            if(targetVertices[j] == sourceVertices[i]) {
              targetWeights[j] = sourceWeights[i];
              isFound = true;
              break;
            }
          }
          if(!isFound)
            // the source simplex is not a face of the target one
            return -1;
        }

        return 0;
      }

      /**
       * @brief Triangulation preconditioning.
       */
      inline void
        preconditionTriangulation(AbstractTriangulation *triangulation) {

        if(triangulation == nullptr)
          return;

        // precondition simplex2face
        triangulation->preconditionEdges();
        triangulation->preconditionCellEdges();

        // precondition face2cofacets
        triangulation->preconditionVertexEdges();
        triangulation->preconditionVertexStars();
        triangulation->preconditionEdgeStars();

        // precondition boundary
        triangulation->preconditionBoundaryVertices();
        if(triangulation->getDimensionality() > 1)
          // in 1D, edges are cells (and the boundary is made of vertices)
          triangulation->preconditionBoundaryEdges();

        if(triangulation->getDimensionality() == 3) {
          triangulation->preconditionTriangles();
          triangulation->preconditionTriangleEdges();
          triangulation->preconditionCellTriangles();
          triangulation->preconditionVertexTriangles();
          triangulation->preconditionEdgeTriangles();
          triangulation->preconditionTriangleStars();
          triangulation->preconditionBoundaryTriangles();
        }
      }

      inline void setInputScalarField(const void *const scalars) {
        scalars_ = scalars;
      }

    protected:
      /// Below this value, a barycentric weight is considered as null.
      static constexpr float barycentricEpsilon_{1e-6};
      /// Relative tolerance used for the null tests on the velocity.
      static constexpr float relativeEpsilon_{1e-6};

      int maximumIterationNumber_{1000000000};
      /// Number of consecutive steps without any motion after which the
      /// advection is considered as arbitrarily close to an extremum.
      int maximumStalledStepNumber_{8};
      const void *scalars_{};
    };
  } // namespace nil
} // namespace ttk

template <class dataType, class triangulationType>
int ttk::nil::NumericalIntegralLines::computeIntegralLine(
  const triangulationType *triangulation,
  const std::pair<SimplexId, int> &seed,
  const std::vector<float> &startBarycentricWeights,
  std::vector<PathPoint> &output,
  const bool &isForward) const {

  output.clear();

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation == nullptr)
    return -1;
  if(scalars_ == nullptr)
    return -2;
  if((seed.second < 0) || (seed.second > triangulation->getDimensionality()))
    return -3;
#endif

  PathPoint current;
  current.simplexId_ = seed.first;
  current.simplexDimension_ = seed.second;
  current.barycentricWeights_ = startBarycentricWeights;

  if((int)current.barycentricWeights_.size() != seed.second + 1)
    // no valid input coordinates: start from the barycenter of the seed
    current.barycentricWeights_.assign(
      seed.second + 1, 1.0 / (seed.second + 1));
  ttk::Geometry::normalizeBarycentricWeights(current.barycentricWeights_);

  // the seed simplex is where the integral line starts, not necessarily the
  // simplex within which the advection takes place. generically, a point in
  // the interior of a face is immediately advected within one of its cofaces
  // (for instance, the mid-point of an edge is taken over by one of the
  // triangles of its star). hand the flow over right away: this does not move
  // the point, it only re-expresses it in the coface. if no coface admits the
  // advection, the seed simplex constrains the flow and is kept as is.
  PathPoint flowPoint;
  if(getFlowSimplex<dataType, triangulationType>(
       triangulation, current, isForward, flowPoint)
     == 0)
    current = flowPoint;

  output.push_back(current);

  std::array<float, 3> previousCoordinates{}, currentCoordinates{};
  getPointCoordinates(triangulation, current, previousCoordinates);

  // an advection step may legitimately not move the current point (it can
  // simply update the simplex supporting it, for instance when leaving a
  // vertex for one of the tetrahedra of its star). however, a point which no
  // longer moves is arbitrarily close to an extremum.
  int stalledStepNumber = 0;
  float pathLength = 0;

  int step = 0, status = REGULAR_STEP;

  for(step = 0; step < maximumIterationNumber_; step++) {

    PathPoint next;

    status = doGradientStep<dataType, triangulationType>(
      triangulation, current, isForward, next);

    if(status < 0)
      return status;

    getPointCoordinates(triangulation, next, currentCoordinates);

    const float stepLength = Geometry::distance(
      previousCoordinates.data(), currentCoordinates.data());

    if(stepLength > pathLength * std::numeric_limits<float>::epsilon()) {
      output.push_back(next);
      pathLength += stepLength;
      stalledStepNumber = 0;
    } else {
      // the point did not move: only update the simplex supporting it
      output.back() = next;
      stalledStepNumber++;
    }

    current = next;
    previousCoordinates = currentCoordinates;

    if(status != REGULAR_STEP)
      // the advection either left the domain or reached an extremum
      break;

    if(stalledStepNumber > maximumStalledStepNumber_) {
      // the advection is arbitrarily close to an extremum
      status = EXTREMUM_REACHED;
      break;
    }
  }

  if(step == maximumIterationNumber_) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
    printWrn("Maximum iteration number reached for seed-#"
             + std::to_string(seed.first)
             + " (dim: " + std::to_string(seed.second) + ").");
  }

  return 0;
}

// TODO
// move that function to the geometry class

template <class dataType, class triangulationType>
int ttk::nil::NumericalIntegralLines::computeNumericalGradient(
  const triangulationType *triangulation,
  const int &simplexDimension,
  const SimplexId &simplexId,
  std::vector<float> &gradient,
  std::vector<float> *barycentricGradient) const {

  gradient = {0, 0, 0};

  if(barycentricGradient)
    barycentricGradient->assign(simplexDimension, 0);

  if(!simplexDimension)
    return -1;

  std::vector<SimplexId> vertexIdentifiers;

  getVertexIdentifiers(
    triangulation, simplexDimension, simplexId, vertexIdentifiers);

  const int vertexNumber = vertexIdentifiers.size();

  std::vector<std::array<float, 3>> vertexPoints(vertexNumber);
  std::vector<dataType> vertexScalars(vertexNumber);

  for(int i = 0; i < (int)vertexNumber; i++) {
    triangulation->getVertexPoint(vertexIdentifiers[i], vertexPoints[i][0],
                                  vertexPoints[i][1], vertexPoints[i][2]);
    vertexScalars[i] = ((const dataType *)scalars_)[vertexIdentifiers[i]];
  }

  // build edge vectors and corresponding differences, wrt v0
  std::vector<std::array<float, 3>> edgeVectors(simplexDimension);
  std::vector<float> edgeDifferences(simplexDimension);

  for(int i = 0; i < simplexDimension; i++) {
    for(int c = 0; c < 3; c++)
      edgeVectors[i][c] = vertexPoints[i + 1][c] - vertexPoints[0][c];
    edgeDifferences[i]
      = ((float)vertexScalars[i + 1]) - ((float)vertexScalars[0]);
  }

  // Gram matrix  gramMatrix[i][j] = edgeVectors[i] . edgeVectors[j]
  std::vector<std::vector<float>> gramMatrix(
    simplexDimension, std::vector<float>(simplexDimension, 0));

  float maximumDiagonalEntry = 0;

  for(int i = 0; i < simplexDimension; i++) {
    for(int j = 0; j < simplexDimension; j++)
      gramMatrix[i][j] = ttk::Geometry::dotProduct(
        edgeVectors[i].data(), edgeVectors[j].data());

    if(gramMatrix[i][i] > maximumDiagonalEntry)
      maximumDiagonalEntry = gramMatrix[i][i];
  }

  if(!(maximumDiagonalEntry > 0))
    // degenerated simplex
    return -2;

  // Gaussian elimintation
  std::vector<std::vector<float>> augmentedMatrix(
    simplexDimension, std::vector<float>(simplexDimension + 1));
  for(int i = 0; i < simplexDimension; ++i) {
    for(int j = 0; j < simplexDimension; ++j)
      augmentedMatrix[i][j] = gramMatrix[i][j];
    augmentedMatrix[i][simplexDimension] = edgeDifferences[i];
  }

  for(int col = 0; col < simplexDimension; col++) {
    // Partial pivot
    int pivot = col;
    for(int row = col + 1; row < simplexDimension; row++)
      if(std::abs(augmentedMatrix[row][col])
         > std::abs(augmentedMatrix[pivot][col]))
        pivot = row;
    std::swap(augmentedMatrix[col], augmentedMatrix[pivot]);

    const float diagVal = augmentedMatrix[col][col];
    if(std::abs(diagVal) < powf(10, -FLT_DIG) * maximumDiagonalEntry)
      // degenerated simplex
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

  if(barycentricGradient)
    *barycentricGradient = alpha;

  return 0;
}

template <class dataType, class triangulationType>
int ttk::nil::NumericalIntegralLines::doGradientStep(
  const triangulationType *triangulation,
  const PathPoint &current,
  const bool &isForward,
  PathPoint &next) const {

  next = current;

  // 1) advection within the current simplex.
  // the gradient of a piecewise linear scalar field is constant within a
  // simplex. hence, the integral line is a straight segment there, which can
  // be integrated exactly: the point is advected until it reaches the boundary
  // of the simplex.
  std::vector<float> weights = current.barycentricWeights_;
  std::vector<float> velocity;
  bool hasMoved = false;

  if(getBarycentricVelocity<dataType, triangulationType>(
       triangulation, current.simplexDimension_, current.simplexId_, isForward,
       velocity)
     == 0) {

    float maximumVelocity = 0;
    for(int i = 0; i <= current.simplexDimension_; i++)
      if(std::abs(velocity[i]) > maximumVelocity)
        maximumVelocity = std::abs(velocity[i]);

    // largest arc length which maintains the point within the simplex
    float travelDistance = std::numeric_limits<float>::infinity();
    for(int i = 0; i <= current.simplexDimension_; i++) {
      if(velocity[i] < -relativeEpsilon_ * maximumVelocity) {
        const float distance = current.barycentricWeights_[i] / (-velocity[i]);
        if(distance < travelDistance)
          travelDistance = distance;
      }
    }

    if((travelDistance > 0)
       && (travelDistance < std::numeric_limits<float>::infinity())) {

      for(int i = 0; i <= current.simplexDimension_; i++)
        weights[i]
          = current.barycentricWeights_[i] + travelDistance * velocity[i];
      ttk::Geometry::normalizeBarycentricWeights(weights);

      hasMoved = true;
    }
  }

  // the advected point is now supported by a face of the current simplex
  PathPoint exitPoint;
  if(getSubSimplex(triangulation, current.simplexDimension_, current.simplexId_,
                   weights, exitPoint)
     < 0)
    return -1;

  // 2) identify the simplex in which the advection carries on.
  if(getFlowSimplex<dataType, triangulationType>(
       triangulation, exitPoint, isForward, next)
     == 0)
    return REGULAR_STEP;

  // no coface takes the flow over.
  const bool isOnBoundary = isOnDomainBoundary(
    triangulation, exitPoint.simplexDimension_, exitPoint.simplexId_);

  if((!isOnBoundary) && (exitPoint.simplexDimension_ > 0)
     && ((exitPoint.simplexDimension_ != current.simplexDimension_)
         || (exitPoint.simplexId_ != current.simplexId_))) {

    // in the interior of the domain, the flow is then constrained to the exit
    // face itself (typically, two cells whose gradients both point towards
    // their common face).
    if(getBarycentricVelocity<dataType, triangulationType>(
         triangulation, exitPoint.simplexDimension_, exitPoint.simplexId_,
         isForward, velocity)
       == 0) {

      if(isMotionAdmissible(exitPoint.barycentricWeights_, velocity)) {
        next = exitPoint;
        return REGULAR_STEP;
      }
    }
  }

  // the advection stops here: report the advected point within the current
  // simplex (same simplex, different barycentric weights).
  next.simplexId_ = current.simplexId_;
  next.simplexDimension_ = current.simplexDimension_;
  next.barycentricWeights_ = weights;

  if(hasMoved && isOnBoundary)
    return BOUNDARY_REACHED;

  return EXTREMUM_REACHED;
}

template <class dataType, class triangulationType>
int ttk::nil::NumericalIntegralLines::execute(
  const triangulationType *triangulation,
  const std::vector<std::pair<SimplexId, int>> &seeds,
  std::vector<std::vector<PathPoint>> &output,
  const bool &isForward) const {

  Timer t;

  output.resize(seeds.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif
  for(int i = 0; i < (int)seeds.size(); i++) {
    std::vector<float> barycentricWeights(
      seeds[i].second + 1, 1.0 / (seeds[i].second + 1));
    computeIntegralLine<dataType, triangulationType>(
      triangulation, seeds[i], barycentricWeights, output[i], isForward);

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
    printMsg("  - Seed-#" + std::to_string(seeds[i].first)
               + " (dim: " + std::to_string(seeds[i].second)
               + ", f: " + std::to_string(isForward)
               + "): " + std::to_string(output[i].size()) + " point(s).",
             debug::Priority::DETAIL);
  }

  printMsg("Computed from " + std::to_string(output.size()) + " seed(s)", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}

template <class dataType, class triangulationType>
int ttk::nil::NumericalIntegralLines::getBarycentricVelocity(
  const triangulationType *triangulation,
  const int &simplexDimension,
  const SimplexId &simplexId,
  const bool &isForward,
  std::vector<float> &velocity,
  float *slope) const {

  velocity.clear();

  if(slope)
    *slope = 0;

  if(simplexDimension < 1)
    // no motion within a vertex
    return -1;

  std::vector<float> gradient, barycentricGradient;

  if(computeNumericalGradient<dataType, triangulationType>(
       triangulation, simplexDimension, simplexId, gradient,
       &barycentricGradient)
     < 0)
    return -2;

  // along the unit advection direction, the variation of the scalar field is
  // given by the magnitude of the gradient
  const float magnitude = ttk::Geometry::magnitude(gradient.data());

  if(!(magnitude > 0))
    // uniform scalar field: no motion
    return -3;

  if(slope)
    *slope = magnitude;

  // unit speed advection (the integration variable is the arc length)
  const float scale = (isForward ? 1.0 : -1.0) / magnitude;

  velocity.resize(simplexDimension + 1, 0);
  for(int i = 0; i < simplexDimension; i++) {
    velocity[i + 1] = scale * barycentricGradient[i];
    // the barycentric weights sum up to 1
    velocity[0] -= scale * barycentricGradient[i];
  }

  return 0;
}

template <class dataType, class triangulationType>
int ttk::nil::NumericalIntegralLines::getFlowSimplex(
  const triangulationType *triangulation,
  const PathPoint &point,
  const bool &isForward,
  PathPoint &next) const {

  // the flow is taken over by the coface of the input simplex which maximizes
  // the slope of the scalar field (the steepest one) among the cofaces which
  // admit the advection. note that the gradient restricted to a face of a
  // simplex is the projection of the gradient of the simplex: the slope of a
  // coface is therefore always larger than (or equal to) that of its faces.
  std::vector<std::pair<SimplexId, int>> cofaces;
  getCofaces(triangulation, point.simplexDimension_, point.simplexId_, cofaces);

  std::vector<SimplexId> pointVertices, cofaceVertices;
  getVertexIdentifiers(
    triangulation, point.simplexDimension_, point.simplexId_, pointVertices);

  std::vector<float> velocity;

  PathPoint bestPoint;
  bestPoint.simplexDimension_ = -1;
  float bestSlope = 0;

  for(int i = 0; i < (int)cofaces.size(); i++) {

    PathPoint candidate;
    candidate.simplexId_ = cofaces[i].first;
    candidate.simplexDimension_ = cofaces[i].second;

    getVertexIdentifiers(triangulation, candidate.simplexDimension_,
                         candidate.simplexId_, cofaceVertices);

    // express the advected point in the barycentric basis of the coface
    if(mapBarycentricWeights(pointVertices, point.barycentricWeights_,
                             cofaceVertices, candidate.barycentricWeights_)
       < 0)
      continue;

    float slope = 0;
    if(getBarycentricVelocity<dataType, triangulationType>(
         triangulation, candidate.simplexDimension_, candidate.simplexId_,
         isForward, velocity, &slope)
       < 0)
      continue;

    if(!isMotionAdmissible(candidate.barycentricWeights_, velocity))
      // the advection would immediately leave this coface
      continue;

    // steepest slope: along a unit direction, the variation of the scalar
    // field is given by the magnitude of the gradient.
    // ties (the gradient of the coface is aligned with one of its faces) are
    // settled in favor of the coface of highest dimension (i.e. the least
    // constrained advection).
    if((slope > bestSlope)
       || ((slope > bestSlope * (1 - relativeEpsilon_))
           && (candidate.simplexDimension_ > bestPoint.simplexDimension_))) {
      bestSlope = slope;
      bestPoint = candidate;
    }
  }

  if(bestPoint.simplexDimension_ < 0)
    // no coface takes the flow over
    return -1;

  next = bestPoint;

  return 0;
}

template <class triangulationType>
int ttk::nil::NumericalIntegralLines::getFaceIdentifier(
  const triangulationType *triangulation,
  const int &simplexDimension,
  const SimplexId &simplexId,
  const std::vector<SimplexId> &faceVertices,
  const int &faceDimension,
  SimplexId &faceIdentifier) const {

  faceIdentifier = -1;

  if(faceDimension == simplexDimension) {
    faceIdentifier = simplexId;
    return 0;
  }

  if(!faceDimension) {
    faceIdentifier = faceVertices[0];
    return 0;
  }

  const int cellDimension = triangulation->getDimensionality();

  int faceNumber = 0;
  if(faceDimension == 1) {
    if(simplexDimension == cellDimension)
      faceNumber = triangulation->getCellEdgeNumber(simplexId);
    else
      // edges of a triangle
      faceNumber = 3;
  } else if(faceDimension == 2)
    // triangles of a tetrahedron
    faceNumber = triangulation->getCellTriangleNumber(simplexId);
  else
    return -1;

  // identify the face of the simplex which spans the input vertices
  std::vector<SimplexId> candidateVertices;

  for(int i = 0; i < faceNumber; i++) {

    SimplexId candidateId = -1;

    if(faceDimension == 1) {
      if(simplexDimension == cellDimension)
        triangulation->getCellEdge(simplexId, i, candidateId);
      else
        triangulation->getTriangleEdge(simplexId, i, candidateId);
    } else
      triangulation->getCellTriangle(simplexId, i, candidateId);

    getVertexIdentifiers(
      triangulation, faceDimension, candidateId, candidateVertices);

    bool isMatching = true;
    for(int j = 0; j < (int)faceVertices.size(); j++) {
      bool isFound = false;
      for(int k = 0; k < (int)candidateVertices.size(); k++) {
        if(candidateVertices[k] == faceVertices[j]) {
          isFound = true;
          break;
        }
      }
      if(!isFound) {
        isMatching = false;
        break;
      }
    }

    if(isMatching) {
      faceIdentifier = candidateId;
      return 0;
    }
  }

  return -2;
}

template <class triangulationType>
int ttk::nil::NumericalIntegralLines::getPointCoordinates(
  const triangulationType *triangulation,
  const PathPoint &point,
  std::array<float, 3> &coordinates) const {

  coordinates = {0, 0, 0};

  std::vector<SimplexId> vertexIdentifiers;

  if(getVertexIdentifiers(triangulation, point.simplexDimension_,
                          point.simplexId_, vertexIdentifiers)
     < 0)
    return -1;

  for(int i = 0; i < (int)vertexIdentifiers.size(); i++) {
    std::array<float, 3> vertexPoint;
    triangulation->getVertexPoint(
      vertexIdentifiers[i], vertexPoint[0], vertexPoint[1], vertexPoint[2]);
    for(int c = 0; c < 3; c++)
      coordinates[c] += point.barycentricWeights_[i] * vertexPoint[c];
  }

  return 0;
}

template <class triangulationType>
int ttk::nil::NumericalIntegralLines::getCofaces(
  const triangulationType *triangulation,
  const int &simplexDimension,
  const SimplexId &simplexId,
  std::vector<std::pair<SimplexId, int>> &cofaces) const {

  cofaces.clear();

  const int cellDimension = triangulation->getDimensionality();

  if((simplexDimension < 0) || (simplexDimension >= cellDimension))
    // a top dimensional cell has no coface
    return -1;

  SimplexId cofaceId = -1;

  if(!simplexDimension) {
    if(cellDimension > 1) {
      // the edges of the star of the vertex
      // (in 1D, edges are cells: they are collected below)
      const SimplexId edgeNumber
        = triangulation->getVertexEdgeNumber(simplexId);
      for(SimplexId i = 0; i < edgeNumber; i++) {
        triangulation->getVertexEdge(simplexId, i, cofaceId);
        cofaces.emplace_back(cofaceId, 1);
      }
    }

    if(cellDimension == 3) {
      // the triangles of the star of the vertex
      const SimplexId triangleNumber
        = triangulation->getVertexTriangleNumber(simplexId);
      for(SimplexId i = 0; i < triangleNumber; i++) {
        triangulation->getVertexTriangle(simplexId, i, cofaceId);
        cofaces.emplace_back(cofaceId, 2);
      }
    }

    // the cells of the star of the vertex
    const SimplexId starNumber = triangulation->getVertexStarNumber(simplexId);
    for(SimplexId i = 0; i < starNumber; i++) {
      triangulation->getVertexStar(simplexId, i, cofaceId);
      cofaces.emplace_back(cofaceId, cellDimension);
    }

    return 0;
  }

  if(simplexDimension == 1) {
    if(cellDimension == 3) {
      // the triangles of the star of the edge
      const SimplexId triangleNumber
        = triangulation->getEdgeTriangleNumber(simplexId);
      for(SimplexId i = 0; i < triangleNumber; i++) {
        triangulation->getEdgeTriangle(simplexId, i, cofaceId);
        cofaces.emplace_back(cofaceId, 2);
      }
    }

    // the cells of the star of the edge
    const SimplexId starNumber = triangulation->getEdgeStarNumber(simplexId);
    for(SimplexId i = 0; i < starNumber; i++) {
      triangulation->getEdgeStar(simplexId, i, cofaceId);
      cofaces.emplace_back(cofaceId, cellDimension);
    }

    return 0;
  }

  // the cells of the star of the triangle
  const SimplexId starNumber = triangulation->getTriangleStarNumber(simplexId);
  for(SimplexId i = 0; i < starNumber; i++) {
    triangulation->getTriangleStar(simplexId, i, cofaceId);
    cofaces.emplace_back(cofaceId, cellDimension);
  }

  return 0;
}

template <class triangulationType>
int ttk::nil::NumericalIntegralLines::getSubSimplex(
  const triangulationType *triangulation,
  const int &simplexDimension,
  const SimplexId &simplexId,
  const std::vector<float> &barycentricWeights,
  PathPoint &subSimplex) const {

  std::vector<SimplexId> vertexIdentifiers;

  if(getVertexIdentifiers(
       triangulation, simplexDimension, simplexId, vertexIdentifiers)
     < 0)
    return -1;

  // the point is supported by the face spanned by the vertices of non-zero
  // barycentric weight
  std::vector<SimplexId> faceVertices;
  std::vector<float> faceWeights;

  for(int i = 0; i < (int)vertexIdentifiers.size(); i++) {
    if(barycentricWeights[i] > barycentricEpsilon_) {
      faceVertices.push_back(vertexIdentifiers[i]);
      faceWeights.push_back(barycentricWeights[i]);
    }
  }

  if(faceVertices.empty()) {
    // degenerated weights: fall back on the closest vertex
    int closestVertex = 0;
    for(int i = 1; i < (int)vertexIdentifiers.size(); i++)
      if(barycentricWeights[i] > barycentricWeights[closestVertex])
        closestVertex = i;
    faceVertices = {vertexIdentifiers[closestVertex]};
    faceWeights = {1};
  }

  const int faceDimension = faceVertices.size() - 1;

  if(faceDimension == simplexDimension) {
    subSimplex.simplexId_ = simplexId;
    subSimplex.simplexDimension_ = simplexDimension;
    subSimplex.barycentricWeights_ = barycentricWeights;
    ttk::Geometry::normalizeBarycentricWeights(subSimplex.barycentricWeights_);
    return 0;
  }

  SimplexId faceIdentifier = -1;

  if(getFaceIdentifier(triangulation, simplexDimension, simplexId, faceVertices,
                       faceDimension, faceIdentifier)
     < 0)
    return -2;

  subSimplex.simplexId_ = faceIdentifier;
  subSimplex.simplexDimension_ = faceDimension;

  // the vertices of the face are not necessarily ordered as in the simplex
  std::vector<SimplexId> subVertexIdentifiers;
  getVertexIdentifiers(
    triangulation, faceDimension, faceIdentifier, subVertexIdentifiers);

  if(mapBarycentricWeights(faceVertices, faceWeights, subVertexIdentifiers,
                           subSimplex.barycentricWeights_)
     < 0)
    return -3;

  ttk::Geometry::normalizeBarycentricWeights(subSimplex.barycentricWeights_);

  return 0;
}

template <class triangulationType>
int ttk::nil::NumericalIntegralLines::getVertexIdentifiers(
  const triangulationType *triangulation,
  const int &simplexDimension,
  const SimplexId &simplexId,
  std::vector<SimplexId> &vertexIdentifiers) const {

  if((simplexDimension < 0) || (simplexDimension > 3))
    return -1;

  vertexIdentifiers.resize(simplexDimension + 1);

  if(simplexDimension == triangulation->getDimensionality()) {
    // top dimensional simplex: use the (faster) cell accessors
    // (in 2D, triangles are cells)
    for(int i = 0; i < simplexDimension + 1; i++)
      triangulation->getCellVertex(simplexId, i, vertexIdentifiers[i]);
    return 0;
  }

  switch(simplexDimension) {
    case 0:
      vertexIdentifiers[0] = simplexId;
      break;
    case 1:
      triangulation->getEdgeVertex(simplexId, 0, vertexIdentifiers[0]);
      triangulation->getEdgeVertex(simplexId, 1, vertexIdentifiers[1]);
      break;
    case 2:
      triangulation->getTriangleVertex(simplexId, 0, vertexIdentifiers[0]);
      triangulation->getTriangleVertex(simplexId, 1, vertexIdentifiers[1]);
      triangulation->getTriangleVertex(simplexId, 2, vertexIdentifiers[2]);
      break;
    default:
      return -1;
      break;
  }

  return 0;
}

template <class triangulationType>
bool ttk::nil::NumericalIntegralLines::isOnDomainBoundary(
  const triangulationType *triangulation,
  const int &simplexDimension,
  const SimplexId &simplexId) const {

  if(simplexDimension == triangulation->getDimensionality())
    // a top dimensional cell is never on the boundary of the domain
    return false;

  switch(simplexDimension) {
    case 0:
      return triangulation->isVertexOnBoundary(simplexId);
    case 1:
      return triangulation->isEdgeOnBoundary(simplexId);
    case 2:
      return triangulation->isTriangleOnBoundary(simplexId);
  }

  return false;
}
