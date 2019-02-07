/// ingroup base
/// \class ttk::HarmonicFieldComputation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date February 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
/// Given an input scalar field and a list of critical points to remove, this
/// class minimally edits the scalar field such that the listed critical points
/// disappear. This procedure is useful to speedup subsequent topological data
/// analysis when outlier critical points can be easily identified. It is
/// also useful for data simplification.
///
/// \sa ttkHarmonicFieldComputation.cpp % for a usage example.

#pragma once

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>
#include <cmath>
#include <set>
#include <tuple>
#include <type_traits>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#endif // TTK_ENABLE_EIGEN

namespace ttk {

class HarmonicFieldComputation : public Debug {

public:
  HarmonicFieldComputation();

  ~HarmonicFieldComputation();

protected:
  Triangulation *triangulation_;
  SimplexId vertexNumber_;
  SimplexId constraintNumber_;
  void *inputScalarFieldPointer_;
  void *vertexIdentifierScalarFieldPointer_;
  void *inputOffsetScalarFieldPointer_;
  bool considerIdentifierAsBlackList_;
  bool addPerturbation_;
  void *outputScalarFieldPointer_;
  void *outputOffsetScalarFieldPointer_;
};
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <HarmonicFieldComputation.cpp>

template <typename dataType>
int ttk::HarmonicFieldComputation::getCriticalPoints(
    dataType *scalars, SimplexId *offsets, std::vector<SimplexId> &minima,
    std::vector<SimplexId> &maxima) const {
  std::vector<int> type(vertexNumber_, 0);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
  for (SimplexId k = 0; k < vertexNumber_; ++k)
    type[k] = getCriticalType<dataType>(k, scalars, offsets);

  for (SimplexId k = 0; k < vertexNumber_; ++k) {
    if (type[k] < 0)
      minima.push_back(k);
    else if (type[k] > 0)
      maxima.push_back(k);
  }
  return 0;
}

template <typename dataType>
int ttk::HarmonicFieldComputation::getCriticalPoints(
    dataType *scalars, SimplexId *offsets, std::vector<SimplexId> &minima,
    std::vector<SimplexId> &maxima, std::vector<bool> &extrema) const {
  std::vector<int> type(vertexNumber_);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for (SimplexId k = 0; k < vertexNumber_; ++k) {
    if (considerIdentifierAsBlackList_ xor extrema[k]) {
      type[k] = getCriticalType<dataType>(k, scalars, offsets);
    }
  }

  for (SimplexId k = 0; k < vertexNumber_; ++k) {
    if (type[k] < 0)
      minima.push_back(k);
    else if (type[k] > 0)
      maxima.push_back(k);
  }
  return 0;
}

template <typename dataType>
int ttk::HarmonicFieldComputation::addPerturbation(dataType *scalars,
                                                   SimplexId *offsets) const {
  dataType epsilon{};

  if (std::is_same<dataType, double>::value)
    epsilon = pow10(1 - DBL_DIG);
  else if (std::is_same<dataType, float>::value)
    epsilon = pow10(1 - FLT_DIG);
  else
    return -1;

  std::vector<std::tuple<dataType, SimplexId, SimplexId>> perturbation(
      vertexNumber_);
  for (SimplexId i = 0; i < vertexNumber_; ++i) {
    std::get<0>(perturbation[i]) = scalars[i];
    std::get<1>(perturbation[i]) = offsets[i];
    std::get<2>(perturbation[i]) = i;
  }

  SweepCmp cmp(true);
  sort(perturbation.begin(), perturbation.end(), cmp);

  for (SimplexId i = 0; i < vertexNumber_; ++i) {
    if (i) {
      if (std::get<0>(perturbation[i]) <= std::get<0>(perturbation[i - 1]))
        std::get<0>(perturbation[i]) =
            std::get<0>(perturbation[i - 1]) + epsilon;
    }
    scalars[std::get<2>(perturbation[i])] = std::get<0>(perturbation[i]);
  }

  return 0;
}

template <typename dataType, typename idType>
int ttk::HarmonicFieldComputation::execute() const {

  // get input data
  dataType *inputScalars = static_cast<dataType *>(inputScalarFieldPointer_);
  dataType *scalars = static_cast<dataType *>(outputScalarFieldPointer_);
  idType *identifiers =
      static_cast<idType *>(vertexIdentifierScalarFieldPointer_);
  idType *inputOffsets = static_cast<idType *>(inputOffsetScalarFieldPointer_);
  SimplexId *offsets =
      static_cast<SimplexId *>(outputOffsetScalarFieldPointer_);

  Timer t;

  {
    std::stringstream msg;
    msg << "[HarmonicFieldComputation] Beginnning computation" << std::endl;
    dMsg(std::cout, msg.str(), advancedInfoMsg);
  }

  // get the user extremum list
  std::vector<bool> extrema(vertexNumber_, false);
  for (SimplexId k = 0; k < constraintNumber_; ++k) {
    const SimplexId identifierId = identifiers[k];

#ifndef TTK_ENABLE_KAMIKAZE
    if (identifierId >= 0 and identifierId < vertexNumber_)
#endif
      extrema[identifierId] = true;
  }

  std::vector<SimplexId> authorizedMinima;
  std::vector<SimplexId> authorizedMaxima;
  std::vector<bool> authorizedExtrema(vertexNumber_, false);

  getCriticalPoints<dataType>(scalars, offsets, authorizedMinima,
                              authorizedMaxima, extrema);

  // processing
  int iteration{};
  for (SimplexId i = 0; i < vertexNumber_; ++i) {

    for (int j = 0; j < 2; ++j) {

      bool isIncreasingOrder = !j;

      cmp.setIsIncreasingOrder(isIncreasingOrder);
      std::set<std::tuple<dataType, SimplexId, SimplexId>, decltype(cmp)>
          sweepFront(cmp);
      std::vector<bool> visitedVertices(vertexNumber_, false);
      std::vector<SimplexId> adjustmentSequence(vertexNumber_);

      // add the seeds
      if (isIncreasingOrder) {
        for (SimplexId k : authorizedMinima) {
          authorizedExtrema[k] = true;
          sweepFront.emplace(scalars[k], offsets[k], k);
          visitedVertices[k] = true;
        }
      } else {
        for (SimplexId k : authorizedMaxima) {
          authorizedExtrema[k] = true;
          sweepFront.emplace(scalars[k], offsets[k], k);
          visitedVertices[k] = true;
        }
      }

      // growth by neighborhood of the seeds
      SimplexId adjustmentPos = 0;
      do {
        auto front = sweepFront.begin();
        if (front == sweepFront.end())
          return -1;

        SimplexId vertexId = std::get<2>(*front);
        sweepFront.erase(front);

        SimplexId neighborNumber =
            triangulation_->getVertexNeighborNumber(vertexId);
        for (SimplexId k = 0; k < neighborNumber; ++k) {
          SimplexId neighbor;
          triangulation_->getVertexNeighbor(vertexId, k, neighbor);
          if (!visitedVertices[neighbor]) {
            sweepFront.emplace(scalars[neighbor], offsets[neighbor], neighbor);
            visitedVertices[neighbor] = true;
          }
        }
        adjustmentSequence[adjustmentPos] = vertexId;
        ++adjustmentPos;
      } while (!sweepFront.empty());

      // save offsets and rearrange scalars
      SimplexId offset = (isIncreasingOrder ? 0 : vertexNumber_ + 1);

      for (SimplexId k = 0; k < vertexNumber_; ++k) {

        if (isIncreasingOrder) {
          if (k and scalars[adjustmentSequence[k]] <=
                        scalars[adjustmentSequence[k - 1]])
            scalars[adjustmentSequence[k]] = scalars[adjustmentSequence[k - 1]];
          ++offset;
        } else {
          if (k and scalars[adjustmentSequence[k]] >=
                        scalars[adjustmentSequence[k - 1]])
            scalars[adjustmentSequence[k]] = scalars[adjustmentSequence[k - 1]];
          --offset;
        }
        offsets[adjustmentSequence[k]] = offset;
      }
    }

    // test convergence
    bool needForMoreIterations{false};
    std::vector<SimplexId> minima;
    std::vector<SimplexId> maxima;
    getCriticalPoints<dataType>(scalars, offsets, minima, maxima);

    if (maxima.size() > authorizedMaxima.size())
      needForMoreIterations = true;
    if (minima.size() > authorizedMinima.size())
      needForMoreIterations = true;
  }

  {
    std::stringstream msg;
    msg << "[HarmonicFieldComputation] Ending computation after"
        << t.getElapsedTime() << std::endl;
    dMsg(std::cout, msg.str(), advancedInfoMsg);
  }

  return 0;
}
