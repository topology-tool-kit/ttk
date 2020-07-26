/// \ingroup base
/// \class ttk::TopologicalSimplification
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date February 2016
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
/// \b Related \b publication \n
/// "Generalized Topological Simplification of Scalar Fields on Surfaces" \n
/// Julien Tierny, Valerio Pascucci \n
/// Proc. of IEEE VIS 2012.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2012.
///
/// \sa ttkTopologicalSimplification.cpp %for a usage example.

#pragma once

// base code includes
#include <Debug.h>
#include <Triangulation.h>

#include <cmath>
#include <set>
#include <tuple>
#include <type_traits>

namespace ttk {

  class SweepCmp {
    bool isIncreasingOrder_{};

  public:
    SweepCmp() {
    }

    SweepCmp(bool isIncreasingOrder) : isIncreasingOrder_{isIncreasingOrder} {
    }

    inline void setIsIncreasingOrder(bool isIncreasingOrder) {
      isIncreasingOrder_ = isIncreasingOrder;
    }

    template <typename dataType>
    bool
      operator()(const std::tuple<dataType, SimplexId, SimplexId> &v0,
                 const std::tuple<dataType, SimplexId, SimplexId> &v1) const {
      if(isIncreasingOrder_) {
        return (std::get<0>(v0) < std::get<0>(v1)
                or (std::get<0>(v0) == std::get<0>(v1)
                    and std::get<1>(v0) < std::get<1>(v1)));
      } else {
        return (std::get<0>(v0) > std::get<0>(v1)
                or (std::get<0>(v0) == std::get<0>(v1)
                    and std::get<1>(v0) > std::get<1>(v1)));
      }
    }
  };

  class TopologicalSimplification : virtual public Debug {
  public:
    TopologicalSimplification();

    template <typename dataType>
    bool isLowerThan(SimplexId a,
                     SimplexId b,
                     dataType *scalars,
                     SimplexId *offsets) const;

    template <typename dataType>
    bool isHigherThan(SimplexId a,
                      SimplexId b,
                      dataType *scalars,
                      SimplexId *offsets) const;

    template <typename dataType, typename triangulationType>
    int getCriticalType(SimplexId vertexId,
                        dataType *scalars,
                        SimplexId *offsets,
                        const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int getCriticalPoints(dataType *scalars,
                          SimplexId *offsets,
                          std::vector<SimplexId> &minList,
                          std::vector<SimplexId> &maxList,
                          const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int getCriticalPoints(dataType *scalars,
                          SimplexId *offsets,
                          std::vector<SimplexId> &minList,
                          std::vector<SimplexId> &maxList,
                          std::vector<bool> &blackList,
                          const triangulationType &triangulation) const;

    template <typename dataType>
    int addPerturbation(dataType *scalars, SimplexId *offsets) const;

    template <typename dataType, typename idType, typename triangulationType>
    int execute(const dataType *const inputScalars,
                dataType *const outputScalars,
                const idType *const identifiers,
                const idType *const inputOffsets,
                SimplexId *const offsets,
                const triangulationType &triangulation) const;

    inline int preconditionTriangulation(AbstractTriangulation *triangulation) {
      if(triangulation) {
        vertexNumber_ = triangulation->getNumberOfVertices();
        triangulation->preconditionVertexNeighbors();
      }
      return 0;
    }

    inline void setVertexNumber(SimplexId vertexNumber) {
      vertexNumber_ = vertexNumber;
    }

    inline void setConstraintNumber(SimplexId constraintNumber) {
      constraintNumber_ = constraintNumber;
    }

    inline void setInputScalarFieldPointer(void *data) {
      inputScalarFieldPointer_ = data;
    }

    inline void setVertexIdentifierScalarFieldPointer(void *data) {
      vertexIdentifierScalarFieldPointer_ = data;
    }

    inline void setInputOffsetScalarFieldPointer(void *data) {
      inputOffsetScalarFieldPointer_ = data;
    }

    inline void setConsiderIdentifierAsBlackList(bool onOff) {
      considerIdentifierAsBlackList_ = onOff;
    }

    inline void setAddPerturbation(bool onOff) {
      addPerturbation_ = onOff;
    }

    inline void setOutputScalarFieldPointer(void *data) {
      outputScalarFieldPointer_ = data;
    }

    inline void setOutputOffsetScalarFieldPointer(void *data) {
      outputOffsetScalarFieldPointer_ = data;
    }

  protected:
    SimplexId vertexNumber_{};
    SimplexId constraintNumber_{};
    void *inputScalarFieldPointer_{};
    void *vertexIdentifierScalarFieldPointer_{};
    void *inputOffsetScalarFieldPointer_{};
    bool considerIdentifierAsBlackList_{false};
    bool addPerturbation_{false};
    void *outputScalarFieldPointer_{};
    void *outputOffsetScalarFieldPointer_{};
  };
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <TopologicalSimplification.cpp>

template <typename dataType>
bool ttk::TopologicalSimplification::isLowerThan(SimplexId a,
                                                 SimplexId b,
                                                 dataType *scalars,
                                                 SimplexId *offsets) const {
  return (scalars[a] < scalars[b]
          or (scalars[a] == scalars[b] and offsets[a] < offsets[b]));
}

template <typename dataType>
bool ttk::TopologicalSimplification::isHigherThan(SimplexId a,
                                                  SimplexId b,
                                                  dataType *scalars,
                                                  SimplexId *offsets) const {
  return (scalars[a] > scalars[b]
          or (scalars[a] == scalars[b] and offsets[a] > offsets[b]));
}

template <typename dataType, typename triangulationType>
int ttk::TopologicalSimplification::getCriticalType(
  SimplexId vertex,
  dataType *scalars,
  SimplexId *offsets,
  const triangulationType &triangulation) const {

  bool isMinima{true};
  bool isMaxima{true};
  SimplexId neighborNumber = triangulation.getVertexNeighborNumber(vertex);
  for(SimplexId i = 0; i < neighborNumber; ++i) {
    SimplexId neighbor;
    triangulation.getVertexNeighbor(vertex, i, neighbor);

    if(isLowerThan<dataType>(neighbor, vertex, scalars, offsets))
      isMinima = false;
    if(isHigherThan<dataType>(neighbor, vertex, scalars, offsets))
      isMaxima = false;
    if(!isMinima and !isMaxima) {
      return 0;
    }
  }

  if(isMinima)
    return -1;
  if(isMaxima)
    return 1;

  return 0;
}

template <typename dataType, typename triangulationType>
int ttk::TopologicalSimplification::getCriticalPoints(
  dataType *scalars,
  SimplexId *offsets,
  std::vector<SimplexId> &minima,
  std::vector<SimplexId> &maxima,
  const triangulationType &triangulation) const {

  std::vector<int> type(vertexNumber_, 0);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
  for(SimplexId k = 0; k < vertexNumber_; ++k)
    type[k] = getCriticalType<dataType>(k, scalars, offsets, triangulation);

  for(SimplexId k = 0; k < vertexNumber_; ++k) {
    if(type[k] < 0)
      minima.push_back(k);
    else if(type[k] > 0)
      maxima.push_back(k);
  }
  return 0;
}

template <typename dataType, typename triangulationType>
int ttk::TopologicalSimplification::getCriticalPoints(
  dataType *scalars,
  SimplexId *offsets,
  std::vector<SimplexId> &minima,
  std::vector<SimplexId> &maxima,
  std::vector<bool> &extrema,
  const triangulationType &triangulation) const {
  std::vector<int> type(vertexNumber_);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId k = 0; k < vertexNumber_; ++k) {
    if(considerIdentifierAsBlackList_ xor extrema[k]) {
      type[k] = getCriticalType<dataType>(k, scalars, offsets, triangulation);
    }
  }

  for(SimplexId k = 0; k < vertexNumber_; ++k) {
    if(type[k] < 0)
      minima.push_back(k);
    else if(type[k] > 0)
      maxima.push_back(k);
  }
  return 0;
}

template <typename dataType>
int ttk::TopologicalSimplification::addPerturbation(dataType *scalars,
                                                    SimplexId *offsets) const {
  dataType epsilon{};

  if(std::is_same<dataType, double>::value)
    epsilon = Geometry::powIntTen<dataType>(1 - DBL_DIG);
  else if(std::is_same<dataType, float>::value)
    epsilon = Geometry::powIntTen<dataType>(1 - FLT_DIG);
  else
    return -1;

  std::vector<std::tuple<dataType, SimplexId, SimplexId>> perturbation(
    vertexNumber_);
  for(SimplexId i = 0; i < vertexNumber_; ++i) {
    std::get<0>(perturbation[i]) = scalars[i];
    std::get<1>(perturbation[i]) = offsets[i];
    std::get<2>(perturbation[i]) = i;
  }

  SweepCmp cmp(true);
  sort(perturbation.begin(), perturbation.end(), cmp);

  for(SimplexId i = 0; i < vertexNumber_; ++i) {
    if(i) {
      if(std::get<0>(perturbation[i]) <= std::get<0>(perturbation[i - 1]))
        std::get<0>(perturbation[i])
          = std::get<0>(perturbation[i - 1]) + epsilon;
    }
    scalars[std::get<2>(perturbation[i])] = std::get<0>(perturbation[i]);
  }

  return 0;
}

template <typename dataType, typename idType, typename triangulationType>
int ttk::TopologicalSimplification::execute(
  const dataType *const inputScalars,
  dataType *const outputScalars,
  const idType *const identifiers,
  const idType *const inputOffsets,
  SimplexId *const offsets,
  const triangulationType &triangulation) const {

  Timer t;

  // pre-processing
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId k = 0; k < vertexNumber_; ++k) {
    outputScalars[k] = inputScalars[k];
    if(std::isnan((double)outputScalars[k]))
      outputScalars[k] = 0;

    offsets[k] = inputOffsets[k];
  }

  // get the user extremum list
  std::vector<bool> extrema(vertexNumber_, false);
  for(SimplexId k = 0; k < constraintNumber_; ++k) {
    const SimplexId identifierId = identifiers[k];

#ifndef TTK_ENABLE_KAMIKAZE
    if(identifierId >= 0 and identifierId < vertexNumber_)
#endif
      extrema[identifierId] = true;
  }

  std::vector<SimplexId> authorizedMinima;
  std::vector<SimplexId> authorizedMaxima;
  std::vector<bool> authorizedExtrema(vertexNumber_, false);

  getCriticalPoints(outputScalars, offsets, authorizedMinima, authorizedMaxima,
                    extrema, triangulation);

  this->printMsg("Maintaining " + std::to_string(constraintNumber_)
                   + " constraints (" + std::to_string(authorizedMinima.size())
                   + " minima and " + std::to_string(authorizedMaxima.size())
                   + "maxima)",
                 debug::Priority::DETAIL);

  // declare the tuple-comparison functor
  SweepCmp cmp;

  // processing
  int iteration{};
  for(SimplexId i = 0; i < vertexNumber_; ++i) {

    this->printMsg("Starting simplifying iteration #" + std::to_string(i),
                   debug::Priority::DETAIL);

    for(int j = 0; j < 2; ++j) {

      bool isIncreasingOrder = !j;

      cmp.setIsIncreasingOrder(isIncreasingOrder);
      std::set<std::tuple<dataType, SimplexId, SimplexId>, decltype(cmp)>
        sweepFront(cmp);
      std::vector<bool> visitedVertices(vertexNumber_, false);
      std::vector<SimplexId> adjustmentSequence(vertexNumber_);

      // add the seeds
      if(isIncreasingOrder) {
        for(SimplexId k : authorizedMinima) {
          authorizedExtrema[k] = true;
          sweepFront.emplace(outputScalars[k], offsets[k], k);
          visitedVertices[k] = true;
        }
      } else {
        for(SimplexId k : authorizedMaxima) {
          authorizedExtrema[k] = true;
          sweepFront.emplace(outputScalars[k], offsets[k], k);
          visitedVertices[k] = true;
        }
      }

      // growth by neighborhood of the seeds
      SimplexId adjustmentPos = 0;
      do {
        auto front = sweepFront.begin();
        if(front == sweepFront.end())
          return -1;

        SimplexId vertexId = std::get<2>(*front);
        sweepFront.erase(front);

        SimplexId neighborNumber
          = triangulation.getVertexNeighborNumber(vertexId);
        for(SimplexId k = 0; k < neighborNumber; ++k) {
          SimplexId neighbor;
          triangulation.getVertexNeighbor(vertexId, k, neighbor);
          if(!visitedVertices[neighbor]) {
            sweepFront.emplace(
              outputScalars[neighbor], offsets[neighbor], neighbor);
            visitedVertices[neighbor] = true;
          }
        }
        adjustmentSequence[adjustmentPos] = vertexId;
        ++adjustmentPos;
      } while(!sweepFront.empty());

      // save offsets and rearrange outputScalars
      SimplexId offset = (isIncreasingOrder ? 0 : vertexNumber_ + 1);

      for(SimplexId k = 0; k < vertexNumber_; ++k) {

        if(isIncreasingOrder) {
          if(k
             && outputScalars[adjustmentSequence[k]]
                  <= outputScalars[adjustmentSequence[k - 1]])
            outputScalars[adjustmentSequence[k]]
              = outputScalars[adjustmentSequence[k - 1]];
          ++offset;
        } else {
          if(k
             && outputScalars[adjustmentSequence[k]]
                  >= outputScalars[adjustmentSequence[k - 1]])
            outputScalars[adjustmentSequence[k]]
              = outputScalars[adjustmentSequence[k - 1]];
          --offset;
        }
        offsets[adjustmentSequence[k]] = offset;
      }
    }

    // test convergence
    bool needForMoreIterations{false};
    std::vector<SimplexId> minima;
    std::vector<SimplexId> maxima;
    getCriticalPoints(outputScalars, offsets, minima, maxima, triangulation);

    if(maxima.size() > authorizedMaxima.size())
      needForMoreIterations = true;
    if(minima.size() > authorizedMinima.size())
      needForMoreIterations = true;

    this->printMsg(
      std::vector<std::vector<std::string>>{
        {"#Minima", std::to_string(minima.size())},
        {"#Maxima", std::to_string(maxima.size())},
      },
      debug::Priority::DETAIL);

    if(!needForMoreIterations) {
      for(SimplexId k : minima) {
        if(!authorizedExtrema[k]) {
          needForMoreIterations = true;
          break;
        }
      }
    }
    if(!needForMoreIterations) {
      for(SimplexId k : maxima) {
        if(!authorizedExtrema[k]) {
          needForMoreIterations = true;
          break;
        }
      }
    }

    // optional adding of perturbation
    if(addPerturbation_)
      addPerturbation<dataType>(outputScalars, offsets);

    ++iteration;
    if(!needForMoreIterations)
      break;
  }

  this->printMsg(
    "Simplified scalar field", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}
