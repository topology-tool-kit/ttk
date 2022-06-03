/// \ingroup base
/// \class ttk::ScalarFieldSmoother
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2014.
///
/// \brief TTK processing package for scalar field smoothing.
///
/// This class is a dummy example for the development of TTK classes. It
/// smooths an input scalar field by averaging the scalar values on the link
/// of each vertex.
///
/// \param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa ttkScalarFieldSmoother.cpp %for a usage example.

#pragma once

// base code includes
#include <Triangulation.h>

namespace ttk {

  class ScalarFieldSmoother : virtual public Debug {

  public:
    ScalarFieldSmoother();
    ~ScalarFieldSmoother() override;

    inline void setDimensionNumber(const int &dimensionNumber) {
      dimensionNumber_ = dimensionNumber;
    }

    inline void setInputDataPointer(void *data) {
      inputData_ = data;
    }

    inline void setOutputDataPointer(void *data) {
      outputData_ = data;
    }

    inline void setMaskDataPointer(const char *const mask) {
      mask_ = mask;
    }

    int preconditionTriangulation(AbstractTriangulation *triangulation) {
      // Pre-condition functions.
      if(triangulation) {
        triangulation->preconditionVertexNeighbors();
      }
      return 0;
    }

    template <class dataType, class triangulationType = AbstractTriangulation>
    int smooth(const triangulationType *triangulation,
               const int &numberOfIterations,
               const int *rankArray = nullptr,
               const SimplexId *globalIds = nullptr) const;

  protected:
    int dimensionNumber_{1};
    void *inputData_{nullptr}, *outputData_{nullptr};
    const char *mask_{nullptr};
  };

} // namespace ttk

// template functions
template <class dataType, class triangulationType>
int ttk::ScalarFieldSmoother::smooth(const triangulationType *triangulation,
                                     const int &numberOfIterations,
                                     const int *rankArray,
                                     const SimplexId *globalIds) const {

  Timer t;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation)
    return -1;
  if(!dimensionNumber_)
    return -2;
  if(!inputData_)
    return -3;
  if(!outputData_)
    return -4;
#endif
  bool useMPI = false;
  TTK_FORCE_USE(useMPI);
  TTK_FORCE_USE(rankArray);
  TTK_FORCE_USE(globalIds);

#if TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI() && rankArray != nullptr && globalIds != nullptr)
    useMPI = true;
#endif
  SimplexId vertexNumber = triangulation->getNumberOfVertices();

  std::vector<dataType> tmpData(vertexNumber * dimensionNumber_, 0);

  dataType *outputData = (dataType *)outputData_;
  dataType *inputData = (dataType *)inputData_;
  // init the output
  for(SimplexId i = 0; i < vertexNumber; i++) {
    for(int j = 0; j < dimensionNumber_; j++) {
      outputData[dimensionNumber_ * i + j]
        = inputData[dimensionNumber_ * i + j];
    }
  }

  printMsg("Smoothing " + std::to_string(vertexNumber) + " vertices", 0, 0,
           threadNumber_, ttk::debug::LineMode::REPLACE);

  int timeBuckets = 10;
  if(numberOfIterations < timeBuckets)
    timeBuckets = numberOfIterations;

  for(int it = 0; it < numberOfIterations; it++) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < vertexNumber; i++) {

      // avoid to process masked vertices
      if(mask_ != nullptr && mask_[i] == 0)
        continue;

      for(int j = 0; j < dimensionNumber_; j++) {
        const auto curr{dimensionNumber_ * i + j};
        tmpData[curr] = outputData[curr];

        const auto neighborNumber = triangulation->getVertexNeighborNumber(i);
        for(SimplexId k = 0; k < neighborNumber; k++) {
          SimplexId neighborId = -1;
          triangulation->getVertexNeighbor(i, k, neighborId);
          tmpData[curr] += outputData[dimensionNumber_ * (neighborId) + j];
        }
        tmpData[curr] /= static_cast<double>(neighborNumber + 1);
      }
    }

    if(numberOfIterations) {
      // assign the tmpData back to the output
      for(SimplexId i = 0; i < vertexNumber; i++) {
        for(int j = 0; j < dimensionNumber_; j++) {
          // only set value for unmasked points
          if(mask_ == nullptr || mask_[i] != 0) {
            outputData[dimensionNumber_ * i + j]
              = tmpData[dimensionNumber_ * i + j];
          }
        }
      }
    }
#if TTK_ENABLE_MPI
    if(useMPI) {
      // after each iteration we need to exchange the ghostcell values with our
      // neighbors
      std::unordered_map<SimplexId, SimplexId> map;
      triangulation->getVertexGlobalIdMap(map);
      exchangeGhostCells<dataType, SimplexId>(
        outputData, rankArray, globalIds, map, vertexNumber, MPI_COMM_WORLD);
    }
#endif

    if(debugLevel_ >= (int)(debug::Priority::INFO)) {
      if(!(it % ((numberOfIterations) / timeBuckets))) {
        printMsg("Smoothing " + std::to_string(vertexNumber) + " vertices",
                 (it / (float)numberOfIterations), t.getElapsedTime(),
                 threadNumber_, debug::LineMode::REPLACE);
      }
    }
  }

  printMsg("Smoothed " + std::to_string(vertexNumber) + " vertices", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}
