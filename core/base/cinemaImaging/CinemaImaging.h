#pragma once

#include <Debug.h>

#include <limits>
#include <vector>

namespace ttk {
  class CinemaImaging : virtual public Debug {
  public:
    static const unsigned int INVALID_ID{
      std::numeric_limits<unsigned int>::max()};
    template <typename DT, typename IT>
    int interpolateArray(DT *outputArray,

                         const unsigned int *primitiveIds,
                         const float *barycentricCoordinates,
                         const IT *connectivityList,

                         const DT *inputArray,
                         const size_t &nTuples,
                         const size_t &nComponents = 1,
                         const DT &missingValue
                         = std::numeric_limits<DT>::has_quiet_NaN
                             ? std::numeric_limits<DT>::quiet_NaN()
                             : std::numeric_limits<DT>::max()) const;

    template <typename DT>
    int lookupArray(DT *outputArray,

                    const unsigned int *primitiveIds,

                    const DT *inputArray,
                    const size_t &nTuples,
                    const size_t &nComponents = 1,
                    const DT &missingValue
                    = std::numeric_limits<DT>::has_quiet_NaN
                        ? std::numeric_limits<DT>::quiet_NaN()
                        : std::numeric_limits<DT>::max()) const;
  };
} // namespace ttk

template <typename DT, typename IT>
int ttk::CinemaImaging::interpolateArray(DT *outputArray,

                                         const unsigned int *primitiveIds,
                                         const float *barycentricCoordinates,
                                         const IT *connectivityList,

                                         const DT *inputArray,
                                         const size_t &nTuples,
                                         const size_t &nComponents,
                                         const DT &missingValue) const {

  if(nComponents != 1)
    return 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
  for(size_t i = 0; i < nTuples; i++) {
    const unsigned int &cellId = primitiveIds[i];
    if(cellId == CinemaImaging::INVALID_ID) {
      outputArray[i] = missingValue;
      continue;
    }

    const size_t cellIndex = cellId * 3;
    const IT &v0 = connectivityList[cellIndex + 0];
    const IT &v1 = connectivityList[cellIndex + 1];
    const IT &v2 = connectivityList[cellIndex + 2];

    const size_t bcIndex = i * 2;
    const float &u = barycentricCoordinates[bcIndex + 0];
    const float &v = barycentricCoordinates[bcIndex + 1];
    const float w = 1 - u - v;

    outputArray[i]
      = w * inputArray[v0] + u * inputArray[v1] + v * inputArray[v2];
  }

  return 1;
}

template <typename DT>
int ttk::CinemaImaging::lookupArray(DT *outputArray,

                                    const unsigned int *primitiveIds,

                                    const DT *inputArray,
                                    const size_t &nTuples,
                                    const size_t &nComponents,
                                    const DT &missingValue) const {

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
  for(size_t i = 0; i < nTuples; i++) {
    size_t outputOffset = i * nComponents;
    const unsigned int &cellId = primitiveIds[i];

    if(cellId == CinemaImaging::INVALID_ID) {
      for(size_t j = 0; j < nComponents; j++)
        outputArray[outputOffset++] = missingValue;
      continue;
    } else {
      size_t inputOffset = cellId * nComponents;
      for(size_t j = 0; j < nComponents; j++)
        outputArray[outputOffset++] = inputArray[inputOffset++];
    }
  }

  return 1;
}
