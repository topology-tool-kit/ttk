/// \ingroup base
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.

#pragma once

#include <RipsPersistenceDiagramUtils.h>

namespace ripser {

  using value_t = double;
#if defined(TTK_ENABLE_RIPSER_128BITS_IDS) \
  && (defined(__GNUC__) || defined(__clang__))
  using index_t = __int128;
#else
  using index_t = int64_t;
#endif
  using coefficient_t = uint16_t;

  template <typename PersistenceType>
  void ripser(std::vector<std::vector<value_t>> points,
              PersistenceType &ph,
              value_t threshold,
              index_t dim_max,
              bool distanceMatrix,
              bool criticalEdgesOnly = true,
              bool infinitePairs = true,
              coefficient_t modulus = 2);

  template <typename PersistenceType>
  void ripser(float *data,
              int n,
              int dim,
              PersistenceType &ph,
              value_t threshold,
              index_t dim_max,
              bool criticalEdgesOnly = true,
              bool infinitePairs = true,
              coefficient_t modulus = 2) {

    std::vector<std::vector<value_t>> points(n);
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j < dim; ++j)
        points[i].push_back(data[dim * i + j]);
    }

    ripser(points, ph, threshold, dim_max, false, criticalEdgesOnly,
           infinitePairs, modulus);
  }

} // namespace ripser