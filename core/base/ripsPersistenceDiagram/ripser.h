/// \ingroup base
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.

#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>

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
              coefficient_t modulus = 2);

  template <typename PersistenceType>
  void ripser(float *data, int n, int dim,
              PersistenceType &ph,
              value_t threshold,
              index_t dim_max,
              bool distanceMatrix,
              bool criticalEdgesOnly = true,
              coefficient_t modulus = 2);

} // namespace ripser