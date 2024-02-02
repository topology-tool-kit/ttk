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

namespace ripser {

  using value_t = double;
  using index_t = int64_t;
  using coefficient_t = uint16_t;

  using simplex_t = std::vector<index_t>;
  using simplex_diam_t = std::pair<simplex_t, value_t>;
  using pers_pair_t = std::pair<simplex_diam_t, simplex_diam_t>;

  void ripser(std::vector<std::vector<value_t>> points,
              value_t threshold,
              index_t dim_max,
              bool distanceMatrix,
              std::vector<std::vector<pers_pair_t>> &ph);

} // namespace ripser