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

//#define USE_COEFFICIENTS
//#define INDICATE_PROGRESS
//#define USE_ROBINHOOD_HASHMAP

#ifdef USE_ROBINHOOD_HASHMAP

#include "robin_hood.h"

namespace Ripser {

template <class Key, class T, class H, class E>
using hash_map = robin_hood::unordered_map<Key, T, H, E>;
template <class Key> using hash = robin_hood::hash<Key>;

#else

namespace Ripser {

template <class Key, class T, class H, class E> using hash_map = std::unordered_map<Key, T, H, E>;
template <class Key> using hash = std::hash<Key>;

#endif

using value_t = double;
using index_t = int64_t;
using coefficient_t = uint16_t;

using simplex_t = std::vector<index_t>;
using simplex_diam_t = std::pair<simplex_t, value_t>;
using pers_pair_t = std::pair<simplex_diam_t, simplex_diam_t>;

void Ripser(std::vector<std::vector<value_t> > points, value_t threshold, index_t dim_max, bool distanceMatrix, std::vector<std::vector<pers_pair_t> >& ph);

}