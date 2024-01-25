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
//#define PRINT_PERSISTENCE_PAIRS
//#define USE_ROBINHOOD_HASHMAP

#ifdef USE_ROBINHOOD_HASHMAP

#include "robin_hood.h"

template <class Key, class T, class H, class E>
using hash_map = robin_hood::unordered_map<Key, T, H, E>;
template <class Key> using hash = robin_hood::hash<Key>;

#else

template <class Key, class T, class H, class E> using hash_map = std::unordered_map<Key, T, H, E>;
template <class Key> using hash = std::hash<Key>;

#endif

typedef double value_t;
typedef int64_t index_t;
typedef uint16_t coefficient_t;

typedef std::vector<index_t> simplex_t;
typedef std::pair<simplex_t, value_t> simplex_diam_t;
typedef std::pair<simplex_diam_t,simplex_diam_t> pers_pair_t;

void Ripser(std::vector<std::vector<value_t> > points, value_t threshold, index_t dim_max, std::vector<std::vector<pers_pair_t> >& ph);
