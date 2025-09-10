#pragma once

#include <array>
#include <limits>
#include <vector>

#include <boost/version.hpp>
#if ((BOOST_VERSION / 100) % 1000) >= 81
#include <boost/unordered/unordered_flat_map.hpp>
#else
#include <unordered_map>
#endif

#if defined ENABLE_TBB and defined PARALLEL_SORT
#include <execution>
#define SORT(begin, end, comp) std::sort(std::execution::par_unseq, begin, end, comp)
#else
#define SORT(begin, end, comp) std::sort(begin, end, comp)
#endif

#include <RipsPersistenceDiagramUtils.h>
using namespace ttk::rpd;

namespace gph {
  using id_t = int;

  template <unsigned DIM>
  using PointD = std::array<value_t,DIM>;

  template <unsigned DIM>
  using PointCloud = std::vector<std::array<value_t,DIM>>;

  using Facet = std::array<id_t, 3>;

  struct FiltratedFacet {
    Facet f;
    double d;
  };

  inline FiltratedFacet max(FiltratedFacet a, FiltratedFacet b) {
    if (a.d > b.d)
      return a;
    else
      return b;
  }

  struct FiltratedQuadFacet {
    Facet f;
    int c1;
    int c2;
    double d;
    double a;
  };

  using Generator1 = std::pair<std::vector<Edge>,  std::pair<value_t,value_t>>;
  using Generator2 = std::pair<std::vector<Facet>, std::pair<value_t,value_t>>;

}
