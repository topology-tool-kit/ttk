#pragma once

#include <RipsPersistenceDiagramUtils.h>

#include <boost/version.hpp>
#if ((BOOST_VERSION / 100) % 1000) >= 81
#include <boost/unordered/unordered_flat_map.hpp>
#define GPH_HASHMAP boost::unordered_flat_map
#elif ((BOOST_VERSION / 100) % 1000) >= 36
#include <boost/unordered/unordered_map.hpp>
#define GPH_HASHMAP boost::unordered_map
#else
#include <unordered_map>
#include <boost/container_hash/hash.hpp>
#define GPH_HASHMAP std::unordered_map
#endif

using namespace ttk::rpd;

namespace ttk::gph {
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

}
