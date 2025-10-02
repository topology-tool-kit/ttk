#pragma once

#include <RipsPersistenceDiagramUtils.h>

#include <boost/version.hpp>
#if ((BOOST_VERSION / 100) % 1000) >= 81
#include <boost/unordered/unordered_flat_map.hpp>
#elif ((BOOST_VERSION / 100) % 1000) >= 36
#include <boost/unordered/unordered_map.hpp>
#else
#include <unordered_map>
#include <boost/container_hash/hash.hpp>
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

#if ((BOOST_VERSION / 100) % 1000) >= 81
  template <typename X, typename Y> using HashMap = boost::unordered_flat_map<X, Y>;
#elif ((BOOST_VERSION / 100) % 1000) >= 36
  template <typename X, typename Y> using HashMap = boost::unordered_map<X, Y>;
#else
  template <typename X, typename Y> using HashMap = std::unordered_map<X, Y, boost::hash<X>>;
#endif

}
