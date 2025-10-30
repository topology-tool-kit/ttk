#pragma once

#include <RipsPersistenceDiagramUtils.h>

#include <boost/version.hpp>
#if ((BOOST_VERSION / 100) % 1000) >= 81
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#elif ((BOOST_VERSION / 100) % 1000) >= 36
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>
#else
#include <unordered_map>
#include <unordered_set>
#include <boost/container_hash/hash.hpp>
#endif

using namespace ttk::rpd;

namespace ttk::gph {
  using id_t = int;

  template <unsigned DIM>
  using PointD = std::conditional_t<DIM==0, std::vector<value_t>, std::array<value_t,DIM>>;

  template <unsigned DIM>
  using PointCloud = std::vector<PointD<DIM>>;

  using Facet = std::array<id_t, 3>;

  struct FiltratedFacet {
    Facet f;
    double d;
  };

  inline FiltratedFacet max(FiltratedFacet a, FiltratedFacet b) {
    if (a.d > b.d)
      return a;
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
  template <typename X> using HashSet = boost::unordered_flat_set<X>;
#elif ((BOOST_VERSION / 100) % 1000) >= 36
  template <typename X, typename Y> using HashMap = boost::unordered_map<X, Y>;
  template <typename X> using HashSet = boost::unordered_set<X>;
#else
  template <typename X, typename Y> using HashMap = std::unordered_map<X, Y, boost::hash<X>>;
  template <typename X> using HashSet = std::unordered_set<X, boost::hash<X>>;
#endif

}
