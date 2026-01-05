#pragma once

#include <RipsPersistenceDiagramUtils.h>

#include <boost/version.hpp>
#if((BOOST_VERSION / 100) % 1000) >= 81
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#else
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>
#endif

#if((BOOST_VERSION / 100) % 1000) >= 84
#include <boost/unordered/concurrent_flat_map.hpp>
#define TTK_CONCURRENT_HASHTABLE_AVAILABLE
#endif

#ifdef TTK_ENABLE_TBB
#include <tbb/concurrent_vector.h>
#include <tbb/global_control.h>
#endif

#ifdef TTK_ENABLE_OPENMP
#include <omp.h>
#endif

#ifdef __cpp_lib_execution
#include <execution>
#endif

#if defined(TTK_ENABLE_OPENMP) and defined(TTK_ENABLE_TBB) \
  and defined(TTK_CONCURRENT_HASHTABLE_AVAILABLE)
#define TTK_GPH_PARALLEL
#endif

#include <dset.h>

namespace ttk::gph {
  using id_t = int;

  template <unsigned DIM>
  using PointD = std::conditional_t<DIM == 0,
                                    std::vector<rpd::value_t>,
                                    std::array<rpd::value_t, DIM>>;

  template <unsigned DIM>
  using PointCloud = std::vector<PointD<DIM>>;

#if((BOOST_VERSION / 100) % 1000) >= 81
  template <typename X, typename Y>
  using HashMap = boost::unordered_flat_map<X, Y>;
  template <typename X>
  using HashSet = boost::unordered_flat_set<X>;
#else
  template <typename X, typename Y>
  using HashMap = boost::unordered_map<X, Y>;
  template <typename X>
  using HashSet = boost::unordered_set<X>;
#endif

#if((BOOST_VERSION / 100) % 1000) >= 84
  template <typename X, typename Y>
  using ConcurrentHashMap = boost::concurrent_flat_map<X, Y>;
#endif

} // namespace ttk::gph
