#pragma once

#include <DataTypes.h>
#include <Debug.h>

#include <algorithm>
#include <vector>

#if defined(__GNUC__) && !defined(__clang__)
#include <parallel/algorithm>
#endif

namespace ttk {

  /**
   * @brief Sort vertices according to scalars disambiguated by offsets
   *
   * @param[in] nVerts number of vertices
   * @param[in] scalars array of size nVerts, main vertex comparator
   * @param[in] offsets array of size nVerts, disambiguate scalars on plateaux
   * @param[out] order array of size nVerts, computed order of vertices
   * @param[in] nThreads number of parallel threads
   */
  template <typename scalarType, typename idType>
  void sortVertices(const size_t nVerts,
                    const scalarType *const scalars,
                    const idType *const offsets,
                    SimplexId *const order,
                    const int nThreads) {

    // array of pre-sorted vertices
    std::vector<SimplexId> sortedVertices(nVerts);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < sortedVertices.size(); ++i) {
      sortedVertices[i] = i;
    }

#if defined(_GLIBCXX_PARALLEL_FEATURES_H) && defined(TTK_ENABLE_OPENMP)
#define PSORT(NTHREADS)          \
  omp_set_num_threads(NTHREADS); \
  __gnu_parallel::sort
#else
#define PSORT(NTHREADS) std::sort
#endif // _GLIBCXX_PARALLEL_FEATURES_H && TTK_ENABLE_OPENMP

    if(offsets != nullptr) {
      PSORT(nThreads)
      (sortedVertices.begin(), sortedVertices.end(),
       [&](const SimplexId a, const SimplexId b) {
         return (scalars[a] < scalars[b])
                || (scalars[a] == scalars[b] && offsets[a] < offsets[b]);
       });
    } else {
      PSORT(nThreads)
      (sortedVertices.begin(), sortedVertices.end(),
       [&](const SimplexId a, const SimplexId b) {
         return (scalars[a] < scalars[b])
                || (scalars[a] == scalars[b] && a < b);
       });
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < sortedVertices.size(); ++i) {
      order[sortedVertices[i]] = i;
    }
  }

  /**
   * @brief Precondition an order array to be consumed by the base layer API
   *
   * @param[in] nVerts number of vertices
   * @param[in] scalars pointer to scalar field buffer of size @p nVerts
   * @param[out] order pointer to pre-allocated order buffer of size @p nVerts
   * @param[in] nThreads number of threads to be used
   */
  template <typename scalarType>
  inline void preconditionOrderArray(const size_t nVerts,
                                     const scalarType *const scalars,
                                     SimplexId *const order,
                                     const int nThreads
                                     = ttk::globalThreadNumber_) {
    ttk::sortVertices(
      nVerts, scalars, static_cast<int *>(nullptr), order, nThreads);
  }
} // namespace ttk
