#pragma once

#ifdef TTK_ENABLE_OPENMP
#include <omp.h>
#endif // TTK_ENABLE_OPENMP

#if defined(__GNUC__) && !defined(__clang__)
#include <parallel/algorithm>
#endif

#ifdef TTK_ENABLE_OPENMP
namespace ttk {
  /**
   * @brief RAII wrapper to preserve the OpenMP runtime thread number
   */
  class ParallelGuard {
    int oldThreadNumber_;

  public:
    ParallelGuard(const int nThreads)
      : oldThreadNumber_{omp_get_num_threads()} {
      omp_set_num_threads(nThreads);
    }
    ~ParallelGuard() {
      omp_set_num_threads(oldThreadNumber_);
    }
  };
} // namespace ttk
#endif // TTK_ENABLE_OPENMP

/**
 * @brief Parallel sort macro
 *
 * Example use:
 * `TTK_PSORT(nthreads, container.begin(), container.end())`
 * or
 * `TTK_PSORT(nthreads, container.begin(), container.end(), cmp)`.
 */
#if defined(_GLIBCXX_PARALLEL_FEATURES_H) && defined(TTK_ENABLE_OPENMP)
#define TTK_PSORT(NTHREADS, ...)       \
  {                                    \
    ttk::ParallelGuard pg{NTHREADS};   \
    __gnu_parallel::sort(__VA_ARGS__); \
  }
#else
#define TTK_PSORT(NTHREADS, ...) std::sort(__VA_ARGS__);
#endif // _GLIBCXX_PARALLEL_FEATURES_H && TTK_ENABLE_OPENMP

namespace ttk {
  /**
   * @brief RAII wrapper around OpenMP lock
   */
  class Lock {
#ifdef TTK_ENABLE_OPENMP
  public:
    Lock() {
      omp_init_lock(&this->lock_);
    }
    ~Lock() {
      omp_destroy_lock(&this->lock_);
    }
    inline void lock() {
      omp_set_lock(&this->lock_);
    }
    inline void unlock() {
      omp_unset_lock(&this->lock_);
    }
    Lock(const Lock &) = delete;
    Lock(Lock &&) = delete;
    Lock &operator=(const Lock &) = delete;
    Lock &operator=(Lock &&) = delete;

  private:
    omp_lock_t lock_{};
#else
  public:
    inline void lock() {
    }
    inline void unlock() {
    }
#endif // TTK_ENABLE_OPENMP
  };
} // namespace ttk
