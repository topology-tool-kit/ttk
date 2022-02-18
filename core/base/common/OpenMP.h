#pragma once

#ifdef TTK_ENABLE_OPENMP
#include <omp.h>
#endif // TTK_ENABLE_OPENMP

#if defined(__GNUC__) && !defined(__clang__)
#include <parallel/algorithm>
#endif

/**
 * @brief Parallel sort macro
 *
 * Example use:
 * `TTK_PSORT(nthreads)(container.begin(), container.end())`
 */
#if defined(_GLIBCXX_PARALLEL_FEATURES_H) && defined(TTK_ENABLE_OPENMP)
#define TTK_PSORT(NTHREADS)      \
  omp_set_num_threads(NTHREADS); \
  __gnu_parallel::sort
#else
#define TTK_PSORT(NTHREADS) std::sort
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
