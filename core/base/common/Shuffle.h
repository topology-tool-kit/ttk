#include <cstdlib>
#include <vector>

namespace ttk {
  /**
   * @brief Platform-independant alternative to std::shuffle
   * implementing the Fisher-Yates shuffle algorithm
   *
   * @param[in,out] toShuffle Vector of elements to be shuffled
   * @param[in] rng Random number generator
   */
  template <typename T, typename U>
  void shuffle(std::vector<T> &toShuffle, U &&rng) {
    for(size_t i = toShuffle.size() - 1; i >= 1; i--) {
      const auto j = rng() % i;
      std::swap(toShuffle[i], toShuffle[j]);
    }
  }
} // namespace ttk
