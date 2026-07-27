#include <Statistics.h>

#include <algorithm>

using namespace std;
using namespace ttk;

template <typename T>
T Statistics::mean(const T *v, const int &dimension) {
  T mean = 0;
  for(int i = 0; i < dimension; ++i)
    mean += v[i];
  return mean / dimension;
}

template <typename T>
T Statistics::mean(const std::vector<T> &v) {
  return mean(v.data(), v.size());
}

template <typename T>
T Statistics::var(const T *v, const int &dimension) {
  return cov(v, v, dimension, dimension);
}

template <typename T>
T Statistics::var(const std::vector<T> &v) {
  return var(v.data(), v.size());
}

template <typename T>
T Statistics::cov(const T *v1,
                  const T *v2,
                  const int &dimension1,
                  const int &dimension2) {
  T cov = 0;
  T meanV1 = mean(v1, dimension1);
  T meanV2 = mean(v2, dimension2);
  for(int i = 0; i < dimension1; ++i)
    cov += (v1[i] - meanV1) * (v2[i] - meanV2);
  return cov / (dimension1 - 1);
}

template <typename T>
T Statistics::cov(const std::vector<T> &v1, const std::vector<T> &v2) {
  return cov(v1.data(), v2.data(), v1.size(), v2.size());
}

template <typename T>
T Statistics::corr(const T *v1,
                   const T *v2,
                   const int &dimension1,
                   const int &dimension2) {
  return cov(v1, v2, dimension1, dimension2)
         / (std::sqrt(var(v1, dimension1)) * std::sqrt(var(v2, dimension2)));
}

template <typename T>
T Statistics::corr(const std::vector<T> &v1, const std::vector<T> &v2) {
  return corr(v1.data(), v2.data(), v1.size(), v2.size());
}

#define STATISTICS_SPECIALIZE(TYPE)                                \
  template TYPE Statistics::mean<TYPE>(TYPE const *, int const &); \
  template TYPE Statistics::mean<TYPE>(std::vector<TYPE> const &); \
  template TYPE Statistics::var<TYPE>(TYPE const *, int const &);  \
  template TYPE Statistics::var<TYPE>(std::vector<TYPE> const &);  \
  template TYPE Statistics::cov<TYPE>(                             \
    TYPE const *, TYPE const *, int const &, int const &);         \
  template TYPE Statistics::cov<TYPE>(                             \
    std::vector<TYPE> const &, std::vector<TYPE> const &);         \
  template TYPE Statistics::corr<TYPE>(                            \
    TYPE const *, TYPE const *, int const &, int const &);         \
  template TYPE Statistics::corr<TYPE>(                            \
    std::vector<TYPE> const &, std::vector<TYPE> const &);

// explicit specializations for float and double
STATISTICS_SPECIALIZE(float);
STATISTICS_SPECIALIZE(double);
