/// \ingroup base
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date September 2022.
///
/// \brief Minimalist namespace that handles simple statistics computations
/// (mean, variance, correlation etc.).

#pragma once

#include <Debug.h>
#include <array>

namespace ttk {

  namespace Statistics {

    /// Compute the sample mean of a random variable (e.g. a vector)
    /// \param v vector containing the values.
    /// \param dimension Optional parameter that specifies the dimension of
    /// the vector (by default 3).
    /// \return Returns the mean of the vector.
    template <typename T>
    T mean(const T *v, const int &dimension = 3);

    /// Compute the sample mean of a random variable (e.g. a vector)
    /// \param v vector containing the values.
    /// \return Returns the mean of the vector.
    template <typename T>
    T mean(const std::vector<T> &v);

    /// Compute the sample varuance of a random variable (e.g. a vector)
    /// \param v vector containing the values.
    /// \param dimension Optional parameter that specifies the dimension of
    /// the vector (by default 3).
    /// \return Returns the variance of the vector.
    template <typename T>
    T var(const T *v, const int &dimension = 3);

    /// Compute the sample varuance of a random variable (e.g. a vector)
    /// \param v vector containing the values.
    /// \return Returns the variance of the vector.
    template <typename T>
    T var(const std::vector<T> &v);

    /// Compute the sample covariance of two random variables (e.g. two vectors)
    /// \param v1 vector containing the values of the first variable..
    /// \param v2 vector containing the values of the second variable.
    /// \param dimension1 Optional parameter that specifies the dimension of
    /// the first vector (by default 3).
    /// \param dimension2 Optional parameter that specifies the dimension of
    /// the second vector (by default 3).
    /// \return Returns the covariance of the vectors.
    template <typename T>
    T cov(const T *v1,
          const T *v2,
          const int &dimension1 = 3,
          const int &dimension2 = 3);

    /// Compute the sample covariance of two random variables (e.g. two vectors)
    /// \param v1 vector containing the values of the first variable..
    /// \param v2 vector containing the values of the second variable.
    /// \return Returns the covariance of the vectors.
    template <typename T>
    T cov(const std::vector<T> &v1, const std::vector<T> &v2);

    /// Compute the sample correlation of two random variables (e.g. two
    /// vectors) \param v1 vector containing the values of the first variable..
    /// \param v2 vector containing the values of the second variable.
    /// \param dimension1 Optional parameter that specifies the dimension of
    /// the first vector (by default 3).
    /// \param dimension2 Optional parameter that specifies the dimension of
    /// the second vector (by default 3).
    /// \return Returns the correlation of the vectors.
    template <typename T>
    T corr(const T *v1,
           const T *v2,
           const int &dimension1 = 3,
           const int &dimension2 = 3);

    /// Compute the sample correlation of two random variables (e.g. two
    /// vectors) \param v1 vector containing the values of the first variable..
    /// \param v2 vector containing the values of the second variable.
    /// \return Returns the correlation of the vectors.
    template <typename T>
    T corr(const std::vector<T> &v1, const std::vector<T> &v2);

  } // namespace Statistics
} // namespace ttk
