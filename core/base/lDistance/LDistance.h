/// \ingroup base
/// \class ttk::LDistance
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 26/02/2017
///
/// \brief TTK %lDistance processing package.
///
/// %LDistance is a TTK processing package that takes a scalar field on the
/// input and produces a scalar field on the output.
///
/// \sa ttkLDistance.cpp for a usage example.

#pragma once

// Standard.
#include <cmath>
#include <cstdlib>
#include <string>

// Base code.
#include <Debug.h>
#include <Geometry.h>

namespace ttk {

  class LDistance : virtual public Debug {
  public:
    LDistance();

    template <class dataType>
    int execute(const dataType *const inputData1,
                const dataType *const inputData2,
                dataType *const outputData,
                const std::string &distanceType,
                const SimplexId vertexNumber);

    template <class dataType>
    int computeLn(const dataType *const input1,
                  const dataType *const input2,
                  dataType *const output,
                  const int n,
                  const SimplexId vertexNumber);

    template <class dataType>
    int computeLinf(const dataType *const input1,
                    const dataType *const input2,
                    dataType *const output,
                    const SimplexId vertexNumber);

    inline double getResult() {
      return result;
    }

    inline void setPrintRes(const bool data) {
      this->printRes = data;
    }

    template <typename type>
    static type abs_diff(const type var1, const type var2) {
      return (var1 > var2) ? var1 - var2 : var2 - var1;
    }

  protected:
    double result{};
    bool printRes{true};
  };
} // namespace ttk

// template functions
template <class dataType>
int ttk::LDistance::execute(const dataType *const inputData1,
                            const dataType *const inputData2,
                            dataType *const outputData,
                            const std::string &distanceType,
                            const SimplexId vertexNumber) {

  Timer t;
  int status;

// Check variables consistency
#ifndef TTK_ENABLE_KAMIKAZE
  if(inputData1 == nullptr || inputData2 == nullptr) {
    return -1;
  }
#endif

  if(distanceType == "inf") {
    status = computeLinf(inputData1, inputData2, outputData, vertexNumber);
  } else {
    int n = stoi(distanceType);
    if(n < 1)
      return -4;

    status = computeLn(inputData1, inputData2, outputData, n, vertexNumber);
  }

  if(this->printRes) {
    this->printMsg(
      "Data-set processed", 1.0, t.getElapsedTime(), this->threadNumber_);
  }

  return status;
}

template <class dataType>
int ttk::LDistance::computeLn(const dataType *const input1,
                              const dataType *const input2,
                              dataType *const output,
                              const int n,
                              const ttk::SimplexId vertexNumber) {
  dataType sum = 0;

// Compute difference for each point.
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) reduction(+ : sum)
#endif
  for(ttk::SimplexId i = 0; i < vertexNumber; ++i) {
    const dataType diff = abs_diff<dataType>(input1[i], input2[i]);
    const dataType power = Geometry::pow(diff, n);

    // Careful: huge dataset + huge values
    // may exceed double capacity.
    sum += power;

    // Store difference.
    if(output)
      output[i] = power;
  }

  sum = Geometry::pow(sum, 1.0 / (double)n);

  // Affect result.
  result = (double)sum;
  if(this->printRes) {
    this->printMsg("L" + std::to_string(n)
                   + "-distance: " + std::to_string(result));
  }

  return 0;
}

template <class dataType>
int ttk::LDistance::computeLinf(const dataType *const input1,
                                const dataType *const input2,
                                dataType *const output,
                                const ttk::SimplexId vertexNumber) {
  if(vertexNumber < 1)
    return 0;

  dataType maxValue = abs_diff<dataType>(input1[0], input2[0]);

// Compute difference for each point.
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) reduction(max : maxValue)
#endif
  for(ttk::SimplexId i = 1; i < vertexNumber; ++i) {
    const dataType iter = abs_diff<dataType>(input1[i], input2[i]);
    if(iter > maxValue)
      maxValue = iter;

    // Store absolute difference in output.
    if(output)
      output[i] = iter;
  }

  // Affect result.
  result = (double)maxValue;
  if(this->printRes) {
    this->printMsg("Linf-distance: " + std::to_string(result));
  }

  return 0;
}
