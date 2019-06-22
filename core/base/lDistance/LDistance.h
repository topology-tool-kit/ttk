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

#ifndef _LDISTANCE_H
#define _LDISTANCE_H

// Standard.
#include <cmath>
#include <cstdlib>
#include <string>

// Base code.
#include <Wrapper.h>

namespace ttk {

  class LDistance : public Debug {

  public:
    LDistance();

    ~LDistance();

    template <class dataType>
    int execute(const std::string &distanceType);

    template <class dataType>
    int computeLn(dataType *input1,
                  dataType *input2,
                  dataType *output,
                  const int n,
                  const ttk::SimplexId vertexNumber);

    template <class dataType>
    int computeLinf(dataType *input1,
                    dataType *input2,
                    dataType *output,
                    const ttk::SimplexId vertexNumber);

    /// Pass a pointer to an input array representing a scalarfield.
    /// The expected format for the array is the following:
    /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
    /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
    /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
    /// The array is expected to be correctly allocated.
    /// \param data Pointer to the data array.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa setVertexNumber() and setDimensionNumber().
    inline int setInputDataPointer1(void *data) {
      inputData1_ = data;
      return 0;
    }

    inline int setInputDataPointer2(void *data) {
      inputData2_ = data;
      return 0;
    }

    /// Pass a pointer to an output array representing a scalar field.
    /// The expected format for the array is the following:
    /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
    /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
    /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
    /// The array is expected to be correctly allocated.
    /// \param data Pointer to the data array.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa setVertexNumber() and setDimensionNumber().
    inline int setOutputDataPointer(void *data) {
      outputData_ = data;
      return 0;
    }

    inline int setNumberOfPoints(ttk::SimplexId numberOfPoints) {
      numberOfPoints_ = numberOfPoints;
      return 0;
    }

    inline double getResult() {
      return result;
    }

    template <typename type>
    static type abs_diff(const type var1, const type var2) {
      return (var1 > var2) ? var1 - var2 : var2 - var1;
    }

  protected:
    void *inputData1_, *inputData2_, *outputData_;
    double result;
    ttk::SimplexId numberOfPoints_;
  };
} // namespace ttk

// template functions
template <class dataType>
int ttk::LDistance::execute(const std::string &distanceType) {

  Timer t;
  int status;

// Check variables consistency
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputData1_ || !inputData2_ || !outputData_)
    return -1;
#endif

  dataType *outputData = (dataType *)outputData_;
  dataType *inputData1 = (dataType *)inputData1_;
  dataType *inputData2 = (dataType *)inputData2_;

  ttk::SimplexId vertexNumber = numberOfPoints_;

  if(distanceType == "inf") {
    status = computeLinf(inputData1, inputData2, outputData, vertexNumber);
  } else {
    int n = stoi(distanceType);
    if(n < 1)
      return -4;

    status = computeLn(inputData1, inputData2, outputData, n, vertexNumber);
  }

  {
    std::stringstream msg;
    msg << "[LDistance] Data-set (" << vertexNumber << " points) processed in "
        << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return status;
}

template <class dataType>
int ttk::LDistance::computeLn(dataType *input1,
                              dataType *input2,
                              dataType *output,
                              const int n,
                              const ttk::SimplexId vertexNumber) {
  dataType sum = 0;

// Compute difference for each point.
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) reduction(+ : sum)
#endif
  for(ttk::SimplexId i = 0; i < vertexNumber; ++i) {
    const dataType diff = abs_diff<dataType>(input1[i], input2[i]);
    const dataType power = pow(diff, (double)n);

    // Careful: huge dataset + huge values
    // may exceed double capacity.
    sum += power;

    // Store difference.
    if(output)
      output[i] = power;
  }

  sum = pow(sum, 1.0 / (double)n);

  // Affect result.
  result = (double)sum;
  {
    std::stringstream msg;
    msg << "[LDistance] Distance: " << result << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <class dataType>
int ttk::LDistance::computeLinf(dataType *input1,
                                dataType *input2,
                                dataType *output,
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
  {
    std::stringstream msg;
    msg << "[LDistance] Distance: " << result << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // LDISTANCE_H
