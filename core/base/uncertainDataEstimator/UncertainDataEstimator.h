/// \ingroup base
/// \class ttk::UncertainDataEstimator
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \date August 2016.
///
/// \brief TTK processing package that takes an input ensemble data set
/// (represented by a list of scalar fields) and which computes various
/// vertexwise statistics (PDF estimation, bounds, moments, etc.)
///
/// \sa ttkUncertainDataEstimator.cpp %for a usage example.

#pragma once

// base code includes
#include <Wrapper.h>

namespace ttk {

  template <class dataType>
  class PDFBounds : virtual public Debug {
  public:
    int evaluateRealization(const void *voidPointer) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!(numberOfVertices_ > 0)) {
        return -1; // Number of vertices not defined
      }
#endif
      const dataType *inputData
        = reinterpret_cast<const dataType *>(voidPointer);
      auto numberOfVertices = static_cast<size_t>(numberOfVertices_);
      /* Initialize if first call since a change */
      if(!(upperBound_.size() == numberOfVertices)
         || !(lowerBound_.size() == numberOfVertices)) {
        upperBound_.resize(numberOfVertices);
        lowerBound_.resize(numberOfVertices);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        for(size_t i = 0; i < numberOfVertices; i++) {
          upperBound_[i] = inputData[i];
          lowerBound_[i] = inputData[i];
        }
      } else { /* Update the two fields with the new input */
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        for(size_t i = 0; i < numberOfVertices; i++) {
          // Upper Bound
          if(inputData[i] > upperBound_[i]) {
            upperBound_[i] = inputData[i];
          }
          // Lower Bound
          if(inputData[i] < lowerBound_[i]) {
            lowerBound_[i] = inputData[i];
          }
        }
      }
      return 0;
    }

    std::pair<dataType, dataType> getRange() const {
      return {getRangeMin(), getRangeMax()};
    }

    dataType getRangeMax() const {
      if(upperBound_.size()) {
        dataType maxValue = upperBound_[0];
        for(size_t i = 1; i < upperBound_.size(); i++) {
          if(upperBound_[i] > maxValue) {
            maxValue = upperBound_[i];
          }
        }
        return maxValue;
      } else {
        return 0;
      }
    }

    dataType getRangeMin() const {
      if(lowerBound_.size()) {
        dataType minValue = lowerBound_[0];
        for(size_t i = 1; i < lowerBound_.size(); i++) {
          if(lowerBound_[i] < minValue) {
            minValue = lowerBound_[i];
          }
        }
        return minValue;
      } else {
        return 0;
      }
    }

    inline dataType *getLowerBoundPointer() {
      return lowerBound_.data();
    }

    inline dataType *getUpperBoundPointer() {
      return upperBound_.data();
    }

    inline void setNumberOfVertices(const SimplexId number) {
      numberOfVertices_ = number;
    }

  protected:
    SimplexId numberOfVertices_{0};
    std::vector<dataType> upperBound_{};
    std::vector<dataType> lowerBound_{};
  };

  class PDFHistograms : virtual public Debug {
  public:
    template <class dataType>
    int evaluateRealization(const dataType *inputData) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!(rangeMin_ < rangeMax_)) {
        return -1; // Range error
      }
      if(!(numberOfBins_ > 0)) {
        return -2; // Number of bins not defined
      }
      if(!(numberOfVertices_ > 0)) {
        return -3; // Number of vertices not defined
      }
#endif
      if(numberOfInputs_ == 0) {
        /* Initialize */
        probability_.resize(numberOfBins_);
        double dx
          = (rangeMax_ - rangeMin_) / static_cast<double>(numberOfBins_);
        for(size_t i = 0; i < static_cast<size_t>(numberOfBins_); i++) {
          probability_[i].resize(numberOfVertices_);
          binValue_[i] = rangeMin_ + (dx / 2.0) + (static_cast<double>(i) * dx);
        }
      }
      /* Add input datas */
      for(SimplexId i = 0; i < numberOfVertices_; i++) {
        int bin
          = static_cast<int>(floor((inputData[i] - rangeMin_) * numberOfBins_
                                   / (rangeMax_ - rangeMin_)));
        bin = (bin == numberOfBins_) ? numberOfBins_ - 1 : bin;
        probability_[bin][i] += 1.0;
      }
      numberOfInputs_++;
      return 0;
    }

    inline double *getBinFieldPointer(const int binId) {
      if(binId < numberOfBins_) {
        return probability_[binId].data();
      } else {
        return nullptr;
      }
    }

    void getVertexHistogram(const SimplexId vertexId,
                            std::vector<double> &histogram) const;

    void normalize();

    inline void setNumberOfBins(const int number) {
      numberOfBins_ = number;
    }

    inline void setNumberOfVertices(const SimplexId number) {
      numberOfVertices_ = number;
    }

    inline void setRange(const double min, const double max) {
      rangeMin_ = min;
      rangeMax_ = max;
    }

  protected:
    std::vector<double> binValue_{};
    std::vector<std::vector<double>> probability_{};
    int numberOfBins_{0};
    int numberOfInputs_{0};
    SimplexId numberOfVertices_{0};
    double rangeMin_{0.0};
    double rangeMax_{0.0};
    // std::vector<int> selection_; // TODO : selection support
  };

  class UncertainDataEstimator : virtual public Debug {
  public:
    UncertainDataEstimator();

    /// Execute the package.
    /// \return Returns 0 upon success, negative values otherwise.
    template <class dataType>
    int execute();

    /// Pass a pointer to an input array representing a scalarfield.
    /// The array is expected to be correctly allocated. idx in
    /// [0,numberOfInputs_[ \param idx Index of the input scalar field. \param
    /// data Pointer to the data array. \return Returns 0 upon success, negative
    /// values otherwise. \sa setNumberOfInputs() and setVertexNumber().
    inline int setInputDataPointer(const int idx, void *const data) {
      if(idx < numberOfInputs_) {
        inputData_[idx] = data;
      } else {
        return -1;
      }
      return 0;
    }

    /// Pass a pointer to an output array representing the lower bound scalar
    /// field. The array is expected to be correctly allocated. \param data
    /// Pointer to the data array. \return Returns 0 upon success, negative
    /// values otherwise. \sa setVertexNumber()
    inline void setOutputLowerBoundField(void *const data) {
      outputLowerBoundField_ = data;
    }

    /// Pass a pointer to an output array representing the upper bound scalar
    /// field. The array is expected to be correctly allocated. \param data
    /// Pointer to the data array. \return Returns 0 upon success, negative
    /// values otherwise. \sa setVertexNumber()
    inline void setOutputUpperBoundField(void *const data) {
      outputUpperBoundField_ = data;
    }

    inline void setOutputProbability(const int idx, double *const data) {
      if(idx < BinCount) {
        outputProbability_[idx] = data;
      }
    }

    inline void setOutputMeanField(void *const data) {
      outputMeanField_ = data;
    }

    inline void setComputeLowerBound(const bool state) {
      ComputeLowerBound = state;
    }

    inline void setComputeUpperBound(const bool state) {
      ComputeUpperBound = state;
    }

    /// Set the number of vertices in the scalar field.
    /// \param vertexNumber Number of vertices in the data-set.
    /// \return Returns 0 upon success, negative values otherwise.
    inline void setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
    }

    inline void setBinCount(const int binCount) {
      BinCount = binCount;
      outputProbability_.clear();
      outputProbability_.resize(binCount, nullptr);
      binValues_.clear();
      binValues_.resize(binCount);
    }

    /// Set the number of input scalar fields
    /// \param numberOfInputs Number of input scalar fields.
    /// \return Returns 0 upon success, negative values otherwise
    inline void setNumberOfInputs(const int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
      inputData_.clear();
      inputData_.resize(numberOfInputs, nullptr);
    }

    inline double getBinValue(int b) {
      if(b < BinCount)
        return binValues_[b];
      return 0.0;
    }

  protected:
    SimplexId vertexNumber_{0};
    int numberOfInputs_{0};
    int BinCount{0};
    std::vector<double> binValues_{};
    bool ComputeLowerBound{};
    bool ComputeUpperBound{};
    std::vector<void *> inputData_{};
    void *outputLowerBoundField_{};
    void *outputUpperBoundField_{};
    std::vector<double *> outputProbability_{};
    void *outputMeanField_{};
  };
} // namespace ttk

template <class dataType>
int ttk::UncertainDataEstimator::execute() {

  Timer t;

  // Check the consistency of the variables
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfInputs_)
    return -1;
  if(!vertexNumber_)
    return -2;

  for(int i = 0; i < numberOfInputs_; i++) {
    if(!inputData_[i])
      return -4;
  }
  if(!outputLowerBoundField_)
    return -5;
  if(!outputUpperBoundField_)
    return -6;
#endif

  SimplexId count = 0;

  // Pointers type casting
  dataType *outputLowerBoundField = (dataType *)outputLowerBoundField_;
  dataType *outputUpperBoundField = (dataType *)outputUpperBoundField_;
  dataType **inputData = (dataType **)inputData_.data();
  double *outputMeanField = static_cast<double *>(outputMeanField_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId v = 0; v < (SimplexId)vertexNumber_; v++) {

    // Avoid any processing if the abort signal is sent
    if((!wrapper_) || ((wrapper_) && (!wrapper_->needsToAbort()))) {

      // For the lower bound scalar field
      if(ComputeLowerBound) {
        // Initialisation : values of the first input
        outputLowerBoundField[v] = inputData[0][v];
        // Loop over the inputs
        for(int inp = 1; inp < numberOfInputs_; inp++) {
          // Minimum value
          if(ComputeLowerBound)
            if(inputData[inp][v] < outputLowerBoundField[v])
              outputLowerBoundField[v] = inputData[inp][v];
        }
      }

      // For the upper bound scalar field
      if(ComputeUpperBound) {
        // Initialisation : values of the first input
        outputUpperBoundField[v] = inputData[0][v];
        // Loop over the inputs
        for(int inp = 1; inp < numberOfInputs_; inp++) {
          // Maximum value
          if(ComputeUpperBound)
            if(inputData[inp][v] > outputUpperBoundField[v])
              outputUpperBoundField[v] = inputData[inp][v];
        }
      }

      // Update the progress bar of the wrapping code -- to adapt
      if(debugLevel_ > static_cast<int>(debug::Priority::DETAIL)) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
        {
          if((wrapper_) && (!(count % ((vertexNumber_) / 10)))) {
            wrapper_->updateProgress((count + 1.0) / vertexNumber_);
          }

          count++;
        }
      }
    }
  }

  // Histogram
  if(ComputeUpperBound && ComputeLowerBound) {
    // Range
    double range[2];
    range[0] = outputLowerBoundField[0];
    range[1] = outputUpperBoundField[0];

    for(SimplexId v = 0; v < vertexNumber_; v++) {
      if(outputLowerBoundField[v] < range[0])
        range[0] = outputLowerBoundField[v];
      if(outputUpperBoundField[v] > range[1])
        range[1] = outputUpperBoundField[v];
    }

    // Interval between bins
    double dx = (range[1] - range[0]) / (double)BinCount;

    // Bin values
    for(int b = 0; b < BinCount; b++) {
      binValues_[b] = range[0] + (dx / 2.0) + (double)b * dx;
    }

    int idx;
    double increment = 1.0 / (double)numberOfInputs_;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for private(idx) num_threads(threadNumber_)
#endif
    for(SimplexId v = 0; v < vertexNumber_; v++) {
      for(int i = 0; i < numberOfInputs_; i++) {
        idx = (int)floor((inputData[i][v] - range[0]) * BinCount
                         / (range[1] - range[0]));
        idx = (idx == BinCount) ? BinCount - 1 : idx;
        outputProbability_[idx][v] += increment;
      }
    }
  }

  // Mean field
  for(SimplexId v = 0; v < vertexNumber_; v++) {
    double sum = 0.0;
    for(int i = 0; i < numberOfInputs_; i++) {
      sum += static_cast<double>(inputData[i][v]);
    }
    outputMeanField[v] = sum / static_cast<double>(numberOfInputs_);
  }

  this->printMsg(std::vector<std::vector<std::string>>{
    {"#Vertices", std::to_string(vertexNumber_)}});
  this->printMsg(
    "Data-set processed", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}
