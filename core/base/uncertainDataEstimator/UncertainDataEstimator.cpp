#include <UncertainDataEstimator.h>

ttk::UncertainDataEstimator::UncertainDataEstimator() {
  this->setDebugMsgPrefix("UncertainDataEstimator");
}

void ttk::PDFHistograms::getVertexHistogram(
  const ttk::SimplexId vertexId, std::vector<double> &histogram) const {

  histogram.resize(numberOfBins_);
  if(vertexId < numberOfVertices_) {
#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp parallel for num_threads(threadNumber_)
#else
#pragma omp parallel for num_threads(threadNumber_) \
  schedule(static, numberOfBins_ / threadNumber_)
#endif
#endif
    for(int i = 0; i < (int)numberOfBins_; i++) {
      if((SimplexId)probability_[i].size() == numberOfVertices_) {
        histogram[i] = probability_[i][vertexId];
      } else {
        histogram[i] = 0.0;
      }
    }
  } else {
    fill(histogram.begin(), histogram.end(), 0.0);
  }
}

void ttk::PDFHistograms::normalize() {
  const double normalization = 1.0 / static_cast<double>(numberOfInputs_);
#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp parallel for num_threads(threadNumber_)
#else
#pragma omp parallel for num_threads(threadNumber_) collapse(2) \
  schedule(static, (numberOfBins_ * numberOfVertices_) / threadNumber_)
#endif
#endif
  for(int i = 0; i < numberOfBins_; i++) {
    for(SimplexId j = 0; j < numberOfVertices_; j++) {
      probability_[i][j] *= normalization;
    }
  }
}
