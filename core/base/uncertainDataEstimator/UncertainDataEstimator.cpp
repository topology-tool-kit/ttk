#include <UncertainDataEstimator.h>

using namespace std;
using namespace ttk;

UncertainDataEstimator::UncertainDataEstimator() {

  outputLowerBoundField_ = NULL;
  outputUpperBoundField_ = NULL;
  vertexNumber_ = 0;
  numberOfInputs_ = 0;
  binCount_ = 0;
}

UncertainDataEstimator::~UncertainDataEstimator() {
}
