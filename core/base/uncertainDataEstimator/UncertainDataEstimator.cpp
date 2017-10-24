#include                  <UncertainDataEstimator.h>

UncertainDataEstimator::UncertainDataEstimator(){

  inputData_ = NULL;
  outputLowerBoundField_ = NULL;
  outputUpperBoundField_ = NULL;
  outputProbability_ = NULL;
  binValues_ = NULL;
  vertexNumber_ = 0;
  numberOfInputs_ = 0;
  binCount_ = 0;
}

UncertainDataEstimator::~UncertainDataEstimator(){
  if(inputData_)
    free(inputData_);
  if(binValues_)
    free(binValues_);
  if(outputProbability_)
    free(outputProbability_);
}

