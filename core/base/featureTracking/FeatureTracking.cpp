#include                  <FeatureTracking.h>

ttk::FeatureTracking::FeatureTracking()
{
  inputData_ = nullptr;
  numberOfInputs_ = 0;
}

ttk::FeatureTracking::~FeatureTracking()
{
  if (inputData_)
    free(inputData_);
}

