#include <ScalarFieldSmoother.h>

ttk::ScalarFieldSmoother::ScalarFieldSmoother() {
  this->inputData_ = nullptr;
  this->outputData_ = nullptr;
  this->dimensionNumber_ = 1;
  this->mask_ = nullptr;

  this->setDebugMsgPrefix("ScalarFieldSmoother");
}

ttk::ScalarFieldSmoother::~ScalarFieldSmoother() {
}
