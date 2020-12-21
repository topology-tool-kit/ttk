#include <TrackingFromPersistenceDiagrams.h>

ttk::TrackingFromPersistenceDiagrams::TrackingFromPersistenceDiagrams() {
  inputData_ = nullptr;
  numberOfInputs_ = 0;
  this->setDebugMsgPrefix("TrackingFromPersistenceDiagrams");
}

ttk::TrackingFromPersistenceDiagrams::~TrackingFromPersistenceDiagrams() {
  if(inputData_)
    free(inputData_);
}
