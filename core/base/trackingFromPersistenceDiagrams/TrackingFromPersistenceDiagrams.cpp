#include <TrackingFromPersistenceDiagrams.h>

ttk::TrackingFromPersistenceDiagrams::TrackingFromPersistenceDiagrams() {
  inputData_ = nullptr;
  numberOfInputs_ = 0;
}

ttk::TrackingFromPersistenceDiagrams::~TrackingFromPersistenceDiagrams() {
  if(inputData_)
    free(inputData_);
}
