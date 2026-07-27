//
// Created by max on 24/05/18.
//

// Files to be modified:
// ttkTopologicalCompressionReader.cpp
// ttkTopologicalCompressionWriter.cpp
// OtherCompression.h

#pragma once

#include <TopologicalCompression.h>

template <typename dataType>
int ttk::TopologicalCompression::compressForOther(
  int ttkNotUsed(vertexNumber),
  const dataType *const ttkNotUsed(inputData),
  const SimplexId *const ttkNotUsed(inputOffsets),
  dataType *ttkNotUsed(outputData),
  const double &ttkNotUsed(tol)) const {

  ttk::Timer t;
  // Code me

  this->printMsg(
    "Other computed", 1.0, t.getElapsedTime(), this->threadNumber_);
  t.reStart();

  // Code me

  this->printMsg(
    "Scalar field compressed", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}
