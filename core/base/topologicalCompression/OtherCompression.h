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
int ttk::TopologicalCompression::ComputeTotalSizeForOther() {
  // Should return the number of bytes to be written on the output file
  // sizeof(char) = 1 (byte)
  // use sizeof(int), sizeof(double) to get the number of bytes of
  // the matching structures.
  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::computeOther() {
  // Code me
  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::compressForOther(int vertexNumber,
                                                  dataType *inputData,
                                                  dataType *outputData,
                                                  const double &tol) {
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

template <typename dataType>
int ttk::TopologicalCompression::WriteOtherTopology(FILE *fm) {
  this->printWrn("Writing Other index / topology.");
  // Code me
  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::WriteOtherGeometry(FILE *fm) {
  this->printWrn("Writing Other buffer / geometry.");
  // Code me
  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::ReadOtherTopology(FILE *fm) {
  this->printWrn("Reading Other index / topology.");
  // Code me
  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::ReadOtherGeometry(FILE *fm) {
  this->printWrn("Reading Other buffer / geometry.");
  // Code me
  return 0;
}
