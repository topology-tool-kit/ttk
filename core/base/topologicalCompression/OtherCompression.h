//
// Created by max on 24/05/18.
//

// Files to be modified:
// ttkTopologicalCompressionReader.cpp
// ttkTopologicalCompressionWriter.cpp
// OtherCompression.h

#ifndef TTK_OTHERCOMPRESSION_H
#define TTK_OTHERCOMPRESSION_H

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

  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Other computed in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
    t.reStart();
  }

  // Code me

  {
    std::stringstream msg;
    msg << "[TopologicalCompression] Scalar field compressed in "
        << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::WriteOtherTopology(FILE *fm) {
  std::cout << "[TopologicalCompression] Writing Other index / topology."
            << std::endl;

  // Code me

  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::WriteOtherGeometry(FILE *fm) {
  std::cout
    << "[ttkTopologicalCompressionReader] Writing Other buffer / geometry."
    << std::endl;

  // Code me

  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::ReadOtherTopology(FILE *fm) {
  std::cout
    << "[ttkTopologicalCompressionReader] Reading Other index / topology."
    << std::endl;
  // Code me

  return 0;
}

template <typename dataType>
int ttk::TopologicalCompression::ReadOtherGeometry(FILE *fm) {
  std::cout
    << "[ttkTopologicalCompressionReader] Reading Other buffer / geometry."
    << std::endl;
  // Code me
  return 0;
}

#endif // TTK_OTHERCOMPRESSION_H
