/// \ingroup base
/// \class ttk::TrackingFromPersistenceDiagrams
/// \author Maxime Soler <soler.maxime@total.com>
/// \date August 2018.

#ifndef _TRACKINGFROMF_H
#define _TRACKINGFROMF_H

// base code includes
#include <BottleneckDistance.h>
#include <PersistenceDiagram.h>
#include <Wrapper.h>

namespace ttk {

  class TrackingFromFields : public Debug {

  public:
    TrackingFromFields() {
    }

    ~TrackingFromFields() {
    }

    /// Execute the package.
    /// \return Returns 0 upon success, negative values otherwise.
    template <class dataType>
    int execute();

    template <typename dataType>
    int performDiagramComputation(
      int fieldNumber,
      std::vector<std::vector<diagramTuple>> &persistenceDiagrams,
      const ttk::Wrapper *wrapper);

    /// Pass a pointer to an input array representing a scalarfield.
    /// The array is expected to be correctly allocated. idx in
    /// [0,numberOfInputs_[ \param idx Index of the input scalar field. \param
    /// data Pointer to the data array. \return Returns 0 upon success, negative
    /// values otherwise. \sa setNumberOfInputs() and setVertexNumber().
    inline int setInputDataPointer(int idx, void *data) {
      if(idx < numberOfInputs_)
        inputData_[idx] = data;
      else
        return -1;
      return 0;
    }

    /// Set the number of input scalar fields
    /// \param numberOfInputs Number of input scalar fields.
    /// \return Returns 0 upon success, negative values otherwise
    inline int setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
      return 0;
    }

    inline int setTriangulation(ttk::Triangulation *t) {
      triangulation_ = t;
      return 0;
    }

    inline int setInputScalars(std::vector<void *> &is) {
      inputData_ = is;
      return 0;
    }

    inline int setInputOffsets(void *io) {
      inputOffsets_ = io;
      return 0;
    }

  protected:
    int numberOfInputs_;
    std::vector<void *> inputData_;
    void *inputOffsets_;
    ttk::Triangulation *triangulation_; // 1 triangulation for everyone
  };
} // namespace ttk

// template functions
template <class dataType>
int ttk::TrackingFromFields::execute() {
  ttk::Timer t;

  {
    std::stringstream msg;
    msg << "[TrackingFromFields] Data-set "
        << "processed in " << t.getElapsedTime() << " s. (" << threadNumber_
        << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int ttk::TrackingFromFields::performDiagramComputation(
  int fieldNumber,
  std::vector<std::vector<diagramTuple>> &persistenceDiagrams,
  const ttk::Wrapper *wrapper) {

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < fieldNumber; ++i) {
    ttk::PersistenceDiagram persistenceDiagram_;
    persistenceDiagram_.setWrapper(wrapper);
    persistenceDiagram_.setupTriangulation(triangulation_);
    persistenceDiagram_.setThreadNumber(1);
    // should have been done before

    std::vector<std::tuple<ttk::dcg::Cell, ttk::dcg::Cell>> dmt_pairs;
    persistenceDiagram_.setDMTPairs(&dmt_pairs);
    persistenceDiagram_.setInputScalars(inputData_[i]);
    persistenceDiagram_.setInputOffsets(inputOffsets_);
    persistenceDiagram_.setComputeSaddleConnectors(false);
    std::vector<std::tuple<int, CriticalType, int, CriticalType, dataType, int>>
      CTDiagram;

    persistenceDiagram_.setOutputCTDiagram(&CTDiagram);
    persistenceDiagram_.execute<dataType, int>();

    // Copy diagram into augmented diagram.
    persistenceDiagrams[i] = std::vector<diagramTuple>(CTDiagram.size());

    for(int j = 0; j < (int)CTDiagram.size(); ++j) {
      float p[3];
      float q[3];
      auto currentTuple = CTDiagram[j];
      const int a = std::get<0>(currentTuple);
      const int b = std::get<2>(currentTuple);
      triangulation_->getVertexPoint(a, p[0], p[1], p[2]);
      triangulation_->getVertexPoint(b, q[0], q[1], q[2]);
      const double sa = ((double *)inputData_[i])[a];
      const double sb = ((double *)inputData_[i])[b];
      diagramTuple dt
        = std::make_tuple(std::get<0>(currentTuple), std::get<1>(currentTuple),
                          std::get<2>(currentTuple), std::get<3>(currentTuple),
                          std::get<4>(currentTuple), std::get<5>(currentTuple),
                          sa, p[0], p[1], p[2], sb, q[0], q[1], q[2]);

      persistenceDiagrams[i][j] = dt;
    }
  }

  return 0;
}

#endif // _TRACKINGFROMP_H
