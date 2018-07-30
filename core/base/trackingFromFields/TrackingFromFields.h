/// \ingroup base
/// \class ttk::TrackingFromPersistenceDiagrams
/// \author Maxime Soler <soler.maxime@total.com>
/// \date August 2018.

#ifndef _TRACKINGFROMF_H
#define _TRACKINGFROMF_H

// base code includes
#include                  <Wrapper.h>
#include                  <PersistenceDiagram.h>
#include                  <BottleneckDistance.h>

namespace ttk
{

  class TrackingFromFields : public Debug {

    using dataType = double;

    public:

      TrackingFromFields() {}

      ~TrackingFromFields() {}

      /// Execute the package.
      /// \return Returns 0 upon success, negative values otherwise.
      template <class dataType>
      int execute();

      int performDiagramComputation(
        int fieldNumber,
        std::vector<std::vector<diagramTuple>>& persistenceDiagrams,
        const ttk::Wrapper *wrapper);

      /// Pass a pointer to an input array representing a scalarfield.
      /// The array is expected to be correctly allocated. idx in [0,numberOfInputs_[
      /// \param idx Index of the input scalar field.
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa setNumberOfInputs() and setVertexNumber().
      inline int setInputDataPointer(int idx, void *data) {
        if (idx < numberOfInputs_)
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

      inline int setInputScalars(std::vector<void*> &is) {
        inputData_ = is;
        return 0;
      }

      inline int setInputOffsets(void* io) {
        inputOffsets_ = io;
        return 0;
      }

    protected:

      int                   numberOfInputs_;
      std::vector<void*>    inputData_;
      void*                 inputOffsets_;
      ttk::Triangulation    *triangulation_; // 1 triangulation for everyone
  };
}

// template functions
template <class dataType>
int ttk::TrackingFromFields::execute()
{
  ttk::Timer t;

  {
    std::stringstream msg;
    msg << "[TrackingFromFields] Data-set "
      << "processed in " << t.getElapsedTime()
      << " s. (" << threadNumber_ << " thread(s))."
      << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // _TRACKINGFROMP_H
