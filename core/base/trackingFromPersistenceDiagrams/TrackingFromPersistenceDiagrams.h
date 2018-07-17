/// \ingroup base
/// \class ttk::TrackingFromPersistenceDiagrams
/// \author Maxime Soler <soler.maxime@total.com>
/// \date August 2018.

#ifndef _TRACKINGFROMP_H
#define _TRACKINGFROMP_H

// base code includes
#include                  <Wrapper.h>
#include                  <PersistenceDiagram.h>
#include                  <BottleneckDistance.h>

namespace ttk
{

  class TrackingFromFields : public Debug {

    using dataType = double;

    public:

      TrackingFromFields();

      ~TrackingFromFields();

      /// Execute the package.
      /// \return Returns 0 upon success, negative values otherwise.
      template <class dataType>
      int execute();

      int performSingleMatching(
        int i,
        std::vector<std::vector<diagramTuple>*>* inputPersistenceDiagrams,
        std::vector<std::vector<matchingTuple>*>* outputMatchings,
        std::string algorithm,
        std::string wasserstein,
        double tolerance,
        bool is3D,
        double alpha,
        double px, double py, double pz, double ps, double pe,
        const ttk::Wrapper *wrapper);

      int performMatchings(
        int numInputs,
        std::vector<std::vector<diagramTuple>*>* inputPersistenceDiagrams,
        std::vector<std::vector<matchingTuple>*>* outputMatchings,
        std::string algorithm,
        std::string wasserstein,
        double tolerance,
        bool is3D,
        double alpha,
        double px, double py, double pz, double ps, double pe,
        const ttk::Wrapper *wrapper);

      int performTracking(
        std::vector<std::vector<diagramTuple>*>* allDiagrams,
        std::vector<std::vector<matchingTuple>*>* allMatchings,
        std::vector<trackingTuple>* trackings);

      int performPostProcess(
        std::vector<std::vector<diagramTuple>*>* allDiagrams,
        std::vector<trackingTuple>* trackings,
        std::vector<std::set<int>>* trackingTupleToMerged,
        double postProcThresh);

      /// Pass a pointer to an input array representing a scalarfield.
      /// The array is expected to be correctly allocated. idx in [0,numberOfInputs_[
      /// \param idx Index of the input scalar field.
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa setNumberOfInputs() and setVertexNumber().
      inline int setInputDataPointer(int idx, void *data){
        if (idx < numberOfInputs_)
          inputData_[idx] = data;
        else
          return -1;
        return 0;
      }

      /// Set the number of input scalar fields
      /// \param numberOfInputs Number of input scalar fields.
      /// \return Returns 0 upon success, negative values otherwise
      inline int setNumberOfInputs(int numberOfInputs){
        numberOfInputs_ = numberOfInputs;
        return 0;
      }

    protected:

      int                   numberOfInputs_;
      void                  **inputData_; //TODO : vector<void*>
  };
}

// template functions
template <class dataType>
int ttk::TrackingFromFields::execute()
{
  ttk::Timer t;

  // Check the consistency of the variables
#ifndef TTK_ENABLE_KAMIKAZE
  if (!numberOfInputs_)
    return -1;
  if (!inputData_)
    return -3;

  for (int i = 0; i < numberOfInputs_; i++) {
    if (!inputData_[i])
      return -4;
  }
#endif

  {
    std::stringstream msg;
    msg << "[TrackingFromPersistenceDiagrams] Data-set "
      << "processed in " << t.getElapsedTime()
      << " s. (" << threadNumber_ << " thread(s))."
      << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // _TRACKINGFROMP_H
