/// \ingroup base
/// \class ttk::FeatureTracking
/// \author Maxime Soler <soler.maxime@total.com>
/// \date August 2018.

#ifndef _FEATURETRACKING_H
#define _FEATURETRACKING_H

#ifndef diagramTuple
#define diagramTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::NodeType, ttk::ftm::idVertex, \
  ttk::ftm::NodeType, double, ttk::ftm::idVertex, \
  double, float, float, float, double, float, float, float>
#endif

#ifndef matchingTuple
#define matchingTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::idVertex, double>
#endif

#ifndef trackingTuple
#define trackingTuple std::tuple<int, int, std::vector<ttk::ftm::idVertex>>
#endif
// start ts, end ts or -1, list of indices for every ts

#ifndef BNodeType
#define BNodeType ttk::ftm::NodeType
#define BLocalMax ttk::ftm::NodeType::Local_maximum
#define BLocalMin ttk::ftm::NodeType::Local_minimum
#define BSaddle1  ttk::ftm::NodeType::Saddle1
#define BSaddle2  ttk::ftm::NodeType::Saddle2
#define BIdVertex ttk::ftm::idVertex
#endif

// base code includes
#include                  <Wrapper.h>
#include                  <PersistenceDiagram.h>

namespace ttk
{

  class FeatureTracking : public Debug {

    public:

      FeatureTracking();

      ~FeatureTracking();

      /// Execute the package.
      /// \return Returns 0 upon success, negative values otherwise.
      template <class dataType>
      int execute();

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
int ttk::FeatureTracking::execute()
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

#endif // FEATURETRACKING_H
