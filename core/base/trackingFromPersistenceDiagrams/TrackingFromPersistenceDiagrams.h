/// \ingroup base
/// \class ttk::TrackingFromPersistenceDiagrams
/// \author Maxime Soler <soler.maxime@total.com>
/// \date August 2018.
#pragma once

// base code includes
#include <BottleneckDistance.h>
#include <Geometry.h>
#include <set>

namespace ttk {

  using trackingTuple = std::tuple<int, int, std::vector<SimplexId>>;

  class TrackingFromPersistenceDiagrams : virtual public Debug {

  public:
    TrackingFromPersistenceDiagrams();

    /// Execute the package.
    /// \return Returns 0 upon success, negative values otherwise.
    int execute();

    int performSingleMatching(
      int i,
      std::vector<ttk::DiagramType> &inputPersistenceDiagrams,
      std::vector<std::vector<MatchingType>> &outputMatchings,
      const std::string &algorithm,
      const std::string &wasserstein,
      double tolerance,
      double px,
      double py,
      double pz,
      double ps,
      double pe);

    int
      performMatchings(int numInputs,
                       std::vector<ttk::DiagramType> &inputPersistenceDiagrams,
                       std::vector<std::vector<MatchingType>> &outputMatchings,
                       const std::string &algorithm,
                       const std::string &wasserstein,
                       double tolerance,
                       double px,
                       double py,
                       double pz,
                       double ps,
                       double pe);

    int performTracking(std::vector<ttk::DiagramType> &allDiagrams,
                        std::vector<std::vector<MatchingType>> &allMatchings,
                        std::vector<trackingTuple> &trackings);

    int performPostProcess(const std::vector<ttk::DiagramType> &allDiagrams,
                           const std::vector<trackingTuple> &trackings,
                           std::vector<std::set<int>> &trackingTupleToMerged,
                           const double postProcThresh);

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
    inline void setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
    }

  protected:
    int numberOfInputs_{};
    void **inputData_{};
  };
} // namespace ttk
