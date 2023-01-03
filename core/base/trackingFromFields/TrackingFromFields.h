/// \ingroup base
/// \class ttk::TrackingFromPersistenceDiagrams
/// \author Maxime Soler <soler.maxime@total.com>
/// \date August 2018.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/timeTracking/">Time
///   tracking example</a>

#pragma once

// base code includes
#include <BottleneckDistance.h>
#include <PersistenceDiagram.h>
#include <Triangulation.h>

namespace ttk {

  class TrackingFromFields : virtual public Debug {

  public:
    TrackingFromFields() {
      this->setDebugMsgPrefix("TrackingFromFields");
    }

    /// Execute the package.
    /// \return Returns 0 upon success, negative values otherwise.
    // template <class dataType>
    // int execute();

    template <typename dataType,
              typename triangulationType = ttk::AbstractTriangulation>
    int performDiagramComputation(
      int fieldNumber,
      std::vector<ttk::DiagramType> &persistenceDiagrams,
      const triangulationType *triangulation);

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

    inline void
      preconditionTriangulation(AbstractTriangulation *triangulation) {
      ttk::PersistenceDiagram pd{};
      pd.preconditionTriangulation(triangulation);
    }

    inline void setInputScalars(std::vector<void *> &is) {
      inputData_ = is;
    }

    /**
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill every buffer in the @p io vector prior to any
     * computation (the VTK wrapper already includes a mechanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    inline void setInputOffsets(std::vector<SimplexId *> &io) {
      inputOffsets_ = io;
    }

  protected:
    int numberOfInputs_{0};
    std::vector<void *> inputData_{};
    std::vector<SimplexId *> inputOffsets_{};
  };
} // namespace ttk

template <typename dataType, typename triangulationType>
int ttk::TrackingFromFields::performDiagramComputation(
  int fieldNumber,
  std::vector<ttk::DiagramType> &persistenceDiagrams,
  const triangulationType *triangulation) {

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < fieldNumber; ++i) {
    ttk::PersistenceDiagram persistenceDiagram;
    persistenceDiagram.setThreadNumber(1);
    persistenceDiagram.execute(persistenceDiagrams[i],
                               (dataType *)(inputData_[i]), 0, inputOffsets_[i],
                               triangulation);

    // Augment diagram.
    for(auto &pair : persistenceDiagrams[i]) {
      triangulation->getVertexPoint(pair.birth.id, pair.birth.coords[0],
                                    pair.birth.coords[1], pair.birth.coords[2]);
      triangulation->getVertexPoint(pair.death.id, pair.death.coords[0],
                                    pair.death.coords[1], pair.death.coords[2]);
      pair.birth.sfValue = static_cast<double *>(inputData_[i])[pair.birth.id];
      pair.death.sfValue = static_cast<double *>(inputData_[i])[pair.death.id];
    }
  }

  return 0;
}
