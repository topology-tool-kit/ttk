/// \ingroup base
/// \class ttk::BarycentricSubdivision
/// \author Pierre Guillou (pierre.guillou@lip6.fr
/// \date July 2019
///
/// \brief TTK %barycentricSubdivision processing package.
///
/// %BarycentricSubdivision is a TTK processing package that takes a scalar
/// field on the input and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkBarycentricSubdivision.cpp %for a usage example.

#pragma once

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {

  class BarycentricSubdivision : public Debug {

  public:
    BarycentricSubdivision() = default;

    ~BarycentricSubdivision() = default;

    template <class dataType>
    int execute(const int &argument) const;

    inline int setInputDataPointer(void *const data) {
      inputData_ = data;
      return 0;
    }

    inline int setOutputDataPointer(void *const data) {
      outputData_ = data;
      return 0;
    }

    inline int setupTriangulation(Triangulation *const triangulation) {
      inputTriangl_ = triangulation;
      if(inputTriangl_ != nullptr) {
        inputTriangl_->preprocessVertexNeighbors();
        inputTriangl_->preprocessEdges();
        inputTriangl_->preprocessTriangles();
      }

      return 0;
    }

  private:
    int subdiviseTriangulation();

    void *inputData_{}, *outputData_{};
    Triangulation *inputTriangl_{}, *outputTriangl_{};

    LongSimplexId *pointSet_{};
    std::vector<float> points_{};
    std::vector<LongSimplexId> cells_{};
  };
} // namespace ttk

// if the package is a pure template class, uncomment the following line
// #include                  <BarycentricSubdivision.cpp>

// template functions
template <class dataType>
int ttk::BarycentricSubdivision::execute(const int & /*argument*/) const {

  Timer t;

  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputTriangl_)
    return -1;
  if(!inputData_)
    return -2;
  if(!outputData_)
    return -3;
#endif

  dataType *outputData = (dataType *)outputData_;
  dataType *inputData = (dataType *)inputData_;

  SimplexId vertexNumber = inputTriangl_->getNumberOfVertices();

  // init the output -- to adapt
  for(SimplexId i = 0; i < vertexNumber; i++) {
    outputData[i] = inputData[i];
  }

  // the following open-mp processing is only relevant for embarrassingly
  // parallel algorithms (such as smoothing) -- to adapt
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < vertexNumber; i++) {
    // TODO-2
    // processing here!
    // end of TODO-2
  }

  {
    std::stringstream msg;
    msg << "[BarycentricSubdivision] Data-set (" << vertexNumber
        << " points) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}
