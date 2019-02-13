/// ingroup base
/// \class ttk::HarmonicFieldComputation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date February 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
///
///
/// \sa ttkHarmonicFieldComputation.cpp % for a usage example.

#pragma once

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>
#include <cmath>
#include <set>
#include <tuple>
#include <type_traits>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#endif // TTK_ENABLE_EIGEN

using std::cout;
using std::endl;
using std::size_t;
using std::stringstream;

namespace ttk {

class HarmonicFieldComputation : public Debug {

public:
  HarmonicFieldComputation();

  ~HarmonicFieldComputation();

  inline int setVertexNumber(SimplexId vertexNumber) {
    vertexNumber_ = vertexNumber;
    return 0;
  }
  inline int setConstraintNumber(SimplexId constraintNumber) {
    constraintNumber_ = constraintNumber;
    return 0;
  }
  inline int setUseCotanMethod(bool useCotanMethod) {
    useCotanMethod_ = useCotanMethod;
    return 0;
  }
  inline int setupTriangulation(Triangulation *triangulation) {
    triangulation_ = triangulation;
    if (triangulation_) {
      vertexNumber_ = triangulation_->getNumberOfVertices();
      triangulation_->preprocessVertexNeighbors();
    }
    return 0;
  }
  inline int setInputScalarFieldPointer(void *data) {
    inputScalarFieldPointer_ = data;
    return 0;
  }
  inline int setVertexIdentifierScalarFieldPointer(void *data) {
    vertexIdentifierScalarFieldPointer_ = data;
    return 0;
  }
  inline int setOutputScalarFieldPointer(void *data) {
    outputScalarFieldPointer_ = data;
    return 0;
  }

  template <typename scalarFieldType> int execute() const;

protected:
  size_t vertexNumber_;
  size_t constraintNumber_;
  bool useCotanMethod_;
  Triangulation *triangulation_;
  void *inputScalarFieldPointer_;
  void *outputScalarFieldPointer_;
  void *vertexIdentifierScalarFieldPointer_;
};
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <HarmonicFieldComputation.cpp>

// main routine
template <typename scalarFieldType>
int ttk::HarmonicFieldComputation::execute() const {

  // scalar field constraints vertices
  SimplexId *identifiers = static_cast<SimplexId *>(inputScalarFieldPointer_);
  // scalar field
  scalarFieldType *sf =
      static_cast<scalarFieldType *>(outputScalarFieldPointer_);

  Timer t;

  {
    stringstream msg;
    msg << "[HarmonicFieldComputation] Beginnning computation" << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  std::fill(sf, sf + vertexNumber_,
            std::numeric_limits<scalarFieldType>::max());

  std::set<SimplexId> identifiersSet;
  for (size_t i = 0; i < constraintNumber_; i++) {
    identifiersSet.insert(identifiers[i]);
  }
  std::vector<SimplexId> identifiersVec(identifiersSet.begin(),
                                        identifiersSet.end());

  {
    stringstream msg;
    msg << "[HarmonicFieldComputation] Ending computation after"
        << t.getElapsedTime() << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  return 0;
}
