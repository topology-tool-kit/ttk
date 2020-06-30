/// \ingroup base
/// \class ttk::DilateErode
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.07.2020
///
/// \brief TTK %dilateErode processing package.
///
/// %DilateErode is a TTK processing package that either a) dilates a specified
/// label by assigning the label of a corresponding vertex to all its neighbors,
/// or b) erodes a specified label by assigning to a corresponding vertex the
/// largest label among its neighbors.

#pragma once

// base code includes
#include <Debug.h>
#include <Triangulation.h>
#include <limits>

namespace ttk {

  class DilateErode : public virtual Debug {
  public:
    DilateErode() {
      this->setDebugMsgPrefix("DilateErode");
    }

    ~DilateErode() {
    }

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    };

    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int execute(
      // Output
      dataType *outputLabels,

      // Input
      const int &mode,
      const dataType *inputLabels,
      const dataType &pivotLabel,
      triangulationType *triangulation) const;
  };
} // namespace ttk

template <class dataType, class triangulationType = ttk::AbstractTriangulation>
int ttk::DilateErode::execute(
  // Output
  dataType *outputLabels,

  // Input
  const int &mode,
  const dataType *inputLabels,
  const dataType &pivotLabel,
  triangulationType *triangulation) const {

  Timer t;

  SimplexId nVertices = triangulation->getNumberOfVertices();

  std::string msg = std::string(mode == 0 ? "Dilating" : "Eroding") + " value "
                    + std::to_string(pivotLabel);

  this->printMsg(msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

  // NOTE: Directly dilating a vertex value to all its neighbors requires
  // parallel write locks, so instead focusing on vertices that need to
  // update their value optimizes parallel efficiency.
  if(mode == 0) {
// Dilate
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
    for(SimplexId i = 0; i < nVertices; i++) {
      // if current vertex value is not a dilated value
      if(inputLabels[i] != pivotLabel) {
        // check neighbors if they need to be dilated
        const SimplexId nNeighbors = triangulation->getVertexNeighborNumber(i);
        SimplexId nIndex;
        for(SimplexId n = 0; n < nNeighbors; n++) {
          triangulation->getVertexNeighbor(i, n, nIndex);
          if(inputLabels[nIndex] == pivotLabel) {
            outputLabels[i] = inputLabels[nIndex];
            break;
          }
        }
      }
    }
  } else {
// Erode
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
    for(SimplexId i = 0; i < nVertices; i++) {
      // if current vertex value needs to be eroded
      if(inputLabels[i] == pivotLabel) {
        // check neighbors if neighbors have a non-eroded label
        const SimplexId nNeighbors = triangulation->getVertexNeighborNumber(i);
        SimplexId nIndex;
        for(SimplexId n = 0; n < nNeighbors; n++) {
          triangulation->getVertexNeighbor(i, n, nIndex);
          if(inputLabels[nIndex] != pivotLabel) {
            outputLabels[i] = inputLabels[nIndex];
            break;
          }
        }
      }
    }
  }

  this->printMsg(msg, 1, t.getElapsedTime(), this->threadNumber_);

  return 1;
}