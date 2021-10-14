/// \ingroup base
/// \class ttk::MorphologicalOperators
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.07.2020
///
/// \brief TTK %morphologicalOperators processing package.
///
/// %MorphologicalOperators is a TTK processing package that either a) dilates a
/// specified label by assigning the label of a corresponding vertex to all its
/// neighbors, or b) erodes a specified label by assigning to a corresponding
/// vertex the largest label among its neighbors.

#pragma once

// base code includes
#include <Debug.h>
#include <Triangulation.h>

#include <array>
#include <limits>

namespace ttk {

  class MorphologicalOperators : public virtual Debug {
  public:
    MorphologicalOperators() {
      this->setDebugMsgPrefix("MorphologicalOperators");
    }

    ~MorphologicalOperators() {
    }

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    template <class DT, class TT = ttk::AbstractTriangulation>
    int execute(
      // Output
      DT *outputLabels,

      // Input
      const int &mode,
      const int &iterations,
      const bool grayscale,
      const DT *inputLabels,
      const DT &pivotLabel,
      TT *triangulation) const;

    template <class DT, class TT = ttk::AbstractTriangulation>
    int performElementaryMorphoOp(
      // Output
      DT *outputLabels,

      // Input
      const int &mode,
      const int &iterations,
      const bool grayscale,
      const DT *inputLabels,
      const DT &pivotLabel,
      TT *triangulation) const;
  };
} // namespace ttk

template <class DT, class TT>
int ttk::MorphologicalOperators::execute(
  // Output
  DT *outputLabels,

  // Input
  const int &mode,
  const int &iterations,
  const bool grayscale,
  const DT *inputLabels,
  const DT &pivotLabel,
  TT *triangulation) const {

  if(mode < 2) {
    // erosion or dilation
    return this->performElementaryMorphoOp(outputLabels, mode, iterations,
                                           grayscale, inputLabels, pivotLabel,
                                           triangulation);
  } else {
    std::array<int, 2> ops{};
    if(mode == 2) {
      // morphological opening
      ops = {1, 0}; // erosion then dilation
    } else if(mode == 3) {
      // morphological closing
      ops = {0, 1}; // dilation then erosion
    } else {
      this->printErr("Invalid morphological operation requested");
      return 0;
    }

    std::vector<DT> tmp(triangulation->getNumberOfVertices());
    const auto status = this->performElementaryMorphoOp(
      tmp.data(), ops[0], iterations, grayscale, inputLabels, pivotLabel,
      triangulation);
    if(status != 1) {
      return status;
    }
    return this->performElementaryMorphoOp(outputLabels, ops[1], iterations,
                                           grayscale, tmp.data(), pivotLabel,
                                           triangulation);
  }
}

template <class DT, class TT>
int ttk::MorphologicalOperators::performElementaryMorphoOp(
  // Output
  DT *outputLabels,

  // Input
  const int &mode,
  const int &iterations,
  const bool grayscale,
  const DT *inputLabels,
  const DT &pivotLabel,
  TT *triangulation) const {

  SimplexId nVertices = triangulation->getNumberOfVertices();

  std::vector<DT> temp;
  if(iterations > 1) {
    Timer t;
    this->printMsg("Allocating temporary memory", 0, 0, this->threadNumber_,
                   debug::LineMode::REPLACE);
    temp.resize(nVertices);
    this->printMsg("Allocating temporary memory", 1, t.getElapsedTime(),
                   this->threadNumber_);
  }

  std::string msg = std::string(mode == 0 ? "Dilating " : "Eroding ")
                    + std::to_string(iterations) + "x value "
                    + std::to_string(pivotLabel);

  this->printMsg(msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

  Timer t;
  for(int it = 0; it < iterations; it++) {

    const DT *source = it == 0                      ? inputLabels
                       : (iterations + it) % 2 == 0 ? outputLabels
                                                    : temp.data();
    DT *target = (iterations + it) % 2 == 0 ? temp.data() : outputLabels;

    // NOTE: Directly dilating a vertex value to all its neighbors requires
    // parallel write locks, so instead focusing on vertices that need to
    // update their value optimizes parallel efficiency.
    if(!grayscale) {
      if(mode == 0) { // binary dilation
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for(SimplexId i = 0; i < nVertices; i++) {
          // first copy data
          target[i] = source[i];

          // if current vertex value is not a dilated value
          if(source[i] != pivotLabel) {
            // check neighbors if they need to be dilated
            const SimplexId nNeighbors
              = triangulation->getVertexNeighborNumber(i);
            SimplexId nIndex;
            for(SimplexId n = 0; n < nNeighbors; n++) {
              triangulation->getVertexNeighbor(i, n, nIndex);
              if(source[nIndex] == pivotLabel) {
                target[i] = source[nIndex];
                break;
              }
            }
          }
        }
      } else { // binary erosion
        const DT minLabel = std::numeric_limits<DT>::min();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for(SimplexId i = 0; i < nVertices; i++) {
          // first copy data
          target[i] = source[i];

          // if current vertex value needs to be eroded
          if(source[i] == pivotLabel) {
            // check neighbors if neighbors have a non-eroded label
            const SimplexId nNeighbors
              = triangulation->getVertexNeighborNumber(i);
            SimplexId nIndex;
            DT maxNeighborLabel = minLabel;
            for(SimplexId n = 0; n < nNeighbors; n++) {
              triangulation->getVertexNeighbor(i, n, nIndex);
              if(source[nIndex] != pivotLabel
                 && maxNeighborLabel < source[nIndex]) {
                maxNeighborLabel = source[nIndex];
              }
            }
            if(maxNeighborLabel != minLabel)
              target[i] = maxNeighborLabel;
          }
        }
      }
    } else {

      if(mode == 0) { // grayscale dilation
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(SimplexId i = 0; i < nVertices; i++) {
          // first copy data
          target[i] = source[i];
          const auto nNeighs = triangulation->getVertexNeighborNumber(i);
          for(SimplexId n = 0; n < nNeighs; n++) {
            SimplexId neigh{};
            triangulation->getVertexNeighbor(i, n, neigh);
            target[i] = std::max(target[i], source[neigh]);
          }
        }
      } else { // grayscale erosion
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(SimplexId i = 0; i < nVertices; i++) {
          // first copy data
          target[i] = source[i];
          const auto nNeighs = triangulation->getVertexNeighborNumber(i);
          for(SimplexId n = 0; n < nNeighs; n++) {
            SimplexId neigh{};
            triangulation->getVertexNeighbor(i, n, neigh);
            target[i] = std::min(target[i], source[neigh]);
          }
        }
      }
    }

    this->printMsg(msg, (float)it / (float)(iterations - 1), t.getElapsedTime(),
                   this->threadNumber_, debug::LineMode::REPLACE);
  }

  this->printMsg(msg, 1, t.getElapsedTime(), this->threadNumber_);

  return 1;
}
