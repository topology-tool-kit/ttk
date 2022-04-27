/// \ingroup base
/// \class ttk::DistanceField
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date March 2016
///
/// \brief TTK processing package for distance field computation on PL
/// manifolds.
///
/// This package takes a list of sources (a set of points with their global
/// identifiers attached to them) and produces a distance field to the closest
/// source.
///
/// \b Related \b publication \n
/// "A note on two problems in connexion with graphs" \n
/// Edsger W. Dijkstra \n
/// Numerische Mathematik, 1959.
///
/// \sa ttkDistanceField.cpp %for a usage example.
#pragma once

// base code includes
#include "Triangulation.h"
#include <Dijkstra.h>

// std includes
#include <limits>
#include <set>
#include <string>

namespace ttk {

  class DistanceField : virtual public Debug {
  public:
    DistanceField();

    template <typename dataType>
    dataType getDistance(const SimplexId a, const SimplexId b) const;

    template <typename dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int execute(const triangulationType *) const;

    inline void setVertexNumber(SimplexId vertexNumber) {
      vertexNumber_ = vertexNumber;
    }

    inline void setSourceNumber(SimplexId sourceNumber) {
      sourceNumber_ = sourceNumber;
    }

    inline int preconditionTriangulation(AbstractTriangulation *triangulation) {
      return triangulation->preconditionVertexNeighbors();
    }

    inline void setVertexIdentifierScalarFieldPointer(SimplexId *const data) {
      vertexIdentifierScalarFieldPointer_ = data;
    }

    inline void setOutputScalarFieldPointer(void *data) {
      outputScalarFieldPointer_ = data;
    }

    inline void setOutputIdentifiers(void *data) {
      outputIdentifiers_ = data;
    }

    inline void setOutputSegmentation(void *data) {
      outputSegmentation_ = data;
    }

  protected:
    SimplexId vertexNumber_{};
    SimplexId sourceNumber_{};
    SimplexId *vertexIdentifierScalarFieldPointer_{};
    void *outputScalarFieldPointer_{};
    void *outputIdentifiers_{};
    void *outputSegmentation_{};
  };
} // namespace ttk

template <typename dataType, class triangulationType>
int ttk::DistanceField::execute(const triangulationType *triangulation_) const {

  // start global timer
  ttk::Timer globalTimer;

  const SimplexId *identifiers = vertexIdentifierScalarFieldPointer_;
  dataType *dist = static_cast<dataType *>(outputScalarFieldPointer_);
  SimplexId *origin = static_cast<SimplexId *>(outputIdentifiers_);
  SimplexId *seg = static_cast<SimplexId *>(outputSegmentation_);

  Timer t;

  std::fill(dist, dist + vertexNumber_, std::numeric_limits<dataType>::max());
  std::fill(origin, origin + vertexNumber_, -1);

  // get the sources
  std::set<SimplexId> isSource;
  for(SimplexId k = 0; k < sourceNumber_; ++k)
    isSource.insert(identifiers[k]);
  std::vector<SimplexId> sources(isSource.begin(), isSource.end());
  isSource.clear();

  // prepare output
  std::vector<std::vector<dataType>> scalars(sources.size());

  // @PETER This doesn't seem to very work efficient, there's multilple source
  // shortest paths algorithms.
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < sources.size(); ++i) {
    int ret = Dijkstra::shortestPath<dataType>(
      sources[i], *triangulation_, scalars[i]);
    if(ret != 0) {
      this->printErr(
        "Algorithm not successful (error code:  " + std::to_string(ret) + ").");
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId k = 0; k < vertexNumber_; ++k) {
    for(size_t i = 0; i < sources.size(); ++i) {
      if(i == 0 or dist[k] > scalars[i][k]) {
        dist[k] = scalars[i][k];
        origin[k] = sources[i];
        seg[k] = i;
      }
    }
  }

  this->printMsg(
    "Complete", 1.0, globalTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}
