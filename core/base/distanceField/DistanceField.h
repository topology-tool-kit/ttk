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

#ifndef _DISTANCEFIELD_H
#define _DISTANCEFIELD_H

// base code includes
#include <Dijkstra.h>
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

// std includes
#include <limits>
#include <set>

namespace ttk {

  class DistanceField : public Debug {

  public:
    DistanceField();
    ~DistanceField();

    template <typename dataType>
    dataType getDistance(const SimplexId a, const SimplexId b) const;

    template <typename dataType>
    int execute() const;

    inline int setVertexNumber(SimplexId vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }

    inline int setSourceNumber(SimplexId sourceNumber) {
      sourceNumber_ = sourceNumber;
      return 0;
    }

    inline int setupTriangulation(Triangulation *triangulation) {
      triangulation_ = triangulation;
      if(triangulation_) {
        triangulation_->preprocessVertexNeighbors();
      }
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

    inline int setOutputIdentifiers(void *data) {
      outputIdentifiers_ = data;
      return 0;
    }

    inline int setOutputSegmentation(void *data) {
      outputSegmentation_ = data;
      return 0;
    }

  protected:
    SimplexId vertexNumber_;
    SimplexId sourceNumber_;
    Triangulation *triangulation_;
    void *vertexIdentifierScalarFieldPointer_;
    void *outputScalarFieldPointer_;
    void *outputIdentifiers_;
    void *outputSegmentation_;
  };
} // namespace ttk

template <typename dataType>
int ttk::DistanceField::execute() const {
  SimplexId *identifiers
    = static_cast<SimplexId *>(vertexIdentifierScalarFieldPointer_);
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
  std::vector<SimplexId> sources;
  for(auto s : isSource)
    sources.push_back(s);
  isSource.clear();

  // prepare output
  std::vector<std::vector<dataType>> scalars(sources.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)sources.size(); ++i) {
    Dijkstra::shortestPath<dataType>(sources[i], *triangulation_, scalars[i]);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId k = 0; k < vertexNumber_; ++k) {
    for(SimplexId i = 0; i < (SimplexId)sources.size(); ++i) {
      if(i == 0 or dist[k] > scalars[i][k]) {
        dist[k] = scalars[i][k];
        origin[k] = sources[i];
        seg[k] = i;
      }
    }
  }

  {
    std::stringstream msg;
    msg << "[DistanceField] Data-set (" << vertexNumber_
        << " points) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // DISTANCEFIELD_H
