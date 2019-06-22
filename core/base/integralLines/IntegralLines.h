/// \ingroup base
/// \class ttk::IntegralLines
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date March 2016
///
/// \brief TTK processing package for the computation of edge-based integral
/// lines of the gradient of an input scalar field defined on a PL manifold.
///
/// Given a list of sources, the package produces forward or backward integral
/// lines along the edges of the input triangulation.
///
/// \sa ttkIntegralLines.cpp %for a usage example.

#ifndef _DISCRETESTREAMLINE_H
#define _DISCRETESTREAMLINE_H

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

// std includes
#include <unordered_set>

namespace ttk {
  enum Direction { Forward = 0, Backward };

  class IntegralLines : public Debug {

  public:
    IntegralLines();
    ~IntegralLines();

    template <typename dataType>
    inline float getDistance(const SimplexId &a, const SimplexId &b) const {
      float p0[3];
      triangulation_->getVertexPoint(a, p0[0], p0[1], p0[2]);
      float p1[3];
      triangulation_->getVertexPoint(b, p1[0], p1[1], p1[2]);

      return Geometry::distance(p0, p1, 3);
    }

    template <typename dataType>
    inline float getGradient(const SimplexId &a,
                             const SimplexId &b,
                             dataType *scalars) const {
      return fabs(scalars[b] - scalars[a]) / getDistance<dataType>(a, b);
    }

    template <typename dataType, typename idType>
    int execute() const;

    template <typename dataType, typename idType, class Compare>
    int execute(Compare cmp) const;

    inline int setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }

    inline int setSeedNumber(const SimplexId &seedNumber) {
      seedNumber_ = seedNumber;
      return 0;
    }

    inline int setDirection(int direction) {
      direction_ = direction;
      return 0;
    }

    inline int setupTriangulation(Triangulation *triangulation) {
      triangulation_ = triangulation;
      if(triangulation_) {
        triangulation_->preprocessVertexNeighbors();
      }
      return 0;
    }

    inline int setInputScalarField(void *data) {
      inputScalarField_ = data;
      return 0;
    }

    inline int setInputOffsets(void *data) {
      inputOffsets_ = data;
      return 0;
    }

    inline int setVertexIdentifierScalarField(void *data) {
      vertexIdentifierScalarField_ = data;
      return 0;
    }

    inline int
      setOutputTrajectories(std::vector<std::vector<SimplexId>> *trajectories) {
      outputTrajectories_ = trajectories;
      return 0;
    }

  protected:
    SimplexId vertexNumber_;
    SimplexId seedNumber_;
    int direction_;
    Triangulation *triangulation_;
    void *inputScalarField_;
    void *inputOffsets_;
    void *vertexIdentifierScalarField_;
    std::vector<std::vector<SimplexId>> *outputTrajectories_;
  };
} // namespace ttk

template <typename dataType, typename idType>
int ttk::IntegralLines::execute() const {
  idType *offsets = static_cast<idType *>(inputOffsets_);
  SimplexId *identifiers
    = static_cast<SimplexId *>(vertexIdentifierScalarField_);
  dataType *scalars = static_cast<dataType *>(inputScalarField_);
  std::vector<std::vector<SimplexId>> *trajectories = outputTrajectories_;

  Timer t;

  // get the seeds
  std::unordered_set<SimplexId> isSeed;
  for(SimplexId k = 0; k < seedNumber_; ++k)
    isSeed.insert(identifiers[k]);
  std::vector<SimplexId> seeds;
  for(auto k : isSeed)
    seeds.push_back(k);
  isSeed.clear();

  trajectories->resize(seeds.size());
  for(SimplexId i = 0; i < (SimplexId)seeds.size(); ++i) {
    SimplexId v{seeds[i]};
    (*trajectories)[i].push_back(v);

    bool isMax{};
    while(!isMax) {
      SimplexId vnext{-1};
      float fnext = std::numeric_limits<float>::min();
      SimplexId neighborNumber = triangulation_->getVertexNeighborNumber(v);
      bool isLocalMax = true;
      bool isLocalMin = true;
      for(SimplexId k = 0; k < neighborNumber; ++k) {
        SimplexId n;
        triangulation_->getVertexNeighbor(v, k, n);

        if(scalars[n] <= scalars[v])
          isLocalMax = false;
        if(scalars[n] >= scalars[v])
          isLocalMin = false;

        if((direction_ == static_cast<int>(Direction::Forward))
           xor (scalars[n] < scalars[v])) {
          const float f = getGradient<dataType>(v, n, scalars);
          if(f > fnext) {
            vnext = n;
            fnext = f;
          }
        }
      }

      if(vnext == -1 and !isLocalMax and !isLocalMin) {
        idType onext = -1;
        for(SimplexId k = 0; k < neighborNumber; ++k) {
          SimplexId n;
          triangulation_->getVertexNeighbor(v, k, n);

          if(scalars[n] == scalars[v]) {
            const idType o = offsets[n];
            if((direction_ == static_cast<int>(Direction::Forward))
               xor (o < offsets[v])) {
              if(o > onext) {
                vnext = n;
                onext = o;
              }
            }
          }
        }
      }

      if(vnext == -1)
        isMax = true;
      else {
        v = vnext;
        (*trajectories)[i].push_back(v);
      }
    }
  }

  {
    std::stringstream msg;
    msg << "[IntegralLines] Data-set (" << vertexNumber_
        << " points) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType, typename idType, class Compare>
int ttk::IntegralLines::execute(Compare cmp) const {
  idType *offsets = static_cast<idType *>(inputOffsets_);
  SimplexId *identifiers
    = static_cast<SimplexId *>(vertexIdentifierScalarField_);
  dataType *scalars = static_cast<dataType *>(inputScalarField_);
  std::vector<std::vector<SimplexId>> *trajectories = outputTrajectories_;

  Timer t;

  // get the seeds
  std::unordered_set<SimplexId> isSeed;
  for(SimplexId k = 0; k < seedNumber_; ++k)
    isSeed.insert(identifiers[k]);
  std::vector<SimplexId> seeds;
  for(auto k : isSeed)
    seeds.push_back(k);
  isSeed.clear();

  trajectories->resize(seeds.size());
  for(SimplexId i = 0; i < (SimplexId)seeds.size(); ++i) {
    SimplexId v{seeds[i]};
    (*trajectories)[i].push_back(v);

    bool isMax{};
    while(!isMax) {
      SimplexId vnext{-1};
      float fnext = std::numeric_limits<float>::min();
      SimplexId neighborNumber = triangulation_->getVertexNeighborNumber(v);
      bool isLocalMax = true;
      bool isLocalMin = true;
      for(SimplexId k = 0; k < neighborNumber; ++k) {
        SimplexId n;
        triangulation_->getVertexNeighbor(v, k, n);

        if(scalars[n] <= scalars[v])
          isLocalMax = false;
        if(scalars[n] >= scalars[v])
          isLocalMin = false;

        if((direction_ == static_cast<int>(Direction::Forward))
           xor (scalars[n] < scalars[v])) {
          const float f = getGradient<dataType>(v, n, scalars);
          if(f > fnext) {
            vnext = n;
            fnext = f;
          }
        }
      }

      if(vnext == -1 and !isLocalMax and !isLocalMin) {
        SimplexId onext = -1;
        for(SimplexId k = 0; k < neighborNumber; ++k) {
          SimplexId n;
          triangulation_->getVertexNeighbor(v, k, n);

          if(scalars[n] == scalars[v]) {
            const SimplexId o = offsets[n];
            if((direction_ == static_cast<int>(Direction::Forward))
               xor (o < offsets[v])) {
              if(o > onext) {
                vnext = n;
                onext = o;
              }
            }
          }
        }
      }

      if(vnext == -1)
        isMax = true;
      else {
        v = vnext;
        (*trajectories)[i].push_back(v);

        if(cmp(v))
          isMax = true;
      }
    }
  }

  {
    std::stringstream msg;
    msg << "[IntegralLines] Data-set (" << vertexNumber_
        << " points) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // DISCRETESTREAMLINE_H
