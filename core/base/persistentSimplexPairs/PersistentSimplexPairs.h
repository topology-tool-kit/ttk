/// \ingroup base
/// \class ttk::PersistentSimplexPairs
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date November 2021.
///
/// \brief Textbook algorithm to find persistence pairs.
///
/// This algorithm is described in "Algorithm and Theory of
/// Computation Handbook (Second Edition) - Special Topics and
/// Techniques" by Atallah and Blanton on page 97.

#pragma once

#include <AbstractTriangulation.h>
#include <Debug.h>
#include <VisitedMask.h>

#include <algorithm>
#include <string>
#include <vector>

namespace ttk {

  class PersistentSimplexPairs : virtual public Debug {
  public:
    PersistentSimplexPairs();

    struct PersistencePair {
      /** first (lower) vertex id */
      SimplexId birth;
      /** second (higher) vertex id */
      SimplexId death;
      /** pair type (min-saddle: 0, saddle-saddle: 1, saddle-max: 2) */
      SimplexId type;

      PersistencePair(SimplexId b, SimplexId d, SimplexId t)
        : birth{b}, death{d}, type{t} {
      }
    };

    /**
     * @brief Preprocess all the required connectivity requests on the
     * triangulation.
     */
    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      if(data != nullptr) {
        const auto dim = data->getDimensionality();
        data->preconditionEdges();
        if(dim == 2) {
          data->preconditionCellEdges();
        } else if(dim == 3) {
          data->preconditionTriangles();
          data->preconditionTriangleEdges();
          data->preconditionCellTriangles();
        }
        this->nVerts_ = data->getNumberOfVertices();
        this->nEdges_ = data->getNumberOfEdges();
        this->nTri_ = dim > 1 ? data->getNumberOfTriangles() : 0;
        this->nTetra_ = dim > 2 ? data->getNumberOfCells() : 0;
      }
    }

    /**
     * @brief Compute the persistence pairs from the triangulation
     * simplicial complex
     */
    template <typename triangulationType>
    int computePersistencePairs(std::vector<PersistencePair> &pairs,
                                const SimplexId *const orderField,
                                const triangulationType &triangulation) const;

  private:
    struct Simplex {
      SimplexId dim_{-1}; // dimension
      SimplexId id_{-1}; // id in triangulation (overlap between dimensions)
      SimplexId cellId_{-1}; // cell id (unique)
      // face (triangulation) indices
      std::array<SimplexId, 4> faceIds_{-1, -1, -1, -1};
      // order on vertices, sorted in descending order
      std::array<SimplexId, 4> vertsOrder_{-1, -1, -1, -1};

      friend bool operator<(const Simplex &lhs, const Simplex &rhs) {
        return lhs.vertsOrder_ < rhs.vertsOrder_;
      }

      void fillVert(const SimplexId v, const SimplexId *const offset) {
        this->dim_ = 0;
        this->id_ = v;
        this->cellId_ = v;
        this->vertsOrder_[0] = offset[v];
      }

      template <typename triangulationType>
      void fillEdge(const SimplexId e,
                    const SimplexId c,
                    const SimplexId *const offset,
                    const triangulationType &triangulation) {
        this->dim_ = 1;
        this->id_ = e;
        this->cellId_ = c;
        triangulation.getEdgeVertex(e, 0, this->faceIds_[0]);
        triangulation.getEdgeVertex(e, 1, this->faceIds_[1]);
        this->vertsOrder_[0] = offset[this->faceIds_[0]];
        this->vertsOrder_[1] = offset[this->faceIds_[1]];
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }

      template <typename triangulationType>
      void fillTriangle(const SimplexId t,
                        const SimplexId c,
                        const SimplexId *const offset,
                        const triangulationType &triangulation) {
        this->dim_ = 2;
        this->id_ = t;
        this->cellId_ = c;
        triangulation.getTriangleEdge(t, 0, this->faceIds_[0]);
        triangulation.getTriangleEdge(t, 1, this->faceIds_[1]);
        triangulation.getTriangleEdge(t, 2, this->faceIds_[2]);
        triangulation.getTriangleVertex(t, 0, this->vertsOrder_[0]);
        triangulation.getTriangleVertex(t, 1, this->vertsOrder_[1]);
        triangulation.getTriangleVertex(t, 2, this->vertsOrder_[2]);
        this->vertsOrder_[0] = offset[this->vertsOrder_[0]];
        this->vertsOrder_[1] = offset[this->vertsOrder_[1]];
        this->vertsOrder_[2] = offset[this->vertsOrder_[2]];
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }

      template <typename triangulationType>
      void fillTetra(const SimplexId T,
                     const SimplexId c,
                     const SimplexId *const offset,
                     const triangulationType &triangulation) {
        this->dim_ = 3;
        this->id_ = T;
        this->cellId_ = c;
        triangulation.getCellTriangle(T, 0, this->faceIds_[0]);
        triangulation.getCellTriangle(T, 1, this->faceIds_[1]);
        triangulation.getCellTriangle(T, 2, this->faceIds_[2]);
        triangulation.getCellTriangle(T, 3, this->faceIds_[3]);
        triangulation.getCellVertex(T, 0, this->vertsOrder_[0]);
        triangulation.getCellVertex(T, 1, this->vertsOrder_[1]);
        triangulation.getCellVertex(T, 2, this->vertsOrder_[2]);
        triangulation.getCellVertex(T, 3, this->vertsOrder_[3]);
        this->vertsOrder_[0] = offset[this->vertsOrder_[0]];
        this->vertsOrder_[1] = offset[this->vertsOrder_[1]];
        this->vertsOrder_[2] = offset[this->vertsOrder_[2]];
        this->vertsOrder_[3] = offset[this->vertsOrder_[3]];
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }
    };

    template <typename triangulationType>
    std::vector<Simplex>
      computeFiltrationOrder(const SimplexId *const offset,
                             const triangulationType &triangulation) const;

    inline void addCellBoundary(const Simplex &c, VisitedMask &boundary) const {
      for(SimplexId i = 0; i < c.dim_ + 1; ++i) {
        const auto f{c.faceIds_[i]};
        if(!boundary.isVisited_[f]) {
          boundary.insert(f);
        } else {
          boundary.remove(f);
        }
      }
    }

    inline SimplexId getCellId(const SimplexId cdim,
                               const SimplexId cid) const {
      if(cdim == 0) {
        return cid;
      } else if(cdim == 1) {
        return cid + this->nVerts_;
      } else if(cdim == 2) {
        return cid + this->nVerts_ + this->nEdges_;
      } else if(cdim == 3) {
        return cid + this->nVerts_ + this->nEdges_ + this->nTri_;
      }
      return -1;
    }

    SimplexId eliminateBoundaries(const Simplex &c,
                                  VisitedMask &boundary,
                                  const std::vector<SimplexId> &filtOrder,
                                  const std::vector<Simplex> &partners) const;

    int pairCells(std::vector<PersistencePair> &pairs,
                  std::array<std::vector<bool>, 3> &boundaries,
                  const std::vector<Simplex> &filtration,
                  const std::vector<SimplexId> &filtOrder) const;

    SimplexId nVerts_{0};
    SimplexId nEdges_{0};
    SimplexId nTri_{0};
    SimplexId nTetra_{0};
  };

} // namespace ttk

template <typename triangulationType>
int ttk::PersistentSimplexPairs::computePersistencePairs(
  std::vector<ttk::PersistentSimplexPairs::PersistencePair> &pairs,
  const SimplexId *const orderField,
  const triangulationType &triangulation) const {

  Timer tm{};

  // every simplex in the triangulation, sorted by filtration
  const auto filtration
    = this->computeFiltrationOrder(orderField, triangulation);

  std::array<std::vector<bool>, 3> boundaries{};
  boundaries[0].resize(this->nVerts_, false);
  boundaries[1].resize(this->nEdges_, false);
  boundaries[2].resize(this->nTri_, false);

  // simplex id -> filtration order
  std::vector<SimplexId> filtOrder(filtration.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < filtration.size(); ++i) {
    filtOrder[filtration[i].cellId_] = i;
  }

  this->pairCells(pairs, boundaries, filtration, filtOrder);

  this->printMsg(
    "Computed " + std::to_string(pairs.size()) + " persistence pairs", 1.0,
    tm.getElapsedTime(), 1);

  return 0;
}

template <typename triangulationType>
std::vector<ttk::PersistentSimplexPairs::Simplex>
  ttk::PersistentSimplexPairs::computeFiltrationOrder(
    const SimplexId *const offset,
    const triangulationType &triangulation) const {

  Timer tm{};

  std::vector<Simplex> res(this->nVerts_ + this->nEdges_ + this->nTri_
                           + this->nTetra_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < this->nVerts_; ++i) {
    res[i].fillVert(i, offset);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < this->nEdges_; ++i) {
    const auto o = this->nVerts_ + i;
    res[o].fillEdge(i, o, offset, triangulation);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < this->nTri_; ++i) {
    const auto o = this->nVerts_ + this->nEdges_ + i;
    res[o].fillTriangle(i, o, offset, triangulation);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < this->nTetra_; ++i) {
    const auto o = this->nVerts_ + this->nEdges_ + this->nTri_ + i;
    res[o].fillTetra(i, o, offset, triangulation);
  }

  TTK_PSORT(this->threadNumber_, res.begin(), res.end());

  this->printMsg(
    "Computed filtration order", 1.0, tm.getElapsedTime(), this->threadNumber_);

  return res;
}
