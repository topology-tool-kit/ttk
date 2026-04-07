/// \ingroup base
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.

#pragma once

#include <array>
#include <vector>

#include <boost/dynamic_bitset.hpp>

namespace ttk::rpd {
  using id_t = int;
  using value_t = double;
  constexpr value_t inf = std::numeric_limits<value_t>::infinity();

  using PointCloud = std::vector<std::vector<value_t>>;

  using Simplex = std::vector<id_t>;
  using FiltratedSimplex = std::pair<Simplex, value_t>;
  using PersistencePair = std::pair<FiltratedSimplex, FiltratedSimplex>;
  using Diagram = std::vector<PersistencePair>;
  using MultidimensionalDiagram = std::vector<Diagram>;

  using Edge = std::pair<id_t, id_t>;
  using EdgeSet = std::vector<Edge>;
  using EdgeSets3 = std::array<EdgeSet, 3>;
  using EdgeSets4 = std::array<EdgeSet, 4>;
  enum CRIT { DEATH0, BIRTH1, DEATH1, CASC1 };
  using Cascade = EdgeSet;

  using Facet = std::array<id_t, 3>;
  using Generator1 = std::pair<std::vector<Edge>, std::pair<value_t, value_t>>;
  using Generator2 = std::pair<std::vector<Facet>, std::pair<value_t, value_t>>;

  struct FiltratedEdge {
    Edge e;
    value_t d;
  };
  inline FiltratedEdge max(const FiltratedEdge &a, const FiltratedEdge &b) {
    if(a.d > b.d)
      return a;
    return b;
  }

  struct FiltratedQuadEdge {
    Edge e;
    int f1;
    int f2;
    value_t d;
  };

  struct FiltratedTriangle {
    std::tuple<id_t, id_t, id_t> t;
    value_t d;
  };

  inline bool operator<(const FiltratedEdge &e1, const FiltratedEdge &e2) {
    if(e1.d == e2.d)
      return e1.e < e2.e;
    return e1.d < e2.d;
  }

  inline bool operator<(const FiltratedTriangle &f1,
                        const FiltratedTriangle &f2) {
    if(f1.d == f2.d)
      return f1.t < f2.t;
    return f1.d < f2.d;
  }

  class UnionFind {
  private:
    std::vector<int> parent_;
    std::vector<unsigned char> rank_;

  public:
    explicit UnionFind(unsigned n);
    int find(int x);
    void merge(int x, int y);
    int mergeRet(int x, int y);
    [[nodiscard]] bool isRoot(int x) const;
  };

  class BoundaryContainer {
  public:
    BoundaryContainer(std::vector<id_t> &simplices, unsigned size)
      : ids_(simplices) {
      mask_.resize(size, false);
      for(id_t const &id : ids_)
        mask_[id] = true;
    }
    void exclusiveAddBoundary(std::vector<id_t> const &boundary) {
      for(id_t const &id : boundary) {
        if(!mask_[id])
          ids_.emplace_back(id);
        else
          ids_.erase(std::find(ids_.begin(), ids_.end(), id));
        mask_.flip(id);
      }
    }

  private:
    std::vector<id_t> &ids_;
    boost::dynamic_bitset<> mask_;
  };

} // namespace ttk::rpd