/// \ingroup base
/// \class ttk::rpd::FastRipsPersistenceDiagram2
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.
///
/// \brief TTK base class that computes the persistence diagram of a Rips
/// complex of a planar point cloud using a fast, dedicated algorithm.
///
/// This module defines the %FastRipsPersistenceDiagram2 class that takes a
/// planar point cloud and computes the persistence diagram of its Rips complex
/// using a geometric-only algorithm based on the computation of the relative
/// neighborhood graph (RNG) and minmax length (MML) triangulations.
///
/// \b Related \b publication \n
/// "Topological Autoencoders++: Fast and Accurate Cycle-Aware Dimensionality
/// Reduction" \n
/// Mattéo Clémot, Julie Digne, Julien Tierny, \n
/// arXiv preprint, 2025.

#pragma once

#ifdef TTK_ENABLE_CGAL

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_id_2.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>

#include <PairCells.h>
#include <RipsPersistenceDiagramUtils.h>

namespace ttk::rpd {

  class FastRipsPersistenceDiagram2 : virtual public Debug {

    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Vb = CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K>;
    using Fb = CGAL::Triangulation_face_base_with_id_2<K>;
    using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using Delaunay = CGAL::Delaunay_triangulation_2<K, Tds>;
    using Point = Delaunay::Point_2;
    using Traits = CGAL::Search_traits_2<K>;
    using Tree = CGAL::Kd_tree<Traits>;
    using Fuzzy_sphere = CGAL::Fuzzy_sphere<Traits>;

  public:
    explicit FastRipsPersistenceDiagram2(const PointCloud &points);
    explicit FastRipsPersistenceDiagram2(float* data, int n);

    template <typename T>
    void compute0Persistence(T &ph0, bool parallelSort = false);

    template <typename T>
    void computeDelaunayRips0And1Persistence(T &ph, bool parallelSort = false);

    template <typename T>
    void computeRips0And1Persistence(T &ph, bool parallelSort = false, bool parallelMML = false);

    void exportRips1Generators(std::vector<Generator> &generators);

  private:
    Timer tm_ {};
    const unsigned n_;
    unsigned nFaces_;
    std::vector<std::pair<Point,unsigned>> points_;
    Delaunay delaunay_;

    std::vector<FiltratedQuadEdge> urquhart_;
    std::vector<FiltratedQuadEdge> rng_;
    std::vector<FiltratedEdge> deathPoly_;
    std::vector<double> birthPoly_;

    //common to Delaunay-Rips and Rips
    void computeDelaunay();
    void computeUrquhart(UnionFind& UF, std::vector<FiltratedEdge>& maxDelaunay, bool parallelSort);
    void compute1PH(std::vector<FiltratedQuadEdge> const& critical, UnionFind &UF, MultidimensionalDiagram &ph);

    //specific to Rips
    [[nodiscard]] static bool isLensEmpty(Point const& p1, Point const& p2, Tree const& tree, double const& d);
    [[nodiscard]] static bool isRightSemiLensEmpty(Point const& p1, Point const& p2, Tree const& tree);
    void reindexPolygons(UnionFind const& UF, std::vector<FiltratedEdge> const& maxDelaunay, std::vector<int>& indexPolys);
    void computePolygonRipsDeath(bool parallel, UnionFind &UF, std::vector<int> const& indexPolys);
    void pComputePolygonRipsDeath(UnionFind &UF, std::vector<int> const& indexPolys);
    void executePolygonPairCells(bool parallel, UnionFind &UF, std::vector<int> const& indexPolys, EdgeSets4 &ph) const;

    void static add0Pair(FiltratedQuadEdge const& e, Diagram &ph) {
      ph.emplace_back(FiltratedSimplex{{-1}, 0}, FiltratedSimplex{{e.e.first, e.e.second},e.d});
    }

    void static add0Pair(FiltratedQuadEdge const& e, EdgeSet &ph) {
      ph.emplace_back(e.e);
    }

  };

  template <typename T>
  void FastRipsPersistenceDiagram2::compute0Persistence(T &ph0, bool parallelSort) {
    ph0 = T(0);

    //keep only the Urquhart edges, maintain polygon structure
    UnionFind UF(nFaces_);
    std::vector<FiltratedEdge> max_delaunay(nFaces_, FiltratedEdge{{-1,-1},0.});
    computeUrquhart(UF, max_delaunay, parallelSort);

    // compute EMST with Kruskal algorithm
    UnionFind UF_p(n_);
    for (FiltratedQuadEdge const& e : urquhart_) {
      if(UF_p.find(e.e.first) != UF_p.find(e.e.second)) { //we know e is a EMST edge
        UF_p.merge(e.e.first, e.e.second);
        add0Pair(e, ph0);
      }
    }
    if constexpr(std::is_same_v<T, Diagram>)
      ph0.emplace_back(FiltratedSimplex{{-1}, 0.}, FiltratedSimplex{{-1}, inf}); //infinite pair

    printMsg("MST computed", 0., tm_.getElapsedTime());
  }

}

#endif