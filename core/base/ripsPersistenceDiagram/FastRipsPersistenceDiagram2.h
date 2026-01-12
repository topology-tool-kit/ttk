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
/// MattÃ©o ClÃ©mot, Julie Digne, Julien Tierny, \n
/// IEEE Transactions on Visualization and Computer Graphics.
/// Accepted, to be presented at IEEE VIS 2026.

#pragma once

#ifdef TTK_ENABLE_CGAL

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Triangulation_face_base_with_id_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

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
    explicit FastRipsPersistenceDiagram2(float *data, int n);

    template <typename T>
    void compute0Persistence(T &ph0, bool parallelSort = false);

    template <typename T>
    void computeDelaunayRips0And1Persistence(T &ph, bool parallelSort = false);

    template <typename T>
    void computeRips0And1Persistence(T &ph,
                                     bool parallelSort = false,
                                     bool parallelMML = false);

    void exportRips1Generators(std::vector<Generator1> &generators);

  private:
    Timer tm_{};
    const unsigned n_;
    unsigned nFaces_;
    std::vector<std::pair<Point, unsigned>> points_;
    Delaunay delaunay_;

    std::vector<FiltratedQuadEdge> urquhart_;
    std::vector<FiltratedQuadEdge> rng_;
    std::vector<FiltratedEdge> deathPoly_;
    std::vector<double> birthPoly_;

    // common to Delaunay-Rips and Rips
    void computeDelaunay();
    void computeUrquhart(UnionFind &UF,
                         std::vector<FiltratedEdge> &maxDelaunay,
                         bool parallelSort);
    void compute1PH(std::vector<FiltratedQuadEdge> const &critical,
                    UnionFind &UF,
                    MultidimensionalDiagram &ph);

    // specific to Rips
    [[nodiscard]] static bool isLensEmpty(Point const &p1,
                                          Point const &p2,
                                          Tree const &tree,
                                          double const &d);
    [[nodiscard]] static bool
      isRightSemiLensEmpty(Point const &p1, Point const &p2, Tree const &tree);
    void reindexPolygons(UnionFind const &UF,
                         std::vector<FiltratedEdge> const &maxDelaunay,
                         std::vector<int> &indexPolys);
    void computePolygonRipsDeath(bool parallel,
                                 UnionFind &UF,
                                 std::vector<int> const &indexPolys);
    void pComputePolygonRipsDeath(UnionFind &UF,
                                  std::vector<int> const &indexPolys);
    void executePolygonPairCells(bool parallel,
                                 UnionFind &UF,
                                 std::vector<int> const &indexPolys,
                                 EdgeSets4 &ph) const;

    void static add0Pair(FiltratedQuadEdge const &e, Diagram &ph) {
      ph.emplace_back(FiltratedSimplex{{-1}, 0},
                      FiltratedSimplex{{e.e.first, e.e.second}, e.d});
    }

    void static add0Pair(FiltratedQuadEdge const &e, EdgeSet &ph) {
      ph.emplace_back(e.e);
    }
  };

} // namespace ttk::rpd

#endif
