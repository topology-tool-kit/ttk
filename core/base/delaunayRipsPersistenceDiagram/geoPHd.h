#pragma once

#include "geoPHUtils.h"

#ifdef TTK_ENABLE_CGAL

#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Triangulation_vertex.h>
#include <CGAL/Triangulation_full_cell.h>

namespace ttk::gph {

  template <unsigned D>
  using DSimplex = std::array<id_t, D+1>;

  template <unsigned D>
  using ConnectivityHashMap = HashMap<DSimplex<D>, std::vector<std::pair<int,int>>>;

  template <unsigned D>
  struct FiltratedDSimplex {
    DSimplex<D> s;
    double d;
    double a;
  };

  template <unsigned D>
  bool operator< (FiltratedDSimplex<D> const& s1, FiltratedDSimplex<D> const& s2) {
    return std::make_pair(s1.d, s1.a) < std::make_pair(s2.d, s2.a);
  }

  template <unsigned D>
  FiltratedDSimplex<D> max(FiltratedDSimplex<D> s1, FiltratedDSimplex<D> s2) {
    return s1 < s2 ? s2 : s1;
  }

  template <unsigned D>
  using DSimplicialComplex = std::vector<FiltratedDSimplex<D>>;

  template <unsigned DIM>
  class DRPersistenceD {
    using DimTag = CGAL::Dimension_tag<DIM>;
    using K = CGAL::Epick_d<DimTag>;
    using Vb = CGAL::Triangulation_vertex<K, id_t>;
    using FCb = CGAL::Triangulation_full_cell<K, id_t>;
    using Tds = CGAL::Triangulation_data_structure<DimTag, Vb, FCb>;
    using Delaunay = CGAL::Delaunay_triangulation<K, Tds>;
    using Point = typename Delaunay::Point;
    using Facet = typename Delaunay::Facet;
    using VertexHandle = typename Delaunay::Vertex_handle;
    using CellHandle = typename Delaunay::Full_cell_handle;

    struct FiltratedQuadFacet {
      DSimplex<DIM-1> s;
      double d;
      double a;
      int c1;
      int c2;
    };

  public:
    explicit DRPersistenceD(PointCloud<DIM> &points);
    void run(MultidimensionalDiagram &ph);

  private:
    const unsigned N_p;
    unsigned N_c {};
    Delaunay del_;
    PointCloud<DIM> &points_;

    [[nodiscard]] double squaredDistance(const unsigned i1, const unsigned i2) const {
      return std::inner_product(points_[i1].begin(),points_[i1].end(), points_[i2].begin(), 0.,
      std::plus(), [](double u, double v) { return (u - v) * (u - v); });
    }

    template <unsigned D>
    std::pair<double,double> squaredPerturbedDiameter (DSimplex<D> const& s) const {
      double diam = 0.;
      double a = 0.;
      for (unsigned p1 = 0; p1 < D+1; ++p1) {
        for (unsigned p2 = p1+1; p2 < D+1; ++p2) {
          const double dist = squaredDistance(s[p1], s[p2]);
          diam = std::max(diam, dist);
          a += dist;
        }
      }
      return {diam, a};
    }

    template <unsigned D>
    std::array<double, D*(D+1)/2> getLengths (DSimplex<D> const& s) const {
      std::array<double, D*(D+1)/2> lengths;
      int i = 0;
      for (unsigned p1 = 0; p1 < D+1; ++p1) {
        for (unsigned p2 = p1+1; p2 < D+1; ++p2)
          lengths[i++] = squaredDistance(s[p1], s[p2]);
      }
      return lengths;
    }

    template <unsigned D>
    bool checkLinkUrquhart(DSimplex<D> const& s, std::array<double, D*(D+1)/2> const& lengths, std::pair<double,double> diam, id_t linkId) {
      std::array<double, D+1> linkLengths;
      for (unsigned i = 0; i < D+1; ++i)
        linkLengths[i] = squaredDistance(s[i], linkId);

      for (unsigned i = 0; i<D+1; ++i) {
        double d = 0.;
        double a = 0.;
        int k=0;
        for (unsigned p1 = 0; p1 < D+1; ++p1) {
          for (unsigned p2 = p1+1; p2 < D+1; ++p2) {
            double dist;
            if (p1 == i)
              dist = linkLengths[p2];
            else if (p2 == i)
              dist = linkLengths[p1];
            else
              dist = lengths[k];

            d = std::max(d, dist);
            a += dist;
            ++k;
          }
        }
        if (std::make_pair(d, a) > diam)
          return false;
      }

      return true;
    }

    void computeDelaunay();
    void computeDPH(Diagram &ph,
                    DSimplicialComplex<DIM-1> &MSA);

    template <unsigned D>
    void computeNextPH(Diagram &ph,
                       DSimplicialComplex<D> const& MSA,
                       DSimplicialComplex<D-1> &nextMSA);

    template <unsigned D>
    void recurse(MultidimensionalDiagram &ph,
                 DSimplicialComplex<D> const& MSA);
  };

  template <unsigned DIM>
  DRPersistenceD<DIM>::DRPersistenceD(PointCloud<DIM> &points) : N_p(points.size()), del_(DIM), points_(points) {}

  template <unsigned DIM>
  void DRPersistenceD<DIM>::run(MultidimensionalDiagram &ph) {
    ph = MultidimensionalDiagram(DIM);

    computeDelaunay();

    DSimplicialComplex<DIM-1> MSA1;
    computeDPH(ph[DIM-1], MSA1);

    recurse(ph, MSA1);
  }

  template <unsigned DIM>
  template <unsigned D>
  void DRPersistenceD<DIM>::recurse(MultidimensionalDiagram &ph,
                                                DSimplicialComplex<D> const& MSA) {
    if constexpr (D >= 2) {
      DSimplicialComplex<D-1> nextMSA;
      computeNextPH(ph[D-1], MSA, nextMSA);
      recurse(ph, nextMSA);
    }
    else if constexpr (D == 1) { //this is the 0-dimensional homology
      for (FiltratedDSimplex<1> const& e : MSA)
        ph[0].emplace_back(FiltratedSimplex{{}, 0.},
                           FiltratedSimplex{{}, sqrt(e.d)});
    }
  }

  template <unsigned DIM>
  void DRPersistenceD<DIM>::computeDelaunay() {
    del_.insert(points_.begin(), points_.end());

    //index all vertices
    int p = 0;
    for (auto p_it = del_.finite_vertices_begin(); p_it != del_.finite_vertices_end(); ++p_it) {
      points_[p] = p_it->point();
      p_it->data() = p++;
    }

    //index all cells (including infinite ones)
    int k = 0;
    for (auto c_it = del_.full_cells_begin() ; c_it != del_.full_cells_end(); ++c_it)
      c_it->data() = k++;
    N_c = k;
  }

  template<unsigned DIM>
  void DRPersistenceD<DIM>::computeDPH(Diagram &ph,
                                                   DSimplicialComplex<DIM-1> &MSA) {
    DSimplicialComplex<DIM-1> maxDelaunay (N_c);
    UnionFind UF(N_c);
    std::vector<FiltratedQuadFacet> hyperUrquhart;

    /* computing urquhart hypergraph (codimension 1) */
    for (auto f_it = del_.facets_begin(); f_it != del_.facets_end(); ++f_it) {

      const Facet f = *f_it;
      const CellHandle c = f.full_cell();
      const CellHandle c_mirror = c->neighbor(f.index_of_covertex());

      if (del_.is_infinite(f))
        UF.merge(c->data(), c_mirror->data());

      else { //we need to determine whether f is "Urquhart"
        const auto linkPoint1 = c->vertex(f.index_of_covertex());
        const auto linkPoint2 = c_mirror->vertex(c->mirror_index(f.index_of_covertex()));

        DSimplex<DIM-1> facet;
        for (unsigned i = 0; i<DIM; ++i)
          facet[i] = c->vertex((f.index_of_covertex() + i + 1) % (DIM+1))->data();
        std::sort(facet.begin(), facet.end());
        const auto diam = squaredPerturbedDiameter<DIM-1>(facet);

        bool is_urquhart = true;
        for (auto linkPoint : {linkPoint1, linkPoint2}) {
          if (!del_.is_infinite(linkPoint)) {
            const id_t k = linkPoint->data();
            bool largest = true;
            for (unsigned i = 0; i<DIM; ++i) {
              DSimplex<DIM-1> neighbor = facet;
              neighbor[i] = k;
              if (squaredPerturbedDiameter<DIM-1>(neighbor) > diam) {
                largest = false;
                break;
              }
            }
            if (largest) {
              is_urquhart = false;
              break;
            }
          }
        }

        if (is_urquhart)
          hyperUrquhart.push_back({facet, diam.first, diam.second, c->data(), c_mirror->data()});
        else {
          const int poly1 = UF.find(c->data());
          const int poly2 = UF.find(c_mirror->data());
          maxDelaunay[UF.mergeRet(poly1, poly2)] = max(FiltratedDSimplex<DIM-1>{facet, diam.first, diam.second},
                                                           max(maxDelaunay[poly1], maxDelaunay[poly2]));
        }

        //detect infinite polyhedrons (going beyond convex hull)
        if (del_.is_infinite(c))
          maxDelaunay[UF.find(c->data())].d = inf;
        else if (del_.is_infinite(c_mirror))
          maxDelaunay[UF.find(c_mirror->data())].d = inf;
      }
    }

    std::sort(hyperUrquhart.begin(),
              hyperUrquhart.end(),
              [](const FiltratedQuadFacet &f1, const FiltratedQuadFacet &f2) {
      if (f1.d == f2.d)
        return f1.a > f2.a;
      else
        return f1.d > f2.d;
    });

    /* reverse-delete algorithm to determine MSA */
    std::vector<int> latest(maxDelaunay.size());
    std::iota(latest.begin(), latest.end(), 0);

    for (FiltratedQuadFacet const& f : hyperUrquhart) { //sorted by decreasing order
      const int v1 = UF.find(f.c1);
      const int v2 = UF.find(f.c2);
      if (v1 != v2) { // two distinct codimension-1 cavities: merge them by deleting the facet
        UF.merge(v1, v2);

        const int latest1 = latest[v1];
        const int latest2 = latest[v2];
        const FiltratedDSimplex<DIM-1>& death1 = maxDelaunay[latest1];
        const FiltratedDSimplex<DIM-1>& death2 = maxDelaunay[latest2];

        if (death1.d < death2.d) {
          if (f.d < death1.d)
            ph.emplace_back(FiltratedSimplex{{}, sqrt(f.d)},
                            FiltratedSimplex{{}, sqrt(death1.d)});
          latest[UF.find(v1)] = latest2;
        }
        else if (death2.d < death1.d) {
          if (f.d < death2.d)
            ph.emplace_back(FiltratedSimplex{{}, sqrt(f.d)},
                            FiltratedSimplex{{}, sqrt(death2.d)});
          latest[UF.find(v1)] = latest1;
        }
      }
      else // this is a facet from the minimal spanning acycle
        MSA.push_back({f.s, f.d, f.a});
    }

  }

  template <unsigned DIM>
  template <unsigned D>
  void DRPersistenceD<DIM>::computeNextPH(Diagram &ph,
                                                      DSimplicialComplex<D> const& MSA,
                                                      DSimplicialComplex<D-1> &nextMSA) {
    /* Connectivity */

    ConnectivityHashMap<D-1> msa_connectivity;
    msa_connectivity.reserve(MSA.size());
    for (unsigned i = 0; i<MSA.size(); ++i) {
      DSimplex<D-1> face;
      for (unsigned k=0; k<D+1; ++k) {
        for (unsigned j=0; j<D; ++j)
          face[j] = MSA[i].s[j + (j>=k)];
        msa_connectivity[face].reserve(4); //todo adjust guess
        msa_connectivity[face].emplace_back(i, MSA[i].s[k]);
      }
    }

    /* Urquhart-ness and Urquhart-polytopes */

    std::vector<FiltratedDSimplex<D-1>> critical;
    UnionFind UF_msa (MSA.size());
    std::vector<FiltratedDSimplex<D-1>> maxDelaunay (MSA.size());
    for (auto const& [s, neighbors] : msa_connectivity) {

      // first determine whether s is Urquhart
      bool is_urquhart = true;
      const auto lengths = getLengths<D-1>(s);
      std::pair<double,double> diam;
      for (double const& l : lengths)
        diam = {std::max(diam.first, l), diam.second + l};

      for (auto [coface_id, linkPoint_id] : neighbors) {
        if (checkLinkUrquhart<D-1>(s, lengths, diam, linkPoint_id)) {
          is_urquhart = false;
          break;
        }
      }

      // now maintain the polytope structure
      if (is_urquhart)
        critical.push_back({s, diam.first, diam.second});
      else {
        if (neighbors.size() == 1) //set infinite polytope
          maxDelaunay[UF_msa.find(neighbors[0].first)].d = inf;
        else if (neighbors.size() == 2) {
          const int poly1 = UF_msa.find(neighbors[0].first);
          const int poly2 = UF_msa.find(neighbors[1].first);
          maxDelaunay[UF_msa.mergeRet(poly1, poly2)] = max(FiltratedDSimplex<D-1>{s, diam.first, diam.second},
                                                                max(maxDelaunay[poly1], maxDelaunay[poly2]));
        }
        else {
          critical.push_back({s, diam.first, diam.second});
          for (auto const& f_id : neighbors) {
            const int poly = UF_msa.find(f_id.first);
            maxDelaunay[poly] = max(FiltratedDSimplex<D-1>{s, diam.first, diam.second}, maxDelaunay[poly]);
          }
        }
      }
    }

    /* Graph critical -- polytope */

    std::vector<id_t> polytopes;
    for (unsigned x=0; x<MSA.size(); ++x) {
      if (UF_msa.isRoot(x) && maxDelaunay[x].d < inf)
        polytopes.emplace_back(x);
    }

    std::sort(polytopes.begin(), polytopes.end(), [&](const int x1, const int x2) {
      return maxDelaunay[x1] < maxDelaunay[x2];
    });

    std::vector<int> criticalIndices(critical.size());
    std::iota(criticalIndices.begin(), criticalIndices.end(), 0);
    std::sort(criticalIndices.begin(), criticalIndices.end(), [&](const int x1, const int x2) {
          return critical[x1].d < critical[x2].d;
        });
    std::vector<int> criticalOrder(critical.size());
    for (unsigned i=0; i<criticalIndices.size(); ++i)
      criticalOrder[criticalIndices[i]] = i;

    std::vector<std::vector<int>> poly_to_crit(MSA.size());
    for (const int poly : polytopes)
      poly_to_crit[poly].reserve(D+1); //todo adjust guess
    for (unsigned i=0; i<critical.size(); ++i) {
      const FiltratedDSimplex<D-1> c = critical[i];
      for (const auto& [poly,_] : msa_connectivity[c.s]) {
        if (maxDelaunay[UF_msa.find(poly)].d < inf) {
          auto &neighbors = poly_to_crit[UF_msa.find(poly)];
          auto it = std::find(neighbors.begin(), neighbors.end(), criticalOrder[i]);
          if (it == neighbors.end())
            neighbors.push_back(criticalOrder[i]);
          else
            neighbors.erase(it);
        }
      }
    }

    /* PairCells */

    std::vector<int> partner(critical.size(), -1);
    for (const int poly : polytopes) {
      std::set boundary(poly_to_crit[poly].begin(), poly_to_crit[poly].end());
      while (true) {
        const int youngest_id = *boundary.rbegin();
        if (partner[youngest_id] == -1) {
          partner[youngest_id] = poly;
          const FiltratedDSimplex<D-1> &c = critical[criticalIndices[youngest_id]];
          const FiltratedDSimplex<D-1> &death = maxDelaunay[poly];
          if (c.d < death.d)
            ph.emplace_back(FiltratedSimplex{{}, sqrt(c.d)},
                            FiltratedSimplex{{}, sqrt(death.d)});
          break;
        }
        else {
          for (const int crit : poly_to_crit[partner[youngest_id]]) {
            const auto it = boundary.find(crit);
            if (it == boundary.end())
              boundary.insert(crit);
            else
              boundary.erase(it);
          }
        }
      }
      poly_to_crit[poly].assign(boundary.begin(), boundary.end());
    }

    /* Next MSA */
    for (unsigned i=0; i<critical.size(); ++i) {
      if (partner[i] == -1) // unassigned -> go in next MSA
        nextMSA.emplace_back(critical[criticalIndices[i]]);
    }
  }

  template <unsigned DIM>
  void runDelaunayRipsPersistenceDiagram(rpd::PointCloud const& points, MultidimensionalDiagram &diagram) {
    PointCloud<DIM> p(points.size());
    for (unsigned i = 0; i < points.size(); ++i) {
      for (unsigned d = 0; d < DIM; ++d)
        p[i][d] = points[i][d];
    }
    DRPersistenceD<DIM> drpd(p);
    drpd.run(diagram);
  }

  template <unsigned DIM>
  void tryDimension(rpd::PointCloud const& points, MultidimensionalDiagram &diagram) {
    if constexpr (DIM <= TTK_DELAUNAY_MAXIMUM_DIMENSION) {
      if (points[0].size() == DIM)
        runDelaunayRipsPersistenceDiagram<DIM>(points, diagram);
      else
        tryDimension<DIM+1>(points, diagram);
    }
  }

  inline void tryDimensions(rpd::PointCloud const& points, MultidimensionalDiagram &diagram) {
    tryDimension<4>(points, diagram);
  }

}

#endif