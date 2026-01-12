#pragma once

#include <geoPHUtils.h>

#ifdef TTK_ENABLE_CGAL

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation_full_cell.h>
#include <CGAL/Triangulation_vertex.h>

namespace ttk::gph {

  using rpd::Diagram;
  using rpd::inf;
  using rpd::MultidimensionalDiagram;
  using rpd::Simplex;
  using rpd::UnionFind;

  constexpr static unsigned DYN_DIM = 0;

  template <unsigned D>
  using DSimplex = std::
    conditional_t<D == DYN_DIM, std::vector<id_t>, std::array<id_t, D + 1>>;

  template <unsigned D>
  using ValueArray = std::conditional_t<D == DYN_DIM,
                                        std::vector<rpd::value_t>,
                                        std::array<rpd::value_t, D>>;

  template <unsigned D>
  using ConnectivityHashMap
    = HashMap<DSimplex<D>, std::vector<std::pair<int, int>>>;

#ifdef TTK_CONCURRENT_HASHTABLE_AVAILABLE
  template <unsigned D>
  using ConcurrentConnectivityHashMap
    = ConcurrentHashMap<DSimplex<D>, std::vector<std::pair<int, int>>>;
#endif

  template <unsigned D>
  struct FiltratedDSimplex {
    DSimplex<D> s;
    double d;
    double a;
  };

  template <unsigned D>
  bool operator<(FiltratedDSimplex<D> const &s1,
                 FiltratedDSimplex<D> const &s2) {
    return std::make_pair(s1.d, s1.a) < std::make_pair(s2.d, s2.a);
  }

  template <unsigned D>
  FiltratedDSimplex<D> max(FiltratedDSimplex<D> const &s1,
                           FiltratedDSimplex<D> const &s2) {
    return s1 < s2 ? s2 : s1;
  }

  template <unsigned D>
  using DSimplicialComplex = std::vector<FiltratedDSimplex<D>>;

  template <unsigned DIM>
  class DRPersistenceD {
    using DimTag = std::conditional_t<DIM == DYN_DIM,
                                      CGAL::Dynamic_dimension_tag,
                                      CGAL::Dimension_tag<DIM>>;
    using K = CGAL::Epick_d<DimTag>;
    using Vb = CGAL::Triangulation_vertex<K, id_t>;
    using FCb = CGAL::Triangulation_full_cell<K, id_t>;
    using Tds = CGAL::Triangulation_data_structure<DimTag, Vb, FCb>;
    using Delaunay = CGAL::Delaunay_triangulation<K, Tds>;
    using Point = typename Delaunay::Point;
    using Facet = typename Delaunay::Facet;
    using VertexHandle = typename Delaunay::Vertex_handle;
    using CellHandle = typename Delaunay::Full_cell_handle;

    static constexpr unsigned minus1(const unsigned D) {
      return D == DYN_DIM ? DYN_DIM : D - 1;
    }

    static constexpr unsigned plus1(const unsigned D) {
      return D == DYN_DIM ? DYN_DIM : D + 1;
    }

    struct FiltratedQuadFacet {
      DSimplex<minus1(DIM)> s;
      double d;
      double a;
      int c1;
      int c2;
    };

  public:
    explicit DRPersistenceD(PointCloud<DIM> &points, unsigned nThreads = 1);
    void run(MultidimensionalDiagram &ph);

  private:
    const unsigned N_p;
    unsigned N_c{};
    Delaunay del_;
    PointCloud<DIM> &points_;
    const int nThreads_{1};

    [[nodiscard]] double squaredDistance(const unsigned i1,
                                         const unsigned i2) const {
      return std::inner_product(
        points_[i1].begin(), points_[i1].end(), points_[i2].begin(), 0.,
        std::plus(),
        [](const double u, const double v) { return (u - v) * (u - v); });
    }

    template <unsigned D>
    std::pair<double, double>
      squaredPerturbedDiameter(DSimplex<D> const &s) const {
      double diam = 0.;
      double a = 0.;
      for(unsigned p1 = 0; p1 < s.size(); ++p1) {
        for(unsigned p2 = p1 + 1; p2 < s.size(); ++p2) {
          const double dist = squaredDistance(s[p1], s[p2]);
          diam = std::max(diam, dist);
          a += dist;
        }
      }
      return {diam, a};
    }

    template <unsigned D>
    void getLengths(DSimplex<D> const &s,
                    ValueArray<D *(D + 1) / 2> &lengths) const {
      if constexpr(D == DYN_DIM)
        lengths.resize((s.size() - 1) * s.size() / 2);
      else
        lengths.fill(0.); // to mute an uninitialization warning
      int i = 0;
      for(unsigned p1 = 0; p1 < s.size(); ++p1) {
        for(unsigned p2 = p1 + 1; p2 < s.size(); ++p2)
          lengths[i++] = squaredDistance(s[p1], s[p2]);
      }
    }

    template <unsigned D>
    bool checkLinkUrquhart(DSimplex<D> const &s,
                           ValueArray<D *(D + 1) / 2> const &lengths,
                           std::pair<double, double> diam,
                           id_t linkId) const {
      ValueArray<plus1(D)> linkLengths;
      if constexpr(D == DYN_DIM)
        linkLengths.resize(s.size());
      for(unsigned i = 0; i < s.size(); ++i)
        linkLengths[i] = squaredDistance(s[i], linkId);

      for(unsigned i = 0; i < s.size(); ++i) {
        double d = 0.;
        double a = 0.;
        int k = 0;
        for(unsigned p1 = 0; p1 < s.size(); ++p1) {
          for(unsigned p2 = p1 + 1; p2 < s.size(); ++p2) {
            double dist;
            if(p1 == i)
              dist = linkLengths[p2];
            else if(p2 == i)
              dist = linkLengths[p1];
            else
              dist = lengths[k];

            d = std::max(d, dist);
            a += dist;
            ++k;
          }
        }
        if(std::make_pair(d, a) > diam)
          return false;
      }

      return true;
    }

    void computeDelaunay();

    void computeDPH(Diagram &ph, DSimplicialComplex<minus1(DIM)> &MSA);

    void computeDPH_p(Diagram &ph, DSimplicialComplex<minus1(DIM)> &MSA);

    template <unsigned D>
    void computeNextPH(Diagram &ph,
                       DSimplicialComplex<D> &MSA,
                       DSimplicialComplex<minus1(D)> &nextMSA) const;

    template <unsigned D>
    void computeNextPH_p(Diagram &ph,
                         DSimplicialComplex<D> &MSA,
                         DSimplicialComplex<minus1(D)> &nextMSA) const;

    template <unsigned D>
    void recurse(MultidimensionalDiagram &ph, DSimplicialComplex<D> &MSA) const;
  };

  template <unsigned DIM>
  DRPersistenceD<DIM>::DRPersistenceD(PointCloud<DIM> &points,
                                      const unsigned nThreads)
    : N_p(points.size()), del_(DIM), points_(points), nThreads_(nThreads) {
  }

  template <>
  inline DRPersistenceD<DYN_DIM>::DRPersistenceD(PointCloud<DYN_DIM> &points,
                                                 const unsigned nThreads)
    : N_p(points.size()), del_(points[0].size()), points_(points),
      nThreads_(nThreads) {
  }

  template <unsigned DIM>
  void DRPersistenceD<DIM>::run(MultidimensionalDiagram &ph) {
    ph = MultidimensionalDiagram(DIM);

    computeDelaunay();

    DSimplicialComplex<DIM - 1> MSA;
    if(nThreads_ == 1)
      computeDPH(ph[DIM - 1], MSA);
    else
      computeDPH_p(ph[DIM - 1], MSA);

    recurse(ph, MSA);
  }

  template <>
  inline void DRPersistenceD<DYN_DIM>::run(MultidimensionalDiagram &ph) {
    const unsigned DIM = points_[0].size();
    ph = MultidimensionalDiagram(DIM);

    computeDelaunay();

    DSimplicialComplex<DYN_DIM> MSA;
    if(nThreads_ == 1)
      computeDPH(ph[DIM - 1], MSA);
    else
      computeDPH_p(ph[DIM - 1], MSA);

    unsigned D = DIM - 1;
    while(D > 1) {
      DSimplicialComplex<DYN_DIM> nextMSA;
      if(nThreads_ == 1)
        computeNextPH<DYN_DIM>(ph[D - 1], MSA, nextMSA);
      else
        computeNextPH_p<DYN_DIM>(ph[D - 1], MSA, nextMSA);
      MSA = std::move(nextMSA);
      D--;
    }
    for(FiltratedDSimplex<DYN_DIM> const &e : MSA)
      ph[0].emplace_back(
        FiltratedSimplex{{-1}, 0.}, FiltratedSimplex{e.s, sqrt(e.d)});
  }

  template <unsigned DIM>
  template <unsigned D>
  void DRPersistenceD<DIM>::recurse(MultidimensionalDiagram &ph,
                                    DSimplicialComplex<D> &MSA) const {
    if constexpr(D >= 2) {
      DSimplicialComplex<D - 1> nextMSA;
      if(nThreads_ == 1)
        computeNextPH(ph[D - 1], MSA, nextMSA);
      else
        computeNextPH_p(ph[D - 1], MSA, nextMSA);
      recurse(ph, nextMSA);
    } else if constexpr(D == 1) { // this is the 0-dimensional homology
      for(FiltratedDSimplex<1> const &e : MSA)
        ph[0].emplace_back(FiltratedSimplex{{-1}, 0.},
                           FiltratedSimplex{{e.s[0], e.s[1]}, sqrt(e.d)});
    }
  }

  template <unsigned DIM>
  void DRPersistenceD<DIM>::computeDelaunay() {
    del_.insert(points_.begin(), points_.end());

    // index all vertices
    int p = 0;
    for(auto p_it = del_.finite_vertices_begin();
        p_it != del_.finite_vertices_end(); ++p_it) {
      points_[p] = p_it->point();
      p_it->data() = p++;
    }

    // index all cells (including infinite ones)
    int k = 0;
    for(auto c_it = del_.full_cells_begin(); c_it != del_.full_cells_end();
        ++c_it)
      c_it->data() = k++;
    N_c = k;
  }

  template <unsigned DIM>
  void DRPersistenceD<DIM>::computeDPH(Diagram &ph,
                                       DSimplicialComplex<minus1(DIM)> &MSA) {
    DSimplicialComplex<minus1(DIM)> maxDelaunay(N_c);
    UnionFind UF(N_c);
    std::vector<FiltratedQuadFacet> hyperUrquhart;

    /* computing urquhart hypergraph (codimension 1) */
    for(auto f_it = del_.facets_begin(); f_it != del_.facets_end(); ++f_it) {

      const Facet f = *f_it;
      const CellHandle c = f.full_cell();
      const CellHandle c_mirror = c->neighbor(f.index_of_covertex());

      if(del_.is_infinite(f))
        UF.merge(c->data(), c_mirror->data());

      else { // we need to determine whether f is "Urquhart"
        const auto linkPoint1 = c->vertex(f.index_of_covertex());
        const auto linkPoint2
          = c_mirror->vertex(c->mirror_index(f.index_of_covertex()));

        DSimplex<minus1(DIM)> facet;
        if constexpr(DIM == DYN_DIM)
          facet.resize(del_.current_dimension());
        for(unsigned i = 0; i < facet.size(); ++i)
          facet[i]
            = c->vertex((f.index_of_covertex() + i + 1) % (facet.size() + 1))
                ->data();
        std::sort(facet.begin(), facet.end());
        const auto diam = squaredPerturbedDiameter<minus1(DIM)>(facet);

        bool is_urquhart = true;
        for(auto linkPoint : {linkPoint1, linkPoint2}) {
          if(!del_.is_infinite(linkPoint)) {
            const id_t k = linkPoint->data();
            bool largest = true;
            for(unsigned i = 0; i < facet.size(); ++i) {
              DSimplex<minus1(DIM)> neighbor = facet;
              neighbor[i] = k;
              if(squaredPerturbedDiameter<minus1(DIM)>(neighbor) > diam) {
                largest = false;
                break;
              }
            }
            if(largest) {
              is_urquhart = false;
              break;
            }
          }
        }

        if(is_urquhart)
          hyperUrquhart.push_back(
            {facet, diam.first, diam.second, c->data(), c_mirror->data()});
        else {
          const int poly1 = UF.find(c->data());
          const int poly2 = UF.find(c_mirror->data());
          maxDelaunay[UF.mergeRet(poly1, poly2)] = max(
            FiltratedDSimplex<minus1(DIM)>{facet, diam.first, diam.second},
            max(maxDelaunay[poly1], maxDelaunay[poly2]));
        }

        // detect infinite polyhedrons (going beyond convex hull)
        if(del_.is_infinite(c))
          maxDelaunay[UF.find(c->data())].d = inf;
        else if(del_.is_infinite(c_mirror))
          maxDelaunay[UF.find(c_mirror->data())].d = inf;
      }
    }

    std::sort(hyperUrquhart.begin(), hyperUrquhart.end(),
              [](const FiltratedQuadFacet &f1, const FiltratedQuadFacet &f2) {
                if(f1.d == f2.d)
                  return f1.a > f2.a;
                return f1.d > f2.d;
              });

    /* reverse-delete algorithm to determine MSA */
    std::vector<int> latest(maxDelaunay.size());
    std::iota(latest.begin(), latest.end(), 0);

    for(FiltratedQuadFacet const &f :
        hyperUrquhart) { // sorted by decreasing order
      const int v1 = UF.find(f.c1);
      const int v2 = UF.find(f.c2);
      if(v1 != v2) { // two distinct codimension-1 cavities: merge them by
                     // deleting the facet
        UF.merge(v1, v2);

        const int latest1 = latest[v1];
        const int latest2 = latest[v2];
        const FiltratedDSimplex<minus1(DIM)> &death1 = maxDelaunay[latest1];
        const FiltratedDSimplex<minus1(DIM)> &death2 = maxDelaunay[latest2];

        if(death1.d < death2.d) {
          if(f.d < death1.d)
            ph.emplace_back(
              FiltratedSimplex{Simplex(f.s.begin(), f.s.end()), sqrt(f.d)},
              FiltratedSimplex{
                Simplex(death1.s.begin(), death1.s.end()), sqrt(death1.d)});
          latest[UF.find(v1)] = latest2;
        } else if(death2.d < death1.d) {
          if(f.d < death2.d)
            ph.emplace_back(
              FiltratedSimplex{Simplex(f.s.begin(), f.s.end()), sqrt(f.d)},
              FiltratedSimplex{
                Simplex(death2.s.begin(), death2.s.end()), sqrt(death2.d)});
          latest[UF.find(v1)] = latest1;
        }
      } else // this is a facet from the minimal spanning acycle
        MSA.push_back({f.s, f.d, f.a});
    }
  }

  template <unsigned DIM>
  void DRPersistenceD<DIM>::computeDPH_p(Diagram &ph,
                                         DSimplicialComplex<minus1(DIM)> &MSA) {
#ifndef TTK_GPH_PARALLEL
    computeDPH(ph, MSA);
#else
    omp_set_num_threads(nThreads_);

    DisjointSets UF(N_c);
    tbb::concurrent_vector<FiltratedQuadFacet> hyperUrquhart;

    std::vector<FiltratedDSimplex<DIM>> cells(N_c);
    for(auto c_it = del_.finite_full_cells_begin();
        c_it != del_.finite_full_cells_end(); ++c_it) {
      FiltratedDSimplex<DIM> cell;
      if constexpr(DIM == DYN_DIM)
        cell.s.resize(del_.current_dimension() + 1);
      cell.d = -1.;
      for(unsigned i = 0; i < cell.s.size(); ++i)
        cell.s[i] = c_it->vertex(i)->data();
      cells[c_it->data()] = cell;
    }

#pragma omp parallel for
    for(unsigned i = 0; i < cells.size(); ++i) {
      if(cells[i].d == -1.) {
        auto [d, a] = squaredPerturbedDiameter<DIM>(cells[i].s);
        cells[i].d = d;
        cells[i].a = a;
      }
    }

    std::vector<std::tuple<DSimplex<minus1(DIM)>, int, int, int, int>> facets;
    for(auto f_it = del_.facets_begin(); f_it != del_.facets_end(); ++f_it) {
      const Facet f = *f_it;
      const CellHandle c = f.full_cell();
      const CellHandle c_mirror = c->neighbor(f.index_of_covertex());
      if(del_.is_infinite(f))
        UF.unite(c->data(), c_mirror->data());
      else {
        DSimplex<minus1(DIM)> facet;
        if constexpr(DIM == DYN_DIM)
          facet.resize(del_.current_dimension());
        for(unsigned i = 0; i < facet.size(); ++i)
          facet[i]
            = c->vertex((f.index_of_covertex() + i + 1) % (facet.size() + 1))
                ->data();

        const auto linkPoint1 = c->vertex(f.index_of_covertex());
        const auto linkPoint2
          = c_mirror->vertex(c->mirror_index(f.index_of_covertex()));
        const int n1 = del_.is_infinite(linkPoint1) ? -1 : linkPoint1->data();
        const int n2 = del_.is_infinite(linkPoint2) ? -1 : linkPoint2->data();
        const int c1 = c->data();
        const int c2 = c_mirror->data();
        facets.emplace_back(facet, c1, n1, c2, n2);
      }
    }

#pragma omp parallel for
    for(auto &[facet, c1, n1, c2, n2] : facets) {
      // first determine whether s is Urquhart
      std::sort(facet.begin(), facet.end());
      const auto diam = squaredPerturbedDiameter<minus1(DIM)>(facet);

      bool is_urquhart = true;
      for(const id_t k : {n1, n2}) {
        if(k != -1) {
          bool largest = true;
          for(unsigned i = 0; i < facet.size(); ++i) {
            DSimplex<minus1(DIM)> neighbor = facet;
            neighbor[i] = k;
            if(squaredPerturbedDiameter<minus1(DIM)>(neighbor) > diam) {
              largest = false;
              break;
            }
          }
          if(largest) {
            is_urquhart = false;
            break;
          }
        }
      }

      if(is_urquhart)
        hyperUrquhart.push_back({facet, diam.first, diam.second, c1, c2});
      else
        UF.unite(c1, c2);

      if(n1 == -1)
        cells[c1].d = inf;
      else if(n2 == -1)
        cells[c2].d = inf;
    }

    std::vector<std::mutex> maxDelaunayLocks(N_c);
#pragma omp parallel for
    for(unsigned x = 0; x < N_c; ++x) {
      const int poly = UF.find(x);
      std::lock_guard lock(maxDelaunayLocks[poly]);
      cells[poly] = max(cells[poly], cells[x]);
    }

    TTK_PSORT(nThreads_, hyperUrquhart.begin(), hyperUrquhart.end(),
              [](const FiltratedQuadFacet &f1, const FiltratedQuadFacet &f2) {
                if(f1.d == f2.d)
                  return f1.a > f2.a;
                return f1.d > f2.d;
              });

    /* reverse-delete algorithm to determine MSA */
    std::vector<int> latest(N_c);
    std::iota(latest.begin(), latest.end(), 0);

    for(FiltratedQuadFacet const &f :
        hyperUrquhart) { // sorted by decreasing order
      const int v1 = UF.find(f.c1);
      const int v2 = UF.find(f.c2);
      if(v1 != v2) { // two distinct codimension-1 cavities: merge them by
                     // deleting the facet
        UF.unite(v1, v2);

        const int latest1 = latest[v1];
        const int latest2 = latest[v2];
        const FiltratedDSimplex<DIM> &death1 = cells[latest1];
        const FiltratedDSimplex<DIM> &death2 = cells[latest2];

        if(death1.d < death2.d) {
          if(f.d < death1.d)
            ph.emplace_back(
              FiltratedSimplex{Simplex(f.s.begin(), f.s.end()), sqrt(f.d)},
              FiltratedSimplex{
                Simplex(death1.s.begin(), death1.s.end()), sqrt(death1.d)});
          latest[UF.find(v1)] = latest2;
        } else if(death2.d < death1.d) {
          if(f.d < death2.d)
            ph.emplace_back(
              FiltratedSimplex{Simplex(f.s.begin(), f.s.end()), sqrt(f.d)},
              FiltratedSimplex{
                Simplex(death2.s.begin(), death2.s.end()), sqrt(death2.d)});
          latest[UF.find(v1)] = latest1;
        }
      } else // this is a facet from the minimal spanning acycle
        MSA.push_back({f.s, f.d, f.a});
    }
#endif
  }

  template <unsigned DIM>
  template <unsigned D>
  void DRPersistenceD<DIM>::computeNextPH(
    Diagram &ph,
    DSimplicialComplex<D> &MSA,
    DSimplicialComplex<minus1(D)> &nextMSA) const {
    const unsigned N_msa = MSA.size();

    /* Connectivity */

    ConnectivityHashMap<minus1(D)> msa_connectivity;
    msa_connectivity.reserve(N_msa);
    for(unsigned i = 0; i < N_msa; ++i) {
      DSimplex<minus1(D)> face;
      if constexpr(D == DYN_DIM)
        face.resize(MSA[0].s.size() - 1);
      for(unsigned k = 0; k < MSA[0].s.size(); ++k) {
        for(unsigned j = 0; j < MSA[0].s.size() - 1; ++j)
          face[j] = MSA[i].s[j + (j >= k)];
        msa_connectivity[face].reserve(4); // guess
        msa_connectivity[face].emplace_back(i, MSA[i].s[k]);
      }
    }

    /* Urquhart-ness and Urquhart-polytopes */

    std::vector<FiltratedDSimplex<minus1(D)>> critical;
    UnionFind UF_msa(N_msa);

    for(auto const &[s, neighbors] : msa_connectivity) {

      // first determine whether s is Urquhart
      bool is_urquhart = true;
      ValueArray<D * minus1(D) / 2> lengths;
      getLengths<minus1(D)>(s, lengths);
      std::pair diam{0.0, 0.0};
      for(double const &l : lengths)
        diam = {std::max(diam.first, l), diam.second + l};

      for(auto [coface_id, linkPoint_id] : neighbors) {
        if(checkLinkUrquhart<minus1(D)>(s, lengths, diam, linkPoint_id)) {
          is_urquhart = false;
          break;
        }
      }

      // now maintain the polytope structure
      if(is_urquhart)
        critical.push_back({s, diam.first, diam.second});
      else {
        if(neighbors.size() == 1) // set infinite polytope
          MSA[neighbors[0].first].d = inf;
        else if(neighbors.size() == 2)
          UF_msa.merge(neighbors[0].first, neighbors[1].first);
        else
          critical.push_back({s, diam.first, diam.second});
      }
    }

    for(unsigned x = 0; x < N_msa; ++x) {
      const int poly = UF_msa.find(x);
      MSA[poly] = max(MSA[poly], MSA[x]);
    }

    /* Graph critical -- polytope */

    std::vector<id_t> polytopes;
    for(unsigned x = 0; x < N_msa; ++x) {
      if(UF_msa.isRoot(x) && MSA[x].d < inf)
        polytopes.emplace_back(x);
    }

    std::sort(polytopes.begin(), polytopes.end(),
              [&](const int x1, const int x2) { return MSA[x1] < MSA[x2]; });

    std::vector<int> criticalIndices(critical.size());
    std::iota(criticalIndices.begin(), criticalIndices.end(), 0);
    std::sort(criticalIndices.begin(), criticalIndices.end(),
              [&](const int x1, const int x2) {
                return critical[x1].d < critical[x2].d;
              });
    std::vector<int> criticalOrder(critical.size());
    for(unsigned i = 0; i < criticalIndices.size(); ++i)
      criticalOrder[criticalIndices[i]] = i;

    std::vector<std::vector<int>> poly_to_crit(N_msa);
    for(const int poly : polytopes)
      poly_to_crit[poly].reserve(D + 1); // guess
    for(unsigned i = 0; i < critical.size(); ++i) {
      for(const auto &[poly, _] : msa_connectivity[critical[i].s]) {
        if(MSA[UF_msa.find(poly)].d < inf) {
          auto &neighbors = poly_to_crit[UF_msa.find(poly)];
          auto it
            = std::find(neighbors.begin(), neighbors.end(), criticalOrder[i]);
          if(it == neighbors.end())
            neighbors.push_back(criticalOrder[i]);
          else
            neighbors.erase(it);
        }
      }
    }

    /* PairCells */

    std::vector<int> partner(critical.size(), -1);
    for(const int poly : polytopes) {
      HashSet<int> boundary(poly_to_crit[poly].begin(),
                            poly_to_crit[poly].end(),
                            poly_to_crit[poly].size());
      while(true) {
        const int youngest_id
          = *std::max_element(boundary.begin(), boundary.end());
        if(partner[youngest_id] == -1) {
          partner[youngest_id] = poly;
          const FiltratedDSimplex<minus1(D)> &c
            = critical[criticalIndices[youngest_id]];
          const FiltratedDSimplex<D> &death = MSA[poly];
          if(c.d < death.d)
            ph.emplace_back(
              FiltratedSimplex{Simplex(c.s.begin(), c.s.end()), sqrt(c.d)},
              FiltratedSimplex{
                Simplex(death.s.begin(), death.s.end()), sqrt(death.d)});
          break;
        } else {
          for(const int crit : poly_to_crit[partner[youngest_id]]) {
            const auto it = boundary.find(crit);
            if(it == boundary.end())
              boundary.insert(crit);
            else
              boundary.erase(it);
          }
        }
      }
      poly_to_crit[poly].assign(boundary.begin(), boundary.end());
    }

    /* Next MSA */
    for(unsigned i = 0; i < critical.size(); ++i) {
      if(partner[i] == -1) // unassigned -> go in next MSA
        nextMSA.emplace_back(critical[criticalIndices[i]]);
    }
  }

  template <unsigned DIM>
  template <unsigned D>
  void DRPersistenceD<DIM>::computeNextPH_p(
    Diagram &ph,
    DSimplicialComplex<D> &MSA,
    DSimplicialComplex<minus1(D)> &nextMSA) const {
#ifndef TTK_GPH_PARALLEL
    computeNextPH(ph, MSA, nextMSA);
#else
    tbb::global_control gc(
      tbb::global_control::max_allowed_parallelism, nThreads_);
    omp_set_num_threads(nThreads_);
    const unsigned N_msa = MSA.size();

    /* Connectivity */

    ConcurrentConnectivityHashMap<minus1(D)> concurrent_msa_connectivity;
    concurrent_msa_connectivity.reserve(N_msa);

#pragma omp parallel for
    for(unsigned i = 0; i < N_msa; ++i) {
      DSimplex<minus1(D)> face;
      if constexpr(D == DYN_DIM)
        face.resize(MSA[0].s.size() - 1);
      for(unsigned k = 0; k < MSA[0].s.size(); ++k) {
        for(unsigned j = 0; j < MSA[0].s.size() - 1; ++j)
          face[j] = MSA[i].s[j + (j >= k)];
        concurrent_msa_connectivity.emplace_or_visit(
          face,
          std::vector<std::pair<int, int>>{std::make_pair(i, MSA[i].s[k])},
          [&](auto &x) { x.second.emplace_back(i, MSA[i].s[k]); });
      }
    }

    /* Urquhart-ness and Urquhart-polytopes */

    tbb::concurrent_vector<FiltratedDSimplex<minus1(D)>> critical;
    DisjointSets UF_msa(N_msa);

    concurrent_msa_connectivity.cvisit_all(
#ifdef __cpp_lib_execution
      std::execution::par,
#endif
      [&](const auto &x) {
        const auto &[s, neighbors] = x;

        // first determine whether s is Urquhart
        bool is_urquhart = true;
        ValueArray<D * minus1(D) / 2> lengths;
        getLengths<minus1(D)>(s, lengths);
        std::pair diam{0.0, 0.0};
        for(double const &l : lengths)
          diam = {std::max(diam.first, l), diam.second + l};

        for(auto [coface_id, linkPoint_id] : neighbors) {
          if(checkLinkUrquhart<minus1(D)>(s, lengths, diam, linkPoint_id)) {
            is_urquhart = false;
            break;
          }
        }

        // now maintain the polytope structure
        if(is_urquhart)
          critical.push_back({s, diam.first, diam.second});
        else {
          if(neighbors.size() == 1) // set infinite polytope
            MSA[neighbors[0].first].d = inf;
          else if(neighbors.size() == 2)
            UF_msa.unite(neighbors[0].first, neighbors[1].first);
          else
            critical.push_back({s, diam.first, diam.second});
        }
      });

    std::vector<std::mutex> maxDelaunayLocks(N_msa);
#pragma omp parallel for
    for(unsigned x = 0; x < N_msa; ++x) {
      const int poly = UF_msa.find(x);
      std::lock_guard lock(maxDelaunayLocks[poly]);
      MSA[poly] = max(MSA[poly], MSA[x]);
    }

    const ConnectivityHashMap<minus1(D)> msa_connectivity
      = std::move(concurrent_msa_connectivity);

    /* Graph critical -- polytope */

    std::vector<id_t> polytopes;
    for(unsigned x = 0; x < N_msa; ++x) {
      if(UF_msa.isRoot(x) && MSA[x].d < inf)
        polytopes.emplace_back(x);
    }

    TTK_PSORT(nThreads_, polytopes.begin(), polytopes.end(),
              [&](const int x1, const int x2) { return MSA[x1] < MSA[x2]; });

    std::vector<int> criticalIndices(critical.size());
    std::iota(criticalIndices.begin(), criticalIndices.end(), 0);
    TTK_PSORT(nThreads_, criticalIndices.begin(), criticalIndices.end(),
              [&](const int x1, const int x2) {
                return critical[x1].d < critical[x2].d;
              });
    std::vector<int> criticalOrder(critical.size());
    for(unsigned i = 0; i < criticalIndices.size(); ++i)
      criticalOrder[criticalIndices[i]] = i;

    std::vector<std::vector<int>> poly_to_crit(N_msa);
    std::vector<std::mutex> poly_mutex(N_msa);
    std::vector<std::mutex> crit_mutex(critical.size());

#pragma omp parallel for
    for(unsigned i = 0; i < critical.size(); ++i) {
      for(const auto &[poly, _] :
          msa_connectivity.find(critical[i].s)->second) {
        const int uf_poly = UF_msa.find(poly);
        if(MSA[uf_poly].d < inf) {
          auto &neighbors = poly_to_crit[uf_poly];
          std::lock_guard lock(poly_mutex[uf_poly]);
          auto it
            = std::find(neighbors.begin(), neighbors.end(), criticalOrder[i]);
          if(it == neighbors.end())
            neighbors.push_back(criticalOrder[i]);
          else
            neighbors.erase(it);
        }
      }
    }

    /* PairCells */

    std::vector<int> partner(critical.size(), -1);

    auto eliminateBoundary = [&](const int poly) -> int {
      std::lock_guard lock(poly_mutex[poly]);
      HashSet<int> boundary(poly_to_crit[poly].begin(),
                            poly_to_crit[poly].end(),
                            poly_to_crit[poly].size());
      ;
      while(true) {
        const int youngest_id
          = *std::max_element(boundary.begin(), boundary.end());
        std::lock_guard lock_crit(crit_mutex[youngest_id]);
        if(partner[youngest_id] == -1) {
          partner[youngest_id] = poly;
          poly_to_crit[poly].assign(boundary.begin(), boundary.end());
          return -1;
        }
        std::lock_guard lock_partner(poly_mutex[partner[youngest_id]]);
        if(MSA[partner[youngest_id]].d > MSA[poly].d) {
          const int tmp = partner[youngest_id];
          partner[youngest_id] = poly;
          poly_to_crit[poly].assign(boundary.begin(), boundary.end());
          return tmp;
        }
        for(const int crit : poly_to_crit[partner[youngest_id]]) {
          const auto it = boundary.find(crit);
          if(it == boundary.end())
            boundary.insert(crit);
          else
            boundary.erase(it);
        }
      }
    };

#pragma omp parallel for
    for(const int poly : polytopes) {
      int p = poly;
      while(p >= 0)
        p = eliminateBoundary(p);
    }

    /* Next MSA */
    for(unsigned i = 0; i < critical.size(); ++i) {
      const int poly = partner[i];
      if(poly == -1) // unassigned -> go in next MSA
        nextMSA.emplace_back(critical[criticalIndices[i]]);
      else {
        const FiltratedDSimplex<minus1(D)> &c = critical[criticalIndices[i]];
        const FiltratedDSimplex<D> &death = MSA[poly];
        if(c.d < death.d)
          ph.emplace_back(
            FiltratedSimplex{Simplex(c.s.begin(), c.s.end()), sqrt(c.d)},
            FiltratedSimplex{
              Simplex(death.s.begin(), death.s.end()), sqrt(death.d)});
      }
    }
#endif
  }

  template <unsigned DIM>
  void runDelaunayRipsPersistenceDiagram(rpd::PointCloud const &points,
                                         MultidimensionalDiagram &diagram,
                                         const int threads) {
    PointCloud<DIM> p(points.size());
    for(unsigned i = 0; i < points.size(); ++i) {
      for(unsigned d = 0; d < DIM; ++d)
        p[i][d] = points[i][d];
    }
    DRPersistenceD<DIM> drpd(p, threads);
    drpd.run(diagram);
  }

  template <>
  inline void
    runDelaunayRipsPersistenceDiagram<DYN_DIM>(rpd::PointCloud const &points,
                                               MultidimensionalDiagram &diagram,
                                               const int threads) {
    PointCloud<DYN_DIM> p = points;
    DRPersistenceD<DYN_DIM> drpd(p, threads);
    drpd.run(diagram);
  }

  template <unsigned DIM>
  void tryDimension(rpd::PointCloud const &points,
                    MultidimensionalDiagram &diagram,
                    const int threads) {
    if constexpr(DIM <= TTK_DELAUNAY_MAX_COMPILED_DIMENSION) {
      if(points[0].size() == DIM)
        runDelaunayRipsPersistenceDiagram<DIM>(points, diagram, threads);
      else
        tryDimension<DIM + 1>(points, diagram, threads);
    } else
      runDelaunayRipsPersistenceDiagram<DYN_DIM>(points, diagram, threads);
  }

  inline void tryDimensions(rpd::PointCloud const &points,
                            MultidimensionalDiagram &diagram,
                            const int threads) {
    tryDimension<4>(points, diagram, threads);
  }

} // namespace ttk::gph

#endif