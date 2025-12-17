#pragma once

#include <geoPHUtils.h>

#ifdef TTK_ENABLE_CGAL

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

namespace ttk::gph {

  using rpd::Edge;
  using rpd::FiltratedEdge;
  using rpd::FiltratedSimplex;
  using rpd::Generator1;
  using rpd::Generator2;
  using rpd::inf;
  using rpd::MultidimensionalDiagram;
  using rpd::UnionFind;

  using Facet = std::array<id_t, 3>;

  struct FiltratedFacet {
    Facet f;
    double d;
  };

  inline FiltratedFacet max(FiltratedFacet const &a, FiltratedFacet const &b) {
    if(a.d > b.d)
      return a;
    return b;
  }

  struct FiltratedQuadFacet {
    Facet f;
    int c1;
    int c2;
    double d;
    double a;
  };

  class DRPersistence3 {
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Vb = CGAL::Triangulation_vertex_base_with_info_3<int, K>;
    using Fb = CGAL::Triangulation_cell_base_with_info_3<int, K>;
    using Tds = CGAL::Triangulation_data_structure_3<Vb, Fb>;
    using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;
    using Point = Delaunay::Point_3;

    using AdjacencyList = std::vector<int>;
    using ConnectivityHashMap = HashMap<Edge, std::pair<AdjacencyList, bool>>;

  public:
    explicit DRPersistence3(const PointCloud<3> &points)
      : N_p(points.size()), p(points) {
    }

    /**
     * Computes the Delaunay-Rips persistence diagram of the point cloud given
     * to the constructor
     */
    void run(MultidimensionalDiagram &ph) {
      ph = MultidimensionalDiagram(3);

      computeDelaunay();

      UnionFind UF(N_c);
      computeFacetUrquhart(UF);
      compute2PH(UF, ph);

      connectivity();
      computePolygons(ph);
      compute1PH(ph);
    }

    /**
     * Computes the Delaunay-Rips persistence diagram and associated 1- and
     * 2-dimensional generators of the point cloud given to the constructor
     */
    void run(MultidimensionalDiagram &ph,
             std::vector<Generator1> &generators1,
             std::vector<Generator2> &generators2) {
      ph = MultidimensionalDiagram(3);
      generators1.resize(0);
      generators2.resize(0);

      computeDelaunay();

      UnionFind UF(N_c);
      computeFacetUrquhart(UF);
      compute2PH(UF, ph, generators2);

      connectivity();
      computePolygons(ph);
      compute1PH(ph, generators1);
    }

  private:
    const unsigned N_p;
    unsigned N_c{0};
    const PointCloud<3> &p;
    Delaunay del;

    // 1-dimensional
    std::vector<FiltratedEdge> urquhart;
    std::vector<FiltratedEdge> critical1;

    ConnectivityHashMap msa_edges{};

    UnionFind UF_msa{0};

    // 2-dimensional
    std::vector<FiltratedQuadFacet> hyperUrquhart;
    std::vector<FiltratedFacet> msa;
    std::vector<FiltratedFacet> maxDelaunay;

    [[nodiscard]] double squaredDistance(const unsigned i1,
                                         const unsigned i2) const {
      PointD<3> const p1 = p[i1], p2 = p[i2];
      return (p1[0] - p2[0]) * (p1[0] - p2[0])
             + (p1[1] - p2[1]) * (p1[1] - p2[1])
             + (p1[2] - p2[2]) * (p1[2] - p2[2]);
    }
    [[nodiscard]] static std::pair<double, double>
      perturbedDiameter(const double d1, const double d2, const double d3) {
      return {std::max(d1, std::max(d2, d3)), d1 + d2 + d3};
    }

    /**
     * compute Delaunay and initialize tetrahedra IDs (including infinite
     * tetrahedra to deal with the convex hull)
     */
    void computeDelaunay() {
      std::vector<std::pair<Point, unsigned>> points(N_p);
      for(unsigned i = 0; i < N_p; ++i)
        points[i] = std::make_pair(Point(p[i][0], p[i][1], p[i][2]), i);
      del = Delaunay(points.begin(), points.end());
      int k = 0;
      for(auto const &c : del.all_cell_handles())
        c->info() = k++;
      N_c = k;
    }

    void computeFacetUrquhart(UnionFind &UF) {
      maxDelaunay.resize(N_c, FiltratedFacet{{-1, -1, -1}, 0.});
      for(Delaunay::Facet const &f : del.all_facets()) {
        if(del.is_infinite(f))
          UF.merge(f.first->info(), del.mirror_facet(f).first->info());
        else {
          Facet facet{f.first->vertex((f.second + 1) % 4)->info(),
                      f.first->vertex((f.second + 2) % 4)->info(),
                      f.first->vertex((f.second + 3) % 4)->info()};
          const double d01 = squaredDistance(facet[0], facet[1]);
          const double d02 = squaredDistance(facet[0], facet[2]);
          const double d12 = squaredDistance(facet[1], facet[2]);
          const auto diam = perturbedDiameter(d01, d02, d12);
          const Delaunay::Facet f_m = del.mirror_facet(f);

          // check if facet is Urquhart
          bool is_urquhart = true;
          if(!del.is_infinite(f.first->vertex(f.second))) {
            const unsigned k = f.first->vertex(f.second)->info();
            const double d0k = squaredDistance(facet[0], k);
            const double d1k = squaredDistance(facet[1], k);
            const double d2k = squaredDistance(facet[2], k);
            if(perturbedDiameter(d01, d0k, d1k) < diam
               && perturbedDiameter(d02, d0k, d2k) < diam
               && perturbedDiameter(d12, d1k, d2k) < diam)
              is_urquhart = false;
          }
          if(is_urquhart && !del.is_infinite(f_m.first->vertex(f_m.second))) {
            const unsigned l = f_m.first->vertex(f_m.second)->info();
            const double d0l = squaredDistance(facet[0], l);
            const double d1l = squaredDistance(facet[1], l);
            const double d2l = squaredDistance(facet[2], l);
            if(perturbedDiameter(d01, d0l, d1l) < diam
               && perturbedDiameter(d02, d0l, d2l) < diam
               && perturbedDiameter(d12, d1l, d2l) < diam)
              is_urquhart = false;
          }

          if(is_urquhart) { // UH facet
            std::sort(facet.begin(), facet.end());
            hyperUrquhart.push_back({facet, f.first->info(), f_m.first->info(),
                                     sqrt(diam.first), diam.second});
          } else {
            const int poly1 = UF.find(f.first->info());
            const int poly2 = UF.find(f_m.first->info());
            maxDelaunay[UF.mergeRet(poly1, poly2)]
              = max({facet, sqrt(diam.first)},
                    max(maxDelaunay[poly1], maxDelaunay[poly2]));
          }

          // detect infinite polyhedrons (going beyond convex hull)
          if(del.is_infinite(f.first))
            maxDelaunay[UF.find(f.first->info())].d = inf;
          else if(del.is_infinite(f_m.first))
            maxDelaunay[UF.find(f_m.first->info())].d = inf;
        }
      }

      std::sort(hyperUrquhart.begin(), hyperUrquhart.end(),
                [](const FiltratedQuadFacet &f1, const FiltratedQuadFacet &f2) {
                  if(f1.d == f2.d)
                    return f1.a > f2.a;
                  return f1.d > f2.d;
                });
    }

    void compute2PH(UnionFind &UF, MultidimensionalDiagram &ph) {
      std::vector<int> latest(maxDelaunay.size());
      std::iota(latest.begin(), latest.end(), 0);

      for(FiltratedQuadFacet const &f :
          hyperUrquhart) { // sorted by decreasing order
        const int v1 = UF.find(f.c1);
        const int v2 = UF.find(f.c2);
        if(v1
           != v2) { // two distinct cavities: merge them by deleting the facet
          UF.merge(v1, v2);

          const int latest1 = latest[v1];
          const int latest2 = latest[v2];
          const FiltratedFacet &death1 = maxDelaunay[latest1];
          const FiltratedFacet &death2 = maxDelaunay[latest2];

          if(death1.d < death2.d) {
            if(f.d < death1.d)
              ph[2].emplace_back(
                FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                FiltratedSimplex{
                  {death1.f[0], death1.f[1], death1.f[2]}, death1.d});
            latest[UF.find(v1)] = latest2;
          } else if(death2.d < death1.d) {
            if(f.d < death2.d)
              ph[2].emplace_back(
                FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                FiltratedSimplex{
                  {death2.f[0], death2.f[1], death2.f[2]}, death2.d});
            latest[UF.find(v1)] = latest1;
          }
        } else // this is a facet from the minimal spanning acycle
          msa.push_back({f.f, f.d});
      }
    }

    void compute2PH(UnionFind &UF,
                    MultidimensionalDiagram &ph,
                    std::vector<Generator2> &generators2) {
      std::vector<int> latest(maxDelaunay.size());
      std::iota(latest.begin(), latest.end(), 0);

      std::vector<std::vector<Facet>> elementary_generators(N_c);
      for(const FiltratedQuadFacet &f : hyperUrquhart) {
        if(UF.find(f.c1) != UF.find(f.c2)) {
          if(maxDelaunay[UF.find(f.c1)].d != inf)
            elementary_generators[UF.find(f.c1)].push_back(f.f);
          if(maxDelaunay[UF.find(f.c2)].d != inf)
            elementary_generators[UF.find(f.c2)].push_back(f.f);
        }
      }
      std::vector<std::vector<int>> cascade(N_c, {0});
      for(unsigned i = 0; i < N_c; i++)
        cascade[i][0] = i;

      for(FiltratedQuadFacet const &f :
          hyperUrquhart) { // sorted by decreasing order
        const int v1 = UF.find(f.c1);
        const int v2 = UF.find(f.c2);
        if(v1
           != v2) { // two distinct cavities: merge them by deleting the facet
          UF.merge(v1, v2);

          const int latest1 = latest[v1];
          const int latest2 = latest[v2];
          const FiltratedFacet &death1 = maxDelaunay[latest1];
          const FiltratedFacet &death2 = maxDelaunay[latest2];

          if(death1.d < death2.d) {
            if(f.d < death1.d) {
              ph[2].emplace_back(
                FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                FiltratedSimplex{
                  {death1.f[0], death1.f[1], death1.f[2]}, death1.d});
              HashMap<Facet, unsigned> generator;
              for(auto const &c : cascade[latest1]) {
                for(Facet const &f_ : elementary_generators[c])
                  generator[f_]++;
              }
              std::vector<Facet> generator_facets;
              for(auto const &[f_, v] : generator) {
                if(v % 2 == 1)
                  generator_facets.push_back(f_);
              }
              generators2.push_back({generator_facets, {f.d, death1.d}});
            }
            latest[UF.find(v1)] = latest2;
            cascade[latest2].insert(cascade[latest2].end(),
                                    cascade[latest1].begin(),
                                    cascade[latest1].end());
          } else if(death2.d < death1.d) {
            if(f.d < death2.d) {
              ph[2].emplace_back(
                FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                FiltratedSimplex{
                  {death2.f[0], death2.f[1], death2.f[2]}, death2.d});
              HashMap<Facet, unsigned> generator;
              for(auto const &c : cascade[latest2]) {
                for(Facet const &f_ : elementary_generators[c])
                  generator[f_]++;
              }
              std::vector<Facet> generator_facets;
              for(auto const &[f_, v] : generator) {
                if(v % 2 == 1)
                  generator_facets.push_back(f_);
              }
              generators2.push_back({generator_facets, {f.d, death2.d}});
            }
            latest[UF.find(v1)] = latest1;
            cascade[latest1].insert(cascade[latest1].end(),
                                    cascade[latest2].begin(),
                                    cascade[latest2].end());
          }
        } else // this is a facet from the minimal spanning acycle
          msa.push_back({f.f, f.d});
      }
    }

    void connectivity() {
      // edges to facets adjacency
      msa_edges.reserve(
        1.15 * N_c); // this is an estimation of the number of edges in Delaunay
      for(unsigned i = 0; i < msa.size(); ++i) {
        const FiltratedFacet f = msa[i];
        const double d1 = squaredDistance(f.f[0], f.f[1]);
        const double d2 = squaredDistance(f.f[0], f.f[2]);
        const double d3 = squaredDistance(f.f[1], f.f[2]);
        auto &edge1 = msa_edges[std::make_pair(f.f[0], f.f[1])];
        edge1.first.reserve(4);
        edge1.first.push_back(i);
        edge1.second |= d1 > d2 && d1 > d3;
        auto &edge2 = msa_edges[std::make_pair(f.f[0], f.f[2])];
        edge2.first.reserve(4);
        edge2.first.push_back(i);
        edge2.second |= d2 > d1 && d2 > d3;
        auto &edge3 = msa_edges[std::make_pair(f.f[1], f.f[2])];
        edge3.first.reserve(4);
        edge3.first.push_back(i);
        edge3.second |= d3 > d1 && d3 > d2;
      }
    }

    void computePolygons(MultidimensionalDiagram &ph) {
      const unsigned N_msa = msa.size();
      UF_msa = UnionFind(N_msa);

      for(const auto &[e, val] : msa_edges) {
        auto &[adj_f, isNotUrquhart] = val;
        if(isNotUrquhart) { // not UG edge
          if(adj_f.size() == 1)
            msa[UF_msa.find(adj_f[0])].d = inf;
          else if(adj_f.size() == 2)
            UF_msa.merge(adj_f[0], adj_f[1]);
          else // non-manifold junction edge
            critical1.push_back({e, sqrt(squaredDistance(e.first, e.second))});
        } else
          urquhart.push_back({e, sqrt(squaredDistance(e.first, e.second))});
      }

      for(unsigned x = 0; x < N_msa; ++x) {
        const int poly = UF_msa.find(x);
        msa[poly] = max(msa[poly], msa[x]);
      }

      std::sort(urquhart.begin(), urquhart.end(),
                [](const FiltratedEdge &a, const FiltratedEdge &b) {
                  return a.d < b.d;
                });
      UnionFind UF_p(N_p);
      ph[0].reserve(N_p - 1);
      for(FiltratedEdge const &e : urquhart) {
        if(UF_p.find(e.e.first)
           != UF_p.find(e.e.second)) { // we know e is a EMST edge
          UF_p.merge(e.e.first, e.e.second);
          ph[0].emplace_back(FiltratedSimplex{{-1}, 0.},
                             FiltratedSimplex{{e.e.first, e.e.second}, e.d});
        } else
          critical1.emplace_back(e);
      }
    }

    void compute1PH(MultidimensionalDiagram &ph) {
      const unsigned N_msa = msa.size();
      std::vector<int> polys;
      for(unsigned x = 0; x < N_msa; ++x) {
        if(UF_msa.isRoot(x) && msa[x].d < inf)
          polys.push_back(x);
      }

      std::sort(polys.begin(), polys.end(), [&](const int x1, const int x2) {
        return msa[x1].d < msa[x2].d;
      });

      std::vector<int> criticalIndices(critical1.size());
      std::iota(criticalIndices.begin(), criticalIndices.end(), 0);
      std::sort(criticalIndices.begin(), criticalIndices.end(),
                [&](const int e1_id, const int e2_id) {
                  return critical1[e1_id].d < critical1[e2_id].d;
                });
      std::vector<int> criticalOrder(critical1.size());
      for(unsigned i = 0; i < criticalIndices.size(); ++i)
        criticalOrder[criticalIndices[i]] = i;

      std::vector<std::vector<int>> poly_to_crit(N_msa);
      for(const int poly : polys)
        poly_to_crit[poly].reserve(3);
      for(unsigned i = 0; i < critical1.size(); ++i) {
        const FiltratedEdge &e = critical1[i];
        for(const int poly : msa_edges[e.e].first) {
          if(msa[UF_msa.find(poly)].d < inf) {
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

      std::vector<int> partner(critical1.size(), -1);
      for(const int poly : polys) {
        HashSet<int> boundary(poly_to_crit[poly].begin(),
                              poly_to_crit[poly].end(),
                              poly_to_crit[poly].size());
        while(true) {
          const int youngest_id
            = *std::max_element(boundary.begin(), boundary.end());
          if(partner[youngest_id] == -1) {
            partner[youngest_id] = poly;
            const FiltratedEdge &e = critical1[criticalIndices[youngest_id]];
            const FiltratedFacet &death = msa[poly];
            if(e.d < death.d)
              ph[1].emplace_back(
                FiltratedSimplex{{e.e.first, e.e.second}, e.d},
                FiltratedSimplex{
                  {death.f[0], death.f[1], death.f[2]}, death.d});
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
    }

    void compute1PH(MultidimensionalDiagram &ph,
                    std::vector<Generator1> &generators1) {
      const unsigned N_msa = msa.size();
      std::vector<int> polys;
      for(unsigned x = 0; x < N_msa; ++x) {
        if(UF_msa.isRoot(x) && msa[x].d < inf)
          polys.push_back(x);
      }

      std::sort(polys.begin(), polys.end(), [&](const int x1, const int x2) {
        return msa[x1].d < msa[x2].d;
      });

      std::vector<int> criticalIndices(critical1.size());
      std::iota(criticalIndices.begin(), criticalIndices.end(), 0);
      std::sort(criticalIndices.begin(), criticalIndices.end(),
                [&](const int e1_id, const int e2_id) {
                  return critical1[e1_id].d < critical1[e2_id].d;
                });
      std::vector<int> criticalOrder(critical1.size());
      for(unsigned i = 0; i < criticalIndices.size(); ++i)
        criticalOrder[criticalIndices[i]] = i;

      std::vector<std::vector<int>> poly_to_crit(N_msa);
      for(const int poly : polys)
        poly_to_crit[poly].reserve(3);
      for(unsigned i = 0; i < critical1.size(); ++i) {
        const FiltratedEdge &e = critical1[i];
        for(const int poly : msa_edges[e.e].first) {
          if(msa[UF_msa.find(poly)].d < inf) {
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

      std::vector<std::vector<Edge>> elementary_generators(N_msa);
      auto action = [&](const FiltratedEdge &e) {
        auto adj_f = msa_edges[e.e].first;
        for(int &f : adj_f)
          f = UF_msa.find(f);
        std::sort(adj_f.begin(), adj_f.end());
        const auto last = std::unique(adj_f.begin(), adj_f.end());
        for(auto it = adj_f.begin(); it != last; ++it) {
          if(msa[UF_msa.find(*it)].d != inf)
            elementary_generators[UF_msa.find(*it)].push_back(e.e);
        }
      };
      for(auto const &e : urquhart)
        action(e);
      for(auto const &e : critical1) {
        if(msa_edges[e.e]
             .second) // if critical but not urquhart -> non manifold junctions
          action(e);
      }

      std::vector<int> partner(critical1.size(), -1);
      for(const int poly : polys) {
        HashSet<int> boundary(poly_to_crit[poly].begin(),
                              poly_to_crit[poly].end(),
                              poly_to_crit[poly].size());
        std::vector<Edge> &generator = elementary_generators[poly];
        while(true) {
          const int youngest_id
            = *std::max_element(boundary.begin(), boundary.end());
          if(partner[youngest_id] == -1) {
            partner[youngest_id] = poly;
            const FiltratedEdge &e = critical1[criticalIndices[youngest_id]];
            const FiltratedFacet &death = msa[poly];
            if(e.d < death.d) {
              ph[1].emplace_back(
                FiltratedSimplex{{e.e.first, e.e.second}, e.d},
                FiltratedSimplex{
                  {death.f[0], death.f[1], death.f[2]}, death.d});
              generators1.emplace_back(generator, std::make_pair(e.d, death.d));
            }
            break;
          } else {
            for(const int crit : poly_to_crit[partner[youngest_id]]) {
              const auto it = boundary.find(crit);
              if(it == boundary.end())
                boundary.insert(crit);
              else
                boundary.erase(it);
            }
            for(Edge const &e : elementary_generators[partner[youngest_id]]) {
              auto it = std::find(generator.begin(), generator.end(), e);
              if(it == generator.end())
                generator.push_back(e);
              else
                generator.erase(it);
            }
          }
        }
        poly_to_crit[poly].assign(boundary.begin(), boundary.end());
      }
    }
  };

#ifdef TTK_GPH_PARALLEL
  class DRPersistence3_p {
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Vb = CGAL::Triangulation_vertex_base_with_info_3<int, K>;
    using Fb = CGAL::Triangulation_cell_base_with_info_3<int, K>;
    using Tds = CGAL::
      Triangulation_data_structure_3<Vb, Fb, CGAL::Parallel_if_available_tag>;
    using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;
    using Point = Delaunay::Point_3;

    using AdjacencyList = std::vector<int>;
    using ConnectivityHashMap = HashMap<Edge, std::pair<AdjacencyList, bool>>;
    using ConcurrentConnectivityHashMap
      = ConcurrentHashMap<Edge, std::pair<AdjacencyList, bool>>;

  public:
    explicit DRPersistence3_p(const PointCloud<3> &points,
                              const unsigned nThreads = 1)
      : N_p(points.size()), p(points), nThreads_(nThreads) {
    }

    /**
     * Computes the Delaunay-Rips persistence diagram of the point cloud given
     * to the constructor
     */
    void run(MultidimensionalDiagram &ph) {
      ph = MultidimensionalDiagram(3);

      computeDelaunay();

      UnionFind UF(N_c);
      computeFacetUrquhart(UF);
      compute2PH(UF, ph);

      connectivity();
      computePolygons(ph);
      compute1PH(ph);
    }

    /**
     * Computes the Delaunay-Rips persistence diagram and associated 1- and
     * 2-dimensional generators of the point cloud given to the constructor
     */
    void run(MultidimensionalDiagram &ph,
             std::vector<Generator1> &generators1,
             std::vector<Generator2> &generators2) {
      ph = MultidimensionalDiagram(3);
      generators1.resize(0);
      generators2.resize(0);

      computeDelaunay();

      UnionFind UF(N_c);
      computeFacetUrquhart(UF);
      compute2PH(UF, ph, generators2);

      connectivity();
      computePolygons(ph);
      compute1PH(ph, generators1);
    }

  private:
    const unsigned N_p;
    unsigned N_c{0};
    const PointCloud<3> &p;
    Delaunay del;
    const int nThreads_{1};

    // 1-dimensional
    tbb::concurrent_vector<FiltratedEdge> urquhart;
    tbb::concurrent_vector<FiltratedEdge> critical1;

    ConcurrentConnectivityHashMap msa_edges{};

    DisjointSets UF_msa{0};

    // 2-dimensional
    std::vector<FiltratedQuadFacet> hyperUrquhart;
    std::vector<FiltratedFacet> msa;
    std::vector<FiltratedFacet> maxDelaunay;

    [[nodiscard]] double squaredDistance(const unsigned i1,
                                         const unsigned i2) const {
      PointD<3> const p1 = p[i1], p2 = p[i2];
      return (p1[0] - p2[0]) * (p1[0] - p2[0])
             + (p1[1] - p2[1]) * (p1[1] - p2[1])
             + (p1[2] - p2[2]) * (p1[2] - p2[2]);
    }
    [[nodiscard]] static std::pair<double, double>
      perturbedDiameter(const double d1, const double d2, const double d3) {
      return {std::max(d1, std::max(d2, d3)), d1 + d2 + d3};
    }

    /**
     * compute Delaunay and initialize tetrahedra IDs (including infinite
     * tetrahedra to deal with the convex hull)
     */
    void computeDelaunay() {
      tbb::global_control gc(
        tbb::global_control::max_allowed_parallelism, nThreads_);
      std::vector<std::pair<Point, unsigned>> points(N_p);
      CGAL::Bbox_3 bbox(p[0][0], p[0][1], p[0][2], p[0][0], p[0][1], p[0][2]);
      for(unsigned i = 0; i < N_p; ++i) {
        points[i] = std::make_pair(Point(p[i][0], p[i][1], p[i][2]), i);
        bbox
          += CGAL::Bbox_3(p[i][0], p[i][1], p[i][2], p[i][0], p[i][1], p[i][2]);
      }
      Delaunay::Lock_data_structure lock_ds(bbox, 100);
      del = Delaunay(points.begin(), points.end(), &lock_ds);
      int k = 0;
      for(auto const &c : del.all_cell_handles())
        c->info() = k++;
      N_c = k;
    }

    void computeFacetUrquhart(UnionFind &UF) {
      maxDelaunay.resize(N_c, FiltratedFacet{{-1, -1, -1}, 0.});
      for(Delaunay::Facet const &f : del.all_facets()) {
        if(del.is_infinite(f))
          UF.merge(f.first->info(), del.mirror_facet(f).first->info());
        else {
          Facet facet{f.first->vertex((f.second + 1) % 4)->info(),
                      f.first->vertex((f.second + 2) % 4)->info(),
                      f.first->vertex((f.second + 3) % 4)->info()};
          const double d01 = squaredDistance(facet[0], facet[1]);
          const double d02 = squaredDistance(facet[0], facet[2]);
          const double d12 = squaredDistance(facet[1], facet[2]);
          const auto diam = perturbedDiameter(d01, d02, d12);
          const Delaunay::Facet f_m = del.mirror_facet(f);

          // check if facet is Urquhart
          bool is_urquhart = true;
          if(!del.is_infinite(f.first->vertex(f.second))) {
            const unsigned k = f.first->vertex(f.second)->info();
            const double d0k = squaredDistance(facet[0], k);
            const double d1k = squaredDistance(facet[1], k);
            const double d2k = squaredDistance(facet[2], k);
            if(perturbedDiameter(d01, d0k, d1k) < diam
               && perturbedDiameter(d02, d0k, d2k) < diam
               && perturbedDiameter(d12, d1k, d2k) < diam)
              is_urquhart = false;
          }
          if(is_urquhart && !del.is_infinite(f_m.first->vertex(f_m.second))) {
            const unsigned l = f_m.first->vertex(f_m.second)->info();
            const double d0l = squaredDistance(facet[0], l);
            const double d1l = squaredDistance(facet[1], l);
            const double d2l = squaredDistance(facet[2], l);
            if(perturbedDiameter(d01, d0l, d1l) < diam
               && perturbedDiameter(d02, d0l, d2l) < diam
               && perturbedDiameter(d12, d1l, d2l) < diam)
              is_urquhart = false;
          }

          if(is_urquhart) { // UH facet
            std::sort(facet.begin(), facet.end());
            hyperUrquhart.push_back({facet, f.first->info(), f_m.first->info(),
                                     sqrt(diam.first), diam.second});
          } else {
            const int poly1 = UF.find(f.first->info());
            const int poly2 = UF.find(f_m.first->info());
            maxDelaunay[UF.mergeRet(poly1, poly2)]
              = max({facet, sqrt(diam.first)},
                    max(maxDelaunay[poly1], maxDelaunay[poly2]));
          }

          // detect infinite polyhedrons (going beyond convex hull)
          if(del.is_infinite(f.first))
            maxDelaunay[UF.find(f.first->info())].d = inf;
          else if(del.is_infinite(f_m.first))
            maxDelaunay[UF.find(f_m.first->info())].d = inf;
        }
      }

      TTK_PSORT(nThreads_, hyperUrquhart.begin(), hyperUrquhart.end(),
                [](const FiltratedQuadFacet &f1, const FiltratedQuadFacet &f2) {
                  if(f1.d == f2.d)
                    return f1.a > f2.a;
                  return f1.d > f2.d;
                });
    }

    void compute2PH(UnionFind &UF, MultidimensionalDiagram &ph) {
      std::vector<int> latest(maxDelaunay.size());
      std::iota(latest.begin(), latest.end(), 0);

      for(FiltratedQuadFacet const &f :
          hyperUrquhart) { // sorted by decreasing order
        const int v1 = UF.find(f.c1);
        const int v2 = UF.find(f.c2);
        if(v1
           != v2) { // two distinct cavities: merge them by deleting the facet
          UF.merge(v1, v2);

          const int latest1 = latest[v1];
          const int latest2 = latest[v2];
          const FiltratedFacet &death1 = maxDelaunay[latest1];
          const FiltratedFacet &death2 = maxDelaunay[latest2];

          if(death1.d < death2.d) {
            if(f.d < death1.d)
              ph[2].emplace_back(
                FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                FiltratedSimplex{
                  {death1.f[0], death1.f[1], death1.f[2]}, death1.d});
            latest[UF.find(v1)] = latest2;
          } else if(death2.d < death1.d) {
            if(f.d < death2.d)
              ph[2].emplace_back(
                FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                FiltratedSimplex{
                  {death2.f[0], death2.f[1], death2.f[2]}, death2.d});
            latest[UF.find(v1)] = latest1;
          }
        } else // this is a facet from the minimal spanning acycle
          msa.push_back({f.f, f.d});
      }
    }

    void compute2PH(UnionFind &UF,
                    MultidimensionalDiagram &ph,
                    std::vector<Generator2> &generators2) {
      std::vector<int> latest(maxDelaunay.size());
      std::iota(latest.begin(), latest.end(), 0);

      std::vector<std::vector<Facet>> elementary_generators(N_c);
      for(const FiltratedQuadFacet &f : hyperUrquhart) {
        if(UF.find(f.c1) != UF.find(f.c2)) {
          if(maxDelaunay[UF.find(f.c1)].d != inf)
            elementary_generators[UF.find(f.c1)].push_back(f.f);
          if(maxDelaunay[UF.find(f.c2)].d != inf)
            elementary_generators[UF.find(f.c2)].push_back(f.f);
        }
      }
      std::vector<std::vector<int>> cascade(N_c, {0});
      for(unsigned i = 0; i < N_c; i++)
        cascade[i][0] = i;

      for(FiltratedQuadFacet const &f :
          hyperUrquhart) { // sorted by decreasing order
        const int v1 = UF.find(f.c1);
        const int v2 = UF.find(f.c2);
        if(v1
           != v2) { // two distinct cavities: merge them by deleting the facet
          UF.merge(v1, v2);

          const int latest1 = latest[v1];
          const int latest2 = latest[v2];
          const FiltratedFacet &death1 = maxDelaunay[latest1];
          const FiltratedFacet &death2 = maxDelaunay[latest2];

          if(death1.d < death2.d) {
            if(f.d < death1.d) {
              ph[2].emplace_back(
                FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                FiltratedSimplex{
                  {death1.f[0], death1.f[1], death1.f[2]}, death1.d});
              HashMap<Facet, unsigned> generator;
              for(auto const &c : cascade[latest1]) {
                for(Facet const &f_ : elementary_generators[c])
                  generator[f_]++;
              }
              std::vector<Facet> generator_facets;
              for(auto const &[f_, v] : generator) {
                if(v % 2 == 1)
                  generator_facets.push_back(f_);
              }
              generators2.push_back({generator_facets, {f.d, death1.d}});
            }
            latest[UF.find(v1)] = latest2;
            cascade[latest2].insert(cascade[latest2].end(),
                                    cascade[latest1].begin(),
                                    cascade[latest1].end());
          } else if(death2.d < death1.d) {
            if(f.d < death2.d) {
              ph[2].emplace_back(
                FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                FiltratedSimplex{
                  {death2.f[0], death2.f[1], death2.f[2]}, death2.d});
              HashMap<Facet, unsigned> generator;
              for(auto const &c : cascade[latest2]) {
                for(Facet const &f_ : elementary_generators[c])
                  generator[f_]++;
              }
              std::vector<Facet> generator_facets;
              for(auto const &[f_, v] : generator) {
                if(v % 2 == 1)
                  generator_facets.push_back(f_);
              }
              generators2.push_back({generator_facets, {f.d, death2.d}});
            }
            latest[UF.find(v1)] = latest1;
            cascade[latest1].insert(cascade[latest1].end(),
                                    cascade[latest2].begin(),
                                    cascade[latest2].end());
          }
        } else // this is a facet from the minimal spanning acycle
          msa.push_back({f.f, f.d});
      }
    }

    void connectivity() {
      omp_set_num_threads(nThreads_);
      // edges to facets adjacency
      msa_edges.reserve(
        1.15 * N_c); // this is an estimation of the number of edges in Delaunay
#pragma omp parallel for
      for(unsigned i = 0; i < msa.size(); ++i) {
        const FiltratedFacet f = msa[i];
        const double d1 = squaredDistance(f.f[0], f.f[1]);
        const double d2 = squaredDistance(f.f[0], f.f[2]);
        const double d3 = squaredDistance(f.f[1], f.f[2]);
        msa_edges.emplace_or_visit(
          std::make_pair(f.f[0], f.f[1]),
          std::make_pair(AdjacencyList{(int)i}, d1 > d2 && d1 > d3),
          [&](auto &x) {
            x.second.first.push_back(i);
            x.second.second |= d1 > d2 && d1 > d3;
          });
        msa_edges.emplace_or_visit(
          std::make_pair(f.f[0], f.f[2]),
          std::make_pair(AdjacencyList{(int)i}, d2 > d1 && d2 > d3),
          [&](auto &x) {
            x.second.first.push_back(i);
            x.second.second |= d2 > d1 && d2 > d3;
          });
        msa_edges.emplace_or_visit(
          std::make_pair(f.f[1], f.f[2]),
          std::make_pair(AdjacencyList{(int)i}, d3 > d1 && d3 > d2),
          [&](auto &x) {
            x.second.first.push_back(i);
            x.second.second |= d3 > d1 && d3 > d2;
          });
      }
    }

    void computePolygons(MultidimensionalDiagram &ph) {
      omp_set_num_threads(nThreads_);
      const unsigned N_msa = msa.size();
      UF_msa = DisjointSets(N_msa);

      msa_edges.cvisit_all(
#ifdef __cpp_lib_execution
        std::execution::par,
#endif
        [&](const auto &x) {
          auto &[e, val] = x;
          auto &[adj_f, isNotUrquhart] = val;
          if(isNotUrquhart) { // not UG edge
            if(adj_f.size() == 1)
              msa[adj_f[0]].d = inf;
            if(adj_f.size() == 2)
              UF_msa.unite(adj_f[0], adj_f[1]);
            else // non-manifold junction edge
              critical1.push_back(
                {e, sqrt(squaredDistance(e.first, e.second))});
          } else
            urquhart.push_back({e, sqrt(squaredDistance(e.first, e.second))});
        });

      std::vector<std::mutex> maxDelaunayLocks(N_msa);
#pragma omp parallel for
      for(unsigned x = 0; x < N_msa; ++x) {
        const int poly = UF_msa.find(x);
        std::lock_guard lock(maxDelaunayLocks[poly]);
        msa[poly] = max(msa[poly], msa[x]);
      }

      TTK_PSORT(nThreads_, urquhart.begin(), urquhart.end(),
                [](const FiltratedEdge &a, const FiltratedEdge &b) {
                  return a.d < b.d;
                });
      UnionFind UF_p(N_p);
      ph[0].reserve(N_p - 1);
      for(FiltratedEdge const &e : urquhart) {
        if(UF_p.find(e.e.first)
           != UF_p.find(e.e.second)) { // we know e is a EMST edge
          UF_p.merge(e.e.first, e.e.second);
          ph[0].emplace_back(FiltratedSimplex{{-1}, 0.},
                             FiltratedSimplex{{e.e.first, e.e.second}, e.d});
        } else
          critical1.emplace_back(e);
      }
    }

    void compute1PH(MultidimensionalDiagram &ph) {
      omp_set_num_threads(nThreads_);
      const ConnectivityHashMap sequential_msa_edges = std::move(msa_edges);

      const unsigned N_msa = msa.size();
      std::vector<int> polys;
      for(unsigned x = 0; x < N_msa; ++x) {
        if(UF_msa.isRoot(x) && msa[x].d < inf)
          polys.push_back(x);
      }

      TTK_PSORT(
        nThreads_, polys.begin(), polys.end(),
        [&](const int x1, const int x2) { return msa[x1].d < msa[x2].d; });

      std::vector<int> criticalIndices(critical1.size());
      std::iota(criticalIndices.begin(), criticalIndices.end(), 0);
      TTK_PSORT(nThreads_, criticalIndices.begin(), criticalIndices.end(),
                [&](const int e1_id, const int e2_id) {
                  return critical1[e1_id].d < critical1[e2_id].d;
                });
      std::vector<int> criticalOrder(critical1.size());
      for(unsigned i = 0; i < criticalIndices.size(); ++i)
        criticalOrder[criticalIndices[i]] = i;

      std::vector<std::vector<int>> poly_to_crit(N_msa);
      for(const int poly : polys)
        poly_to_crit[poly].reserve(3);
      std::vector<std::mutex> poly_mutex(N_msa);
      std::vector<std::mutex> crit_mutex(critical1.size());

#pragma omp parallel for
      for(unsigned i = 0; i < critical1.size(); ++i) {
        const FiltratedEdge &e = critical1[i];
        for(const int poly : sequential_msa_edges.find(e.e)->second.first) {
          const int uf_poly = UF_msa.find(poly);
          if(msa[uf_poly].d < inf) {
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

      std::vector<int> partner(critical1.size(), -1);

      auto eliminateBoundary = [&](const int poly) -> int {
        std::lock_guard lock(poly_mutex[poly]);
        HashSet<int> boundary(poly_to_crit[poly].begin(),
                              poly_to_crit[poly].end(),
                              poly_to_crit[poly].size());
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
          if(msa[partner[youngest_id]].d > msa[poly].d) {
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
      for(const int poly : polys) {
        int tmp_p = poly;
        while(tmp_p >= 0)
          tmp_p = eliminateBoundary(tmp_p);
      }

      for(unsigned i = 0; i < critical1.size(); ++i) {
        const int poly = partner[i];
        const FiltratedEdge &e = critical1[criticalIndices[i]];
        const FiltratedFacet &death = msa[poly];
        if(e.d < death.d)
          ph[1].emplace_back(
            FiltratedSimplex{{e.e.first, e.e.second}, e.d},
            FiltratedSimplex{{death.f[0], death.f[1], death.f[2]}, death.d});
      }
    }

    void compute1PH(MultidimensionalDiagram &ph,
                    std::vector<Generator1> &generators1) {
      omp_set_num_threads(nThreads_);
      ConnectivityHashMap sequential_msa_edges = std::move(msa_edges);

      const unsigned N_msa = msa.size();
      std::vector<int> polys;
      for(unsigned x = 0; x < N_msa; ++x) {
        if(UF_msa.isRoot(x) && msa[x].d < inf)
          polys.push_back(x);
      }

      TTK_PSORT(
        nThreads_, polys.begin(), polys.end(),
        [&](const int x1, const int x2) { return msa[x1].d < msa[x2].d; });

      std::vector<int> criticalIndices(critical1.size());
      std::iota(criticalIndices.begin(), criticalIndices.end(), 0);
      TTK_PSORT(nThreads_, criticalIndices.begin(), criticalIndices.end(),
                [&](const int e1_id, const int e2_id) {
                  return critical1[e1_id].d < critical1[e2_id].d;
                });
      std::vector<int> criticalOrder(critical1.size());
      for(unsigned i = 0; i < criticalIndices.size(); ++i)
        criticalOrder[criticalIndices[i]] = i;

      std::vector<std::vector<int>> poly_to_crit(N_msa);
      for(const int poly : polys)
        poly_to_crit[poly].reserve(3);
      std::vector<std::mutex> poly_mutex(N_msa);
      std::vector<std::mutex> crit_mutex(critical1.size());

#pragma omp parallel for
      for(unsigned i = 0; i < critical1.size(); ++i) {
        const FiltratedEdge &e = critical1[i];
        for(const int poly : sequential_msa_edges.find(e.e)->second.first) {
          const int uf_poly = UF_msa.find(poly);
          if(msa[uf_poly].d < inf) {
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

      std::vector<std::vector<Edge>> elementary_generators(N_msa);
      auto action = [&](const FiltratedEdge &e) {
        auto adj_f = sequential_msa_edges[e.e].first;
        for(int &f : adj_f)
          f = UF_msa.find(f);
        std::sort(adj_f.begin(), adj_f.end());
        const auto last = std::unique(adj_f.begin(), adj_f.end());
        for(auto it = adj_f.begin(); it != last; ++it) {
          if(msa[UF_msa.find(*it)].d != inf)
            elementary_generators[UF_msa.find(*it)].push_back(e.e);
        }
      };
      for(auto const &e : urquhart)
        action(e);
      for(auto const &e : critical1) {
        if(sequential_msa_edges[e.e]
             .second) // if critical but not urquhart -> non manifold junctions
          action(e);
      }

      std::vector<int> partner(critical1.size(), -1);

      auto eliminateBoundary = [&](const int poly) -> int {
        std::lock_guard lock(poly_mutex[poly]);
        HashSet<int> boundary(poly_to_crit[poly].begin(),
                              poly_to_crit[poly].end(),
                              poly_to_crit[poly].size());
        std::vector<Edge> &generator = elementary_generators[poly];
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
          if(msa[partner[youngest_id]].d > msa[poly].d) {
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
          for(Edge const &e : elementary_generators[partner[youngest_id]]) {
            auto it = std::find(generator.begin(), generator.end(), e);
            if(it == generator.end())
              generator.push_back(e);
            else
              generator.erase(it);
          }
        }
      };

#pragma omp parallel for
      for(const int poly : polys) {
        int tmp_p = poly;
        while(tmp_p >= 0)
          tmp_p = eliminateBoundary(tmp_p);
      }

      for(unsigned i = 0; i < critical1.size(); ++i) {
        const int poly = partner[i];
        const FiltratedEdge &e = critical1[criticalIndices[i]];
        const FiltratedFacet &death = msa[poly];
        const std::vector<Edge> &generator = elementary_generators[poly];
        if(e.d < death.d) {
          ph[1].emplace_back(
            FiltratedSimplex{{e.e.first, e.e.second}, e.d},
            FiltratedSimplex{{death.f[0], death.f[1], death.f[2]}, death.d});
          generators1.emplace_back(generator, std::make_pair(e.d, death.d));
        }
      }
    }
  };
#endif

  inline void
    runDelaunayRipsPersistenceDiagram3(rpd::PointCloud const &points,
                                       MultidimensionalDiagram &diagram,
                                       int threads = 1) {
    PointCloud<3> p(points.size());
    for(unsigned i = 0; i < points.size(); ++i)
      p[i] = {points[i][0], points[i][1], points[i][2]};
    if(threads > 1) {
#ifdef TTK_GPH_PARALLEL
      DRPersistence3_p drpd(p, threads);
#else
      DRPersistence3 drpd(p);
#endif
      drpd.run(diagram);
    } else {
      DRPersistence3 drpd(p);
      drpd.run(diagram);
    }
  }

  inline void
    runDelaunayRipsPersistenceDiagram3(rpd::PointCloud const &points,
                                       MultidimensionalDiagram &diagram,
                                       std::vector<Generator1> &generators1,
                                       std::vector<Generator2> &generators2,
                                       int threads = 1) {
    PointCloud<3> p(points.size());
    for(unsigned i = 0; i < points.size(); ++i)
      p[i] = {points[i][0], points[i][1], points[i][2]};
    if(threads > 1) {
#ifdef TTK_GPH_PARALLEL
      DRPersistence3_p drpd(p, threads);
#else
      DRPersistence3 drpd(p);
#endif
      drpd.run(diagram, generators1, generators2);
    } else {
      DRPersistence3 drpd(p);
      drpd.run(diagram, generators1, generators2);
    }
  }

} // namespace ttk::gph

#endif