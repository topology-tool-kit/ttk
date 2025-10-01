#pragma once

#include "geoPHUtils.h"

#ifdef TTK_ENABLE_CGAL

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/global_control.h>
#endif

namespace ttk::gph {

  class DRPersistence3 {
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Vb = CGAL::Triangulation_vertex_base_with_info_3<int, K>;
    using Fb = CGAL::Triangulation_cell_base_with_info_3<int, K>;
    using Tds = CGAL::Triangulation_data_structure_3<Vb, Fb, CGAL::Parallel_if_available_tag>;
    using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;
    using Point = Delaunay::Point_3;

    using AdjacencyList = std::vector<int>;
    using ConnectivityHashMap = GPH_HASHMAP<Edge, std::pair<AdjacencyList, bool>, boost::hash<Edge>>;

  public:
    explicit DRPersistence3(const PointCloud<3> &points, const int nThreadsSort = 1, const int nThreadsDelaunay = 1) : N_p(points.size()), p(points), nThreadsSort_(nThreadsSort), nThreadsDelaunay_(nThreadsDelaunay) {}

    /**
     * Computes the Delaunay-Rips persistence diagram of the point cloud given to the constructor
     */
    void computeDelaunayRipsPersistence(MultidimensionalDiagram &ph) {
      ph = MultidimensionalDiagram(3);

      computeDelaunay();

      UnionFind UF (N_c);
      computeFacetUrquhart(UF);
      compute2PH(UF, ph);

      connectivity();
      computePolygons(ph);
      compute1PH(ph);
    }

    /**
     * Computes the Delaunay-Rips persistence diagram and associated 1- and 2-dimensional generators
     * of the point cloud given to the constructor
     */
    void computeDelaunayRipsPersistence(MultidimensionalDiagram &ph, std::vector<Generator1> &generators1, std::vector<Generator2> &generators2)  {
      ph = MultidimensionalDiagram(3);
      generators1.resize(0);
      generators2.resize(0);

      computeDelaunay();

      UnionFind UF (N_c);
      computeFacetUrquhart(UF);
      compute2PH(UF, ph, generators2);

      connectivity();
      computePolygons(ph);
      compute1PH(ph, generators1);
    }

  private:
    const unsigned N_p;
    unsigned N_c {0};
    const PointCloud<3>& p;
    Delaunay del;
    const int nThreadsSort_;
    const int nThreadsDelaunay_;

    //1-dimensional
    std::vector<FiltratedEdge> urquhart;
    std::vector<FiltratedEdge> critical1;
    std::vector<FiltratedEdge> maxDelaunay1;

    ConnectivityHashMap msa_edges {};

    UnionFind UF_msa {0}; //necessary for 1-generators, could be local otherwise

    //2-dimensional
    std::vector<FiltratedQuadFacet> hyperUrquhart;
    std::vector<FiltratedQuadFacet> msa;
    std::vector<FiltratedFacet> maxDelaunay;

    [[nodiscard]] double squaredDistance(const unsigned i1, const unsigned i2) const {
      PointD<3> const p1 = p[i1], p2 = p[i2];
      return (p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2]);
    }
    [[nodiscard]] static std::pair<double,double> perturbedDiameter(const double d1, const double d2, const double d3) {
      return {std::max(d1, std::max(d2, d3)), d1 + d2 + d3};
    }

    /**
     * compute Delaunay and initialize tetrahedra IDs (including infinite tetrahedra to deal with the convex hull)
     */
    void computeDelaunay() {
      std::vector<std::pair<Point,unsigned>> points (N_p);
#ifndef CGAL_LINKED_WITH_TBB
      for (unsigned i=0; i<N_p; ++i)
        points[i] = std::make_pair(Point(p[i][0], p[i][1], p[i][2]), i);
      del = Delaunay(points.begin(), points.end());
#else
      tbb::global_control (tbb::global_control::max_allowed_parallelism, nThreadsDelaunay_);
      if (nThreadsDelaunay_ == 1) {
        for (unsigned i=0; i<N_p; ++i)
          points[i] = std::make_pair(Point(p[i][0], p[i][1], p[i][2]), i);
        del = Delaunay(points.begin(), points.end());
      }
      else {
        CGAL::Bbox_3 bbox(p[0][0], p[0][1], p[0][2], p[0][0], p[0][1], p[0][2]);
        for (unsigned i=0; i<N_p; ++i) {
          points[i] = std::make_pair(Point(p[i][0], p[i][1], p[i][2]), i);
          bbox += CGAL::Bbox_3(p[i][0], p[i][1], p[i][2], p[i][0], p[i][1], p[i][2]);
        }
        Delaunay::Lock_data_structure lock_ds(bbox, 100);
        del = Delaunay(points.begin(), points.end(), &lock_ds);
      }
#endif
      int k = 0;
      for(auto const& c : del.all_cell_handles())
        c->info() = k++;
      N_c = k;
    }

    void computeFacetUrquhart(UnionFind &UF) {
      maxDelaunay.resize(N_c, FiltratedFacet{{-1,-1,-1}, 0.});
      for (Delaunay::Facet const& f : del.all_facets()) {
        if (del.is_infinite(f))
          UF.merge(f.first->info(), del.mirror_facet(f).first->info());
        else {
          Facet facet {f.first->vertex((f.second + 1) % 4)->info(),
                       f.first->vertex((f.second + 2) % 4)->info(),
                       f.first->vertex((f.second + 3) % 4)->info()};
          const double d01 = squaredDistance(facet[0], facet[1]);
          const double d02 = squaredDistance(facet[0], facet[2]);
          const double d12 = squaredDistance(facet[1], facet[2]);
          const auto diam = perturbedDiameter(d01, d02, d12);
          const Delaunay::Facet f_m = del.mirror_facet(f);

          // check if facet is Urquhart
          bool is_urquhart = true;
          if (!del.is_infinite(f.first->vertex(f.second))) {
            const unsigned k = f.first->vertex(f.second)->info();
            const double d0k = squaredDistance(facet[0], k);
            const double d1k = squaredDistance(facet[1], k);
            const double d2k = squaredDistance(facet[2], k);
            if (perturbedDiameter(d01, d0k, d1k) < diam
             && perturbedDiameter(d02, d0k, d2k) < diam
             && perturbedDiameter(d12, d1k, d2k) < diam)
              is_urquhart = false;
          }
          if (is_urquhart && !del.is_infinite(f_m.first->vertex(f_m.second))) {
            const unsigned l = f_m.first->vertex(f_m.second)->info();
            const double d0l = squaredDistance(facet[0], l);
            const double d1l = squaredDistance(facet[1], l);
            const double d2l = squaredDistance(facet[2], l);
            if (perturbedDiameter(d01, d0l, d1l) < diam
             && perturbedDiameter(d02, d0l, d2l) < diam
             && perturbedDiameter(d12, d1l, d2l) < diam)
              is_urquhart = false;
          }

          if (is_urquhart) { //UH facet
            std::sort(facet.begin(), facet.end());
            hyperUrquhart.push_back({facet, f.first->info(), f_m.first->info(), sqrt(diam.first), diam.second});
          }
          else {
            const int poly1 = UF.find(f.first->info());
            const int poly2 = UF.find(f_m.first->info());
            maxDelaunay[UF.mergeRet(poly1, poly2)] = max({facet, sqrt(diam.first)},
                                                             max(maxDelaunay[poly1], maxDelaunay[poly2]));
          }

          //detect infinite polyhedrons (going beyond convex hull)
          if (del.is_infinite(f.first))
            maxDelaunay[UF.find(f.first->info())].d = inf;
          else if (del.is_infinite(f_m.first))
            maxDelaunay[UF.find(f_m.first->info())].d = inf;
        }
      }

      TTK_PSORT(nThreadsSort_, hyperUrquhart.begin(), hyperUrquhart.end(), [](const FiltratedQuadFacet &f1, const FiltratedQuadFacet &f2) {
        if (f1.d == f2.d)
          return f1.a > f2.a;
        else
          return f1.d > f2.d;
      });
    }

    void compute2PH(UnionFind &UF, MultidimensionalDiagram &ph)  {
      std::vector<int> latest(maxDelaunay.size());
      std::iota(latest.begin(), latest.end(), 0);

      for (FiltratedQuadFacet const& f : hyperUrquhart) { //sorted by decreasing order
        const int v1 = UF.find(f.c1);
        const int v2 = UF.find(f.c2);
        if (v1 != v2) { // two distinct cavities: merge them by deleting the facet
          UF.merge(v1, v2);

          const int latest1 = latest[v1];
          const int latest2 = latest[v2];
          const FiltratedFacet& death1 = maxDelaunay[latest1];
          const FiltratedFacet& death2 = maxDelaunay[latest2];

          if (death1.d < death2.d) {
            if (f.d < death1.d)
              ph[2].emplace_back(FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                                 FiltratedSimplex{{death1.f[0], death1.f[1], death1.f[2]}, death1.d});
            latest[UF.find(v1)] = latest2;
          }
          else if (death2.d < death1.d) {
            if (f.d < death2.d)
              ph[2].emplace_back(FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                                 FiltratedSimplex{{death2.f[0], death2.f[1], death2.f[2]}, death2.d});
            latest[UF.find(v1)] = latest1;
          }
        }
        else // this is a facet from the minimal spanning acycle
          msa.push_back(f);
      }
    }

    void compute2PH(UnionFind &UF, MultidimensionalDiagram &ph, std::vector<Generator2> &generators2) {
      std::vector<int> latest(maxDelaunay.size());
      std::iota(latest.begin(), latest.end(), 0);

      std::vector<std::vector<Facet>> elementary_generators(N_c);
      for (const FiltratedQuadFacet &f : hyperUrquhart) {
        if (UF.find(f.c1) != UF.find(f.c2)) {
          if (maxDelaunay[UF.find(f.c1)].d != inf)
            elementary_generators[UF.find(f.c1)].push_back(f.f);
          if (maxDelaunay[UF.find(f.c2)].d != inf)
            elementary_generators[UF.find(f.c2)].push_back(f.f);
        }
      }

      for (FiltratedQuadFacet const& f : hyperUrquhart) { //sorted by decreasing order
        const int v1 = UF.find(f.c1);
        const int v2 = UF.find(f.c2);
        if (v1 != v2) { // two distinct cavities: merge them by deleting the facet
          UF.merge(v1, v2);

          const int latest1 = latest[v1];
          const int latest2 = latest[v2];
          const FiltratedFacet& death1 = maxDelaunay[latest1];
          const FiltratedFacet& death2 = maxDelaunay[latest2];

          if (death1.d < death2.d) {
            if (f.d < death1.d) {
              ph[2].emplace_back(FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                                 FiltratedSimplex{{death1.f[0], death1.f[1], death1.f[2]}, death1.d});
              generators2.push_back({elementary_generators[latest1], {f.d, death1.d}});
            }
            latest[UF.find(v1)] = latest2;
            GPH_HASHSET generator (elementary_generators[latest2].begin(), elementary_generators[latest2].end(), elementary_generators[latest2].size());
            for (Facet const& f_ : elementary_generators[latest1]) {
              auto it = generator.find(f_);
              if (it == generator.end())
                generator.insert(f_);
              else
                generator.erase(it);
            }
            elementary_generators[latest2].assign(generator.begin(), generator.end());
          }
          else if (death2.d < death1.d) {
            if (f.d < death2.d) {
              ph[2].emplace_back(FiltratedSimplex{{f.f[0], f.f[1], f.f[2]}, f.d},
                                 FiltratedSimplex{{death2.f[0], death2.f[1], death2.f[2]}, death2.d});
              generators2.push_back({elementary_generators[latest2], {f.d, death2.d}});
            }
            latest[UF.find(v1)] = latest1;
            GPH_HASHSET generator (elementary_generators[latest1].begin(), elementary_generators[latest1].end(), elementary_generators[latest1].size());
            for (Facet const& f_ : elementary_generators[latest2]) {
              auto it = generator.find(f_);
              if (it == generator.end())
                generator.insert(f_);
              else
                generator.erase(it);
            }
            elementary_generators[latest1].assign(generator.begin(), generator.end());
          }
        }
        else // this is a facet from the minimal spanning acycle
          msa.push_back(f);
      }
    }

    void connectivity() {
      //edges to facets adjacency
      msa_edges.reserve(1.15 * N_c); //this is an estimation of the number of edges in Delaunay
      for (unsigned i = 0; i<msa.size(); ++i) {
        const FiltratedQuadFacet f = msa[i];
        const double d1 = squaredDistance(f.f[0], f.f[1]);
        const double d2 = squaredDistance(f.f[0], f.f[2]);
        const double d3 = squaredDistance(f.f[1], f.f[2]);
        auto& edge1 = msa_edges[std::make_pair(f.f[0],f.f[1])];
        edge1.first.reserve(4);
        edge1.first.push_back(i);
        edge1.second |= d1 > d2 && d1 > d3;
        auto& edge2 = msa_edges[std::make_pair(f.f[0],f.f[2])];
        edge2.first.reserve(4);
        edge2.first.push_back(i);
        edge2.second |= d2 > d1 && d2 > d3;
        auto& edge3 = msa_edges[std::make_pair(f.f[1],f.f[2])];
        edge3.first.reserve(4);
        edge3.first.push_back(i);
        edge3.second |= d3 > d1 && d3 > d2;
      }
    }

    void computePolygons(MultidimensionalDiagram &ph)  {
      const unsigned N_msa = msa.size();
      UF_msa = UnionFind(N_msa);
      maxDelaunay1.resize(N_msa, {{},0});

      for (const auto& [e, val] : msa_edges) {
        auto& [adj_f, isNotUrquhart] = val;
        if (isNotUrquhart) { //not UG edge
          if (adj_f.size() == 1)
            maxDelaunay1[UF_msa.find(adj_f[0])].d = inf;
          else if (adj_f.size() == 2) {
            const int poly1 = UF_msa.find(adj_f[0]);
            const int poly2 = UF_msa.find(adj_f[1]);
            const double d = sqrt(squaredDistance(e.first, e.second));
            maxDelaunay1[UF_msa.mergeRet(poly1, poly2)] = max(FiltratedEdge{e, d},
                                                               max(maxDelaunay1[poly1], maxDelaunay1[poly2]));
          }
          else { //non-manifold junction edge
            const double d = sqrt(squaredDistance(e.first, e.second));
            critical1.push_back({e, d});
            for (auto const& f_id : adj_f) {
              const int poly = UF_msa.find(f_id);
              maxDelaunay1[poly] = max(FiltratedEdge{e,d}, maxDelaunay1[poly]);
            }
          }
        }
        else
          urquhart.push_back({e, sqrt(squaredDistance(e.first, e.second))});
      }

      TTK_PSORT(nThreadsSort_, urquhart.begin(), urquhart.end(), [](const FiltratedEdge &a, const FiltratedEdge &b) {
        return a.d < b.d;
      });
      UnionFind UF_p(N_p);
      ph[0].reserve(N_p-1);
      critical1.reserve(urquhart.size()-N_p+1);
      for (FiltratedEdge const& e : urquhart) {
        if(UF_p.find(e.e.first) != UF_p.find(e.e.second)) { //we know e is a EMST edge
          UF_p.merge(e.e.first, e.e.second);
          ph[0].emplace_back(FiltratedSimplex{{-1},0.}, FiltratedSimplex{{e.e.first, e.e.second}, e.d});
        }
        else
          critical1.emplace_back(e);
      }
    }

    void compute1PH(MultidimensionalDiagram &ph) {
      const unsigned N_msa = msa.size();
      std::vector<int> polys;
      for (unsigned x=0; x<N_msa; ++x) {
        if (UF_msa.isRoot(x) && maxDelaunay1[x].d<inf)
          polys.push_back(x);
      }

      TTK_PSORT(nThreadsSort_, polys.begin(), polys.end(), [&](const int x1, const int x2) {
        return maxDelaunay1[x1].d < maxDelaunay1[x2].d;
      });

      std::vector<int> criticalIndices(critical1.size());
      std::iota(criticalIndices.begin(), criticalIndices.end(), 0);
      TTK_PSORT(nThreadsSort_, criticalIndices.begin(), criticalIndices.end(), [&](const int e1_id, const int e2_id) {
        return critical1[e1_id].d < critical1[e2_id].d;
      });
      std::vector<int> criticalOrder(critical1.size());
      for (unsigned i=0; i<criticalIndices.size(); ++i)
        criticalOrder[criticalIndices[i]] = i;

      std::vector<std::vector<int>> poly_to_crit(N_msa);
      poly_to_crit.resize(N_msa);
      for (const int poly : polys)
        poly_to_crit[poly].reserve(3);
      for (unsigned i=0; i<critical1.size(); ++i) {
        const FiltratedEdge &e = critical1[i];
        for (const int poly : msa_edges[e.e].first) {
          if (maxDelaunay1[UF_msa.find(poly)].d<inf) {
            auto &neighbors = poly_to_crit[UF_msa.find(poly)];
            auto it = std::find(neighbors.begin(), neighbors.end(), criticalOrder[i]);
            if (it == neighbors.end())
              neighbors.push_back(criticalOrder[i]);
            else
              neighbors.erase(it);
          }
        }
      }

      std::vector<int> partner_edge(critical1.size(), -1);
      for (const int poly : polys) {
        std::set boundary(poly_to_crit[poly].begin(), poly_to_crit[poly].end());
        while (true) {
          const int youngest_id = *boundary.rbegin();
          if (partner_edge[youngest_id] == -1) {
            partner_edge[youngest_id] = poly;
            const FiltratedEdge& e = critical1[criticalIndices[youngest_id]];
            const FiltratedEdge& death = maxDelaunay1[poly];
            if (e.d < death.d)
              ph[1].emplace_back(FiltratedSimplex{{e.e.first, e.e.second}, e.d},
                                 FiltratedSimplex{{death.e.first, death.e.second}, death.d});
            break;
          }
          else {
            for (const int crit : poly_to_crit[partner_edge[youngest_id]]) {
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
    }

    void compute1PH(MultidimensionalDiagram &ph, std::vector<Generator1> &generators1) {
      const unsigned N_msa = msa.size();
      std::vector<int> polys;
      for (unsigned x=0; x<N_msa; ++x) {
        if (UF_msa.isRoot(x) && maxDelaunay1[x].d<inf)
          polys.push_back(x);
      }

      TTK_PSORT(nThreadsSort_, polys.begin(), polys.end(), [&](const int x1, const int x2) {
        return maxDelaunay1[x1].d < maxDelaunay1[x2].d;
      });

      std::vector<int> criticalIndices(critical1.size());
      std::iota(criticalIndices.begin(), criticalIndices.end(), 0);
      TTK_PSORT(nThreadsSort_, criticalIndices.begin(), criticalIndices.end(), [&](const int e1_id, const int e2_id) {
        return critical1[e1_id].d < critical1[e2_id].d;
      });
      std::vector<int> criticalOrder(critical1.size());
      for (unsigned i=0; i<criticalIndices.size(); ++i)
        criticalOrder[criticalIndices[i]] = i;

      std::vector<std::vector<int>> poly_to_crit(N_msa);
      poly_to_crit.resize(N_msa);
      for (const int poly : polys)
        poly_to_crit[poly].reserve(3);
      for (unsigned i=0; i<critical1.size(); ++i) {
        const FiltratedEdge &e = critical1[i];
        for (const int poly : msa_edges[e.e].first) {
          if (maxDelaunay1[UF_msa.find(poly)].d<inf) {
            auto &neighbors = poly_to_crit[UF_msa.find(poly)];
            auto it = std::find(neighbors.begin(), neighbors.end(), criticalOrder[i]);
            if (it == neighbors.end())
              neighbors.push_back(criticalOrder[i]);
            else
              neighbors.erase(it);
          }
        }
      }

      std::vector<std::vector<Edge>> elementary_generators(maxDelaunay1.size());
      auto action = [&](const FiltratedEdge &e) {
        auto adj_f = msa_edges[e.e].first;
        for (int &f : adj_f)
          f = UF_msa.find(f);
        std::sort(adj_f.begin(), adj_f.end());
        const auto last = std::unique(adj_f.begin(), adj_f.end());
        for (auto it = adj_f.begin(); it != last; ++it) {
          if (maxDelaunay1[UF_msa.find(*it)].d != inf)
            elementary_generators[UF_msa.find(*it)].push_back(e.e);
        }
      };
      for (auto const& e : urquhart)
        action(e);
      for (auto const& e : critical1) {
        if (msa_edges[e.e].second) // if critical but not urquhart -> non manifold junctions
          action(e);
      }

      std::vector<int> partner_edge(critical1.size(), -1);
      for (const int poly : polys) {
        std::set boundary(poly_to_crit[poly].begin(), poly_to_crit[poly].end());
        std::vector<Edge>& generator = elementary_generators[poly];
        while (true) {
          const int youngest_id = *boundary.rbegin();
          if (partner_edge[youngest_id] == -1) {
            partner_edge[youngest_id] = poly;
            const FiltratedEdge& e = critical1[criticalIndices[youngest_id]];
            const FiltratedEdge& death = maxDelaunay1[poly];
            if (e.d < death.d) {
              ph[1].emplace_back(FiltratedSimplex{{e.e.first, e.e.second}, e.d},
                                 FiltratedSimplex{{death.e.first, death.e.second}, death.d});
              generators1.emplace_back(generator, std::make_pair(e.d, death.d));
            }
            break;
          }
          else {
            for (const int crit : poly_to_crit[partner_edge[youngest_id]]) {
              const auto it = boundary.find(crit);
              if (it == boundary.end())
                boundary.insert(crit);
              else
                boundary.erase(it);
            }
            for (Edge const& e : elementary_generators[partner_edge[youngest_id]]) {
              auto it = std::find(generator.begin(), generator.end(), e);
              if (it == generator.end())
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

  inline void runDelaunayRipsPersistenceDiagram3(rpd::PointCloud const& points, MultidimensionalDiagram &diagram, int threads=1) {
    PointCloud<3> p(points.size());
    for (unsigned i = 0; i < points.size(); ++i)
      p[i] = {points[i][0], points[i][1], points[i][2]};
    DRPersistence3 drpd(p, threads, threads);
    drpd.computeDelaunayRipsPersistence(diagram);
  }

  inline void runDelaunayRipsPersistenceDiagram3(rpd::PointCloud const& points, MultidimensionalDiagram &diagram, std::vector<Generator1> &generators1, std::vector<Generator2> &generators2, int threads=1) {
    PointCloud<3> p(points.size());
    for (unsigned i = 0; i < points.size(); ++i)
      p[i] = {points[i][0], points[i][1], points[i][2]};
    DRPersistence3 drpd(p, threads, threads);
    drpd.computeDelaunayRipsPersistence(diagram, generators1, generators2);
  }

}

#endif