#include <FastRipsPersistenceDiagram2.h>

#ifdef TTK_ENABLE_CGAL

using namespace ttk::rpd;

FastRipsPersistenceDiagram2::FastRipsPersistenceDiagram2(
  const PointCloud &points)
  : n_(points.size()) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("GeometricRipsPD");
  tm_.reStart();

  points_ = std::vector<std::pair<Point, unsigned>>(n_);
  for(unsigned i = 0; i < n_; ++i)
    points_[i] = std::make_pair(Point(points[i][0], points[i][1]), i);
  printMsg("Loaded", 0., tm_.getElapsedTime());

  computeDelaunay();
}

FastRipsPersistenceDiagram2::FastRipsPersistenceDiagram2(float *data, int n)
  : n_(n) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("GeometricRipsPD");
  tm_.reStart();

  points_ = std::vector<std::pair<Point, unsigned>>(n_);
  for(unsigned i = 0; i < n_; ++i)
    points_[i] = std::make_pair(Point(data[2 * i], data[2 * i + 1]), i);
  printMsg("Loaded", 0., tm_.getElapsedTime());

  computeDelaunay();
}

template <typename T>
void FastRipsPersistenceDiagram2::compute0Persistence(T &ph0,
                                                      bool parallelSort) {
  ph0 = T(0);

  // keep only the Urquhart edges, maintain polygon structure
  UnionFind UF(nFaces_);
  std::vector<FiltratedEdge> max_delaunay(nFaces_, FiltratedEdge{{-1, -1}, 0.});
  computeUrquhart(UF, max_delaunay, parallelSort);

  // compute EMST with Kruskal algorithm
  UnionFind UF_p(n_);
  for(FiltratedQuadEdge const &e : urquhart_) {
    if(UF_p.find(e.e.first)
       != UF_p.find(e.e.second)) { // we know e is a EMST edge
      UF_p.merge(e.e.first, e.e.second);
      add0Pair(e, ph0);
    }
  }
  if constexpr(std::is_same_v<T, Diagram>)
    ph0.emplace_back(
      FiltratedSimplex{{-1}, 0.}, FiltratedSimplex{{-1}, inf}); // infinite pair

  printMsg("MST computed", 0., tm_.getElapsedTime());
}
template void
  FastRipsPersistenceDiagram2::compute0Persistence(Diagram &ph0,
                                                   bool parallelSort);
template void
  FastRipsPersistenceDiagram2::compute0Persistence(EdgeSet &ph0,
                                                   bool parallelSort);

template <typename T>
void FastRipsPersistenceDiagram2::computeDelaunayRips0And1Persistence(
  T &ph, bool parallelSort) {
  if constexpr(std::is_same_v<T, MultidimensionalDiagram>)
    ph = MultidimensionalDiagram(2);

  // keep only the Urquhart edges, maintain polygon structure
  UnionFind UF(nFaces_);
  deathPoly_.resize(nFaces_, {{-1, -1}, 0.});
  computeUrquhart(UF, deathPoly_, parallelSort);

  // compute EMST with Kruskal algorithm
  UnionFind UF_p(n_);
  std::vector<FiltratedQuadEdge> critical(0);
  for(FiltratedQuadEdge const &e : urquhart_) {
    if(UF_p.find(e.e.first)
       != UF_p.find(e.e.second)) { // we know e is a EMST edge
      UF_p.merge(e.e.first, e.e.second);
      add0Pair(e, ph[0]);
    } else { // we know e is a UG-EMST edge, i.e., it creates a cycle
      critical.push_back(e);
      if constexpr(std::is_same_v<T, EdgeSets3>)
        ph[1].emplace_back(e.e);
    }
  }
  if constexpr(std::is_same_v<T, MultidimensionalDiagram>)
    ph[0].emplace_back(
      FiltratedSimplex{{-1}, 0.}, FiltratedSimplex{{-1}, inf}); // infinite pair

  printMsg("MST computed", 0., tm_.getElapsedTime());

  for(FiltratedQuadEdge &e : urquhart_) {
    e.f1 = UF.find(e.f1); // final path compression
    e.f2 = UF.find(e.f2); // final path compression
  }
  for(FiltratedQuadEdge &e : critical) {
    e.f1 = UF.find(e.f1); // final path compression
    e.f2 = UF.find(e.f2); // final path compression
  }

  if constexpr(std::is_same_v<T, MultidimensionalDiagram>)
    compute1PH(critical, UF, ph);
  else if constexpr(std::is_same_v<T, EdgeSets3>) {
    for(unsigned poly = 0; poly < deathPoly_.size(); ++poly) {
      if(UF.isRoot(poly) && deathPoly_[poly].d != inf)
        ph[2].emplace_back(deathPoly_[poly].e);
    }
  }
}
template void FastRipsPersistenceDiagram2::computeDelaunayRips0And1Persistence(
  MultidimensionalDiagram &ph, bool parallelSort);
template void FastRipsPersistenceDiagram2::computeDelaunayRips0And1Persistence(
  EdgeSets3 &ph, bool parallelSort);

template <typename T>
void FastRipsPersistenceDiagram2::computeRips0And1Persistence(
  T &ph, bool parallelSort, bool parallelMML) {
  if constexpr(std::is_same_v<T, MultidimensionalDiagram>)
    ph = MultidimensionalDiagram(2);

  // compute k-d tree
  Tree tree;
  for(auto const &p : points_)
    tree.insert(p.first);
  printMsg("k-d tree initialized", 0., tm_.getElapsedTime());

  // keep only the Urquhart edges, maintain polygon structure
  UnionFind UF(nFaces_);
  std::vector<FiltratedEdge> max_delaunay(nFaces_, FiltratedEdge{{-1, -1}, 0.});
  computeUrquhart(UF, max_delaunay, parallelSort);

  // compute EMST with Kruskal algorithm
  UnionFind UF_p(n_);
  std::vector<FiltratedQuadEdge> critical(0); // RNG-EMST edges

  for(FiltratedQuadEdge const &e : urquhart_) {
    if(UF_p.find(e.e.first)
       != UF_p.find(e.e.second)) { // we know e is a EMST edge
      UF_p.merge(e.e.first, e.e.second);
      rng_.push_back(e);
      add0Pair(e, ph[0]);
    } else { // we know e is a UG-EMST edge
      // check if e is RNG
      if(isLensEmpty(points_[e.e.first].first, points_[e.e.second].first, tree,
                     e.d)) { // RNG edge
        critical.push_back(e);
        rng_.push_back(e);
        if constexpr(std::is_same_v<T,
                                    EdgeSets3> || std::is_same_v<T, EdgeSets4>)
          ph[1].emplace_back(e.e);
      } else { // not RNG edge : merge neighboring polygons
        const int poly1 = UF.find(e.f1);
        const int poly2 = UF.find(e.f2);
        UF.merge(poly1, poly2);
        max_delaunay[UF.find(poly1)]
          = max(FiltratedEdge{e.e, e.d},
                max(max_delaunay[poly1], max_delaunay[poly2]));
      }
    }
  }
  if constexpr(std::is_same_v<T, MultidimensionalDiagram>)
    ph[0].emplace_back(
      FiltratedSimplex{{-1}, 0.}, FiltratedSimplex{{-1}, inf}); // infinite pair

  printMsg("MST and RNG computed", 0., tm_.getElapsedTime());

  // polygon reindexation
  std::vector<int> index_polys(0);
  reindexPolygons(UF, max_delaunay, index_polys);
  for(FiltratedQuadEdge &e : rng_) {
    e.f1 = index_polys[UF.find(e.f1)]; // final path compression
    e.f2 = index_polys[UF.find(e.f2)]; // final path compression
  }
  for(FiltratedQuadEdge &e : critical) {
    e.f1 = index_polys[UF.find(e.f1)]; // final path compression
    e.f2 = index_polys[UF.find(e.f2)]; // final path compression
  }

  if constexpr(std::is_same_v<T, MultidimensionalDiagram>) {
    computePolygonRipsDeath(parallelMML, UF, index_polys);
    printMsg("MML edges computed", 0., tm_.getElapsedTime());
    UnionFind UF_poly(deathPoly_.size());
    compute1PH(critical, UF_poly, ph);
  } else if constexpr(std::is_same_v<T, EdgeSets3>) {
    computePolygonRipsDeath(parallelMML, UF, index_polys);
    printMsg("MML edges computed", 0., tm_.getElapsedTime());
    for(FiltratedEdge const &poly : deathPoly_) {
      if(poly.d != inf)
        ph[2].emplace_back(poly.e);
    }
  } else if constexpr(std::is_same_v<T, EdgeSets4>) {
    executePolygonPairCells(parallelMML, UF, index_polys, ph);
    printMsg("Local cascades computed", 0., tm_.getElapsedTime());
  }
}
template void FastRipsPersistenceDiagram2::computeRips0And1Persistence(
  MultidimensionalDiagram &ph, bool parallelSort, bool parallelMML);
template void FastRipsPersistenceDiagram2::computeRips0And1Persistence(
  EdgeSets3 &ph, bool parallelSort, bool parallelMML);
template void FastRipsPersistenceDiagram2::computeRips0And1Persistence(
  EdgeSets4 &ph, bool parallelSort, bool parallelMML);

void FastRipsPersistenceDiagram2::computeDelaunay() {
  // compute Delaunay and initialize triangle IDs (including infinite triangles
  // to deal with the convex hull)
  delaunay_ = Delaunay(points_.begin(), points_.end());
  int k = 0;
  for(auto const &f : delaunay_.all_face_handles())
    f->id() = k++;
  nFaces_ = k;

  printMsg("Delaunay triangulation computed", 0., tm_.getElapsedTime());
}

void FastRipsPersistenceDiagram2::computeUrquhart(
  UnionFind &UF, std::vector<FiltratedEdge> &maxDelaunay, bool parallelSort) {
  for(Delaunay::Edge const &e : delaunay_.finite_edges()) {
    const unsigned a = e.first->vertex((e.second + 1) % 3)->info();
    const unsigned b = e.first->vertex((e.second + 2) % 3)->info();
    const double d2
      = CGAL::squared_distance(points_[a].first, points_[b].first);
    const double d = sqrt(d2);
    const Delaunay::Edge e_m = delaunay_.mirror_edge(e);

    // check if edge is Urquhart
    bool is_urquhart = true;
    const Point p_k = e.first->vertex(e.second)->point();
    if(!delaunay_.is_infinite(e.first)
       && CGAL::squared_distance(points_[a].first, p_k) < d2
       && CGAL::squared_distance(points_[b].first, p_k) < d2)
      is_urquhart = false;
    else {
      const Point p_l = e_m.first->vertex(e_m.second)->point();
      if(!delaunay_.is_infinite(e_m.first)
         && CGAL::squared_distance(points_[a].first, p_l) < d2
         && CGAL::squared_distance(points_[b].first, p_l) < d2)
        is_urquhart = false;
    }

    // maintain polygon structure
    if(is_urquhart) // UG edge
      urquhart_.push_back(
        FiltratedQuadEdge{{a, b}, e.first->id(), e_m.first->id(), d});
    else { // nUG edge: maintain UF
      const int poly1 = UF.find(e.first->id());
      const int poly2 = UF.find(e_m.first->id());
      maxDelaunay[UF.mergeRet(poly1, poly2)] = max(
        FiltratedEdge{{a, b}, d}, max(maxDelaunay[poly1], maxDelaunay[poly2]));
    }

    // detect infinite polygons (going beyond convex hull)
    if(delaunay_.is_infinite(e.first))
      maxDelaunay[UF.find(e.first->id())].d = inf;
    else if(delaunay_.is_infinite(e_m.first))
      maxDelaunay[UF.find(e_m.first->id())].d = inf;
  }
  printMsg("Urquhart graph computed", 0., tm_.getElapsedTime());

  if(parallelSort)
    TTK_PSORT(
      globalThreadNumber_, urquhart_.begin(), urquhart_.end(),
      [](FiltratedQuadEdge e1, FiltratedQuadEdge e2) { return e1.d < e2.d; })
  else
    std::sort(
      urquhart_.begin(), urquhart_.end(),
      [](FiltratedQuadEdge e1, FiltratedQuadEdge e2) { return e1.d < e2.d; });
  printMsg("Urquhart graph sorted", 0., tm_.getElapsedTime());
}

void FastRipsPersistenceDiagram2::compute1PH(
  std::vector<FiltratedQuadEdge> const &critical,
  UnionFind &UF,
  MultidimensionalDiagram &ph) {
  std::vector<int> latest(deathPoly_.size());
  std::iota(latest.begin(), latest.end(), 0);
  birthPoly_.resize(deathPoly_.size());

  for(auto it = critical.rbegin(); it != critical.rend(); ++it) {
    const FiltratedQuadEdge e = *it;

    const int v1 = UF.find(e.f1);
    const int v2 = UF.find(e.f2);
    UF.merge(v1, v2);

    const int latest1 = latest[v1];
    const int latest2 = latest[v2];

    if(deathPoly_[latest1].d < deathPoly_[latest2].d) {
      ph[1].emplace_back(FiltratedSimplex{{e.e.first, e.e.second}, e.d},
                         FiltratedSimplex{{deathPoly_[latest1].e.first,
                                           deathPoly_[latest1].e.second},
                                          deathPoly_[latest1].d});
      birthPoly_[latest1] = e.d;
      latest[UF.find(v1)] = latest2;
    } else {
      ph[1].emplace_back(FiltratedSimplex{{e.e.first, e.e.second}, e.d},
                         FiltratedSimplex{{deathPoly_[latest2].e.first,
                                           deathPoly_[latest2].e.second},
                                          deathPoly_[latest2].d});
      birthPoly_[latest2] = e.d;
      latest[UF.find(v1)] = latest1;
    }
  }

  printMsg("Dual Kruskal complete", 0., tm_.getElapsedTime());
}

bool FastRipsPersistenceDiagram2::isLensEmpty(Point const &p1,
                                              Point const &p2,
                                              Tree const &tree,
                                              double const &d) {
  const Fuzzy_sphere fs(p1, d);
  std::vector<Point> ball;
  tree.search(back_inserter(ball), fs);
  for(Point const &p : ball) {
    if(p != p1 && p != p2) {
      if(CGAL::squared_distance(p, p2) < d * d)
        return false;
    }
  }
  return true;
}

bool FastRipsPersistenceDiagram2::isRightSemiLensEmpty(Point const &p1,
                                                       Point const &p2,
                                                       Tree const &tree) {
  const double d2 = CGAL::squared_distance(p1, p2);
  const Fuzzy_sphere fs(p1, sqrt(d2));
  std::vector<Point> ball;
  tree.search(back_inserter(ball), fs);
  for(Point const &p : ball) {
    if(CGAL::squared_distance(p, p2) < d2 && CGAL::right_turn(p, p1, p2))
      return false;
  }
  return true;
}

void FastRipsPersistenceDiagram2::computePolygonRipsDeath(
  bool parallel, UnionFind &UF, std::vector<int> const &indexPolys) {
  const unsigned N_polys = deathPoly_.size();
  std::vector<std::vector<int>> vertices(N_polys);

  // store coarsely vertices in each polygon
  for(const FiltratedQuadEdge &e : rng_) {
    vertices[e.f1].push_back(e.e.first);
    vertices[e.f1].push_back(e.e.second);
    vertices[e.f2].push_back(e.e.first);
    vertices[e.f2].push_back(e.e.second);
  }
  std::vector<Delaunay::Face_handle> repr_face(N_polys);
  for(auto const &fh : delaunay_.finite_face_handles()) {
    if(UF.isRoot(fh->id()))
      repr_face[indexPolys[fh->id()]] = fh;
  }

  // loop on polygons
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for if(parallel)
#else
  TTK_FORCE_USE(parallel);
#endif // TTK_ENABLE_OPENMP
  for(unsigned poly = 0; poly < N_polys; ++poly) {
    if(deathPoly_[poly].d != inf) { // deal only with true polygons
      const double bound_max = deathPoly_[poly].d;
      const double bound_min = sqrt(3) / 2 * bound_max;

      // find involved vertices
      sort(vertices[poly].begin(), vertices[poly].end());
      auto last = unique(vertices[poly].begin(), vertices[poly].end());
      vertices[poly].erase(last, vertices[poly].end());

      std::vector<Point> pts(0);
      std::vector<int> global(0);
      for(int const &v : vertices[poly]) {
        pts.push_back(points_[v].first);
        global.push_back(v);
      }
      const unsigned N_pts = pts.size();
      const Tree loc_tree(pts.begin(), pts.end());

      // enumerate edges and do stuff
      FiltratedEdge death{{-1, -1}, inf};
      for(unsigned i = 0; i < N_pts - 1; ++i) {
        for(unsigned j = i + 1; j < N_pts; ++j) {
          const double d2 = CGAL::squared_distance(pts[i], pts[j]);
          const double d = sqrt(d2);
          if(bound_min <= d && d <= bound_max
             && d < death
                      .d) { // edge length interval allowing to be a candidate
            const Point p_i = pts[i];
            const Point p_j = pts[j];

            const Fuzzy_sphere fs(p_i, d);
            std::vector<Point> ball;
            loc_tree.search(back_inserter(ball), fs);

            // check if edge is a 2-edge (i.e., both semi-lenses are not empty)
            std::vector<Point> lens_l;
            std::vector<Point> lens_r;
            for(auto const &p : ball) {
              if(p != p_i && p != p_j) {
                if(CGAL::squared_distance(p, p_j) < d2) {
                  if(CGAL::right_turn(p_i, p_j, p))
                    lens_r.push_back(p);
                  else
                    lens_l.push_back(p);
                }
              }
            }

            // if 2-edge, check if it is expandable
            if(!lens_l.empty() && !lens_r.empty()) { // 2-edge
              bool l_expandable = false, r_expandable = false;
              for(Point const &p : lens_l) {
                if(!isRightSemiLensEmpty(p_i, p, loc_tree))
                  continue;
                else if(isRightSemiLensEmpty(p, p_j, loc_tree)) {
                  l_expandable = true;
                  break;
                }
              }
              if(!l_expandable)
                continue;
              for(Point const &p : lens_r) {
                if(!isRightSemiLensEmpty(p_j, p, loc_tree))
                  continue;
                else if(isRightSemiLensEmpty(p, p_i, loc_tree)) {
                  r_expandable = true;
                  break;
                }
              }
              if(r_expandable) { // expandable 2-edge, check if it is a diagonal
                                 // of the current polygon (costly)
                if(indexPolys[UF.find(
                     delaunay_
                       .locate(CGAL::midpoint(p_i, p_j), repr_face[poly])
                       ->id())]
                   == int(poly))
                  death = FiltratedEdge{{global[i], global[j]}, d};
              }
            }
          }
        }
      }
      deathPoly_[poly] = death;
    }
  }
}

void FastRipsPersistenceDiagram2::reindexPolygons(
  UnionFind const &UF,
  std::vector<FiltratedEdge> const &maxDelaunay,
  std::vector<int> &indexPolys) {
  // reindex polygons and find those that are beyond convex hull
  indexPolys = std::vector<int>(nFaces_, -1);
  int N_polys = 0;
  for(unsigned x = 0; x < nFaces_; ++x) {
    if(UF.isRoot(x)) {
      indexPolys[x] = N_polys;
      deathPoly_.push_back(maxDelaunay[x]);
      N_polys++;
    }
  }
}

void FastRipsPersistenceDiagram2::pComputePolygonRipsDeath(
  UnionFind &UF, std::vector<int> const &indexPolys) {
  const unsigned N_polys = deathPoly_.size();
  std::vector<std::vector<int>> vertices(N_polys);

  // store coarsely vertices in each polygon
  for(FiltratedQuadEdge const &e : rng_) {
    vertices[e.f1].push_back(e.e.first);
    vertices[e.f1].push_back(e.e.second);
    vertices[e.f2].push_back(e.e.first);
    vertices[e.f2].push_back(e.e.second);
  }
  std::vector<Delaunay::Face_handle> repr_face(N_polys);
  for(auto const &fh : delaunay_.finite_face_handles()) {
    if(UF.isRoot(fh->id()))
      repr_face[indexPolys[fh->id()]] = fh;
  }

  std::vector<std::pair<FiltratedEdge, int>> candidates;
  std::vector<Tree> trees(N_polys);

  // loop on polygons
  for(unsigned poly = 0; poly < N_polys; ++poly) {
    if(deathPoly_[poly].d != inf) { // deal only with true polygons
      const double bound_max = deathPoly_[poly].d;
      const double bound_min = sqrt(3) / 2 * bound_max;

      // find involved vertices
      sort(vertices[poly].begin(), vertices[poly].end());
      auto last = unique(vertices[poly].begin(), vertices[poly].end());
      vertices[poly].erase(last, vertices[poly].end());

      std::vector<Point> pts(0);
      std::vector<int> global(0);
      for(int const &v : vertices[poly]) {
        pts.push_back(points_[v].first);
        global.push_back(v);
      }
      const unsigned N_pts = pts.size();
      trees[poly].insert(pts.begin(), pts.end());

      // enumerate edges and do stuff
      for(unsigned i = 0; i < N_pts - 1; ++i) {
        for(unsigned j = i + 1; j < N_pts; ++j) {
          const double d2 = CGAL::squared_distance(pts[i], pts[j]);
          const double d = sqrt(d2);
          if(bound_min <= d && d <= bound_max) // edge length interval allowing
                                               // to be a candidate
            candidates.emplace_back(
              FiltratedEdge{{global[i], global[j]}, d}, poly);
        }
      }
    }
  }

  std::vector<bool> validCandidate(candidates.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif // TTK_ENABLE_OPENMP
  for(unsigned i = 0; i < candidates.size(); ++i) {
    const Point p_i = points_[candidates[i].first.e.first].first;
    const Point p_j = points_[candidates[i].first.e.second].first;
    const double d = candidates[i].first.d;
    const int poly = candidates[i].second;

    const Fuzzy_sphere fs(p_i, d);
    std::vector<Point> ball;
    trees[poly].search(back_inserter(ball), fs);

    // check if edge is a 2-edge (i.e., both semi-lenses are not empty)
    std::vector<Point> lens_l;
    std::vector<Point> lens_r;
    for(auto const &p : ball) {
      if(p != p_i && p != p_j) {
        if(CGAL::squared_distance(p, p_j) < d * d) {
          if(CGAL::right_turn(p_i, p_j, p))
            lens_r.push_back(p);
          else
            lens_l.push_back(p);
        }
      }
    }

    // if 2-edge, check if it is expandable
    if(!lens_l.empty() && !lens_r.empty()) { // 2-edge
      bool l_expandable = false, r_expandable = false;
      for(Point const &p : lens_l) {
        if(!isRightSemiLensEmpty(p_i, p, trees[poly]))
          continue;
        else if(isRightSemiLensEmpty(p, p_j, trees[poly])) {
          l_expandable = true;
          break;
        }
      }
      if(!l_expandable)
        continue;
      for(Point const &p : lens_r) {
        if(!isRightSemiLensEmpty(p_j, p, trees[poly]))
          continue;
        else if(isRightSemiLensEmpty(p, p_i, trees[poly])) {
          r_expandable = true;
          break;
        }
      }
      if(r_expandable) { // expandable 2-edge, check if it is a diagonal of the
                         // current polygon (costly)
        if(indexPolys[UF.find(
             delaunay_.locate(CGAL::midpoint(p_i, p_j), repr_face[poly])->id())]
           == int(poly))
          validCandidate[i] = true;
      }
    }
  }

  deathPoly_ = std::vector<FiltratedEdge>(N_polys, {{-1, -1}, inf});
  for(unsigned i = 0; i < candidates.size(); ++i) {
    if(validCandidate[i]) {
      const int poly = candidates[i].second;
      if(candidates[i].first.d < deathPoly_[poly].d)
        deathPoly_[poly] = candidates[i].first;
    }
  }
}

void FastRipsPersistenceDiagram2::exportRips1Generators(
  std::vector<Generator1> &generators) {
  const unsigned N_polys = deathPoly_.size();

  std::vector<Generator1> generators_(N_polys);
  for(unsigned i = 0; i < deathPoly_.size(); ++i)
    generators_[i].second = {birthPoly_[i], deathPoly_[i].d};
  for(const FiltratedQuadEdge &e : (N_polys == nFaces_) ? urquhart_ : rng_) {
    if(e.f1 != e.f2) {
      if(deathPoly_[e.f1].d != inf)
        generators_[e.f1].first.push_back(e.e);
      if(deathPoly_[e.f2].d != inf)
        generators_[e.f2].first.push_back(e.e);
    }
  }

  generators.resize(0);
  for(Generator1 const &g : generators_) {
    if(!g.first.empty())
      generators.push_back(g);
  }
}

void FastRipsPersistenceDiagram2::executePolygonPairCells(
  bool parallel,
  UnionFind &UF,
  std::vector<int> const &indexPolys,
  EdgeSets4 &ph) const {
  std::set<Edge> cascadeSet{};

  const unsigned N_polys = deathPoly_.size();
  std::vector<std::vector<int>> vertices(N_polys);

  // store coarsely vertices in each polygon
  for(const FiltratedQuadEdge &e : rng_) {
    vertices[e.f1].push_back(e.e.first);
    vertices[e.f1].push_back(e.e.second);
    vertices[e.f2].push_back(e.e.first);
    vertices[e.f2].push_back(e.e.second);
  }
  std::vector<Delaunay::Face_handle> repr_face(N_polys);
  for(auto const &fh : delaunay_.finite_face_handles()) {
    if(UF.isRoot(fh->id()))
      repr_face[indexPolys[fh->id()]] = fh;
  }

  // loop on polygons
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for if(parallel)
#else
  TTK_FORCE_USE(parallel);
#endif // TTK_ENABLE_OPENMP
  for(unsigned poly = 0; poly < N_polys; ++poly) {
    if(deathPoly_[poly].d != inf) { // deal only with true polygons
      const double bound_max = deathPoly_[poly].d;

      // find involved vertices
      sort(vertices[poly].begin(), vertices[poly].end());
      auto last = unique(vertices[poly].begin(), vertices[poly].end());
      vertices[poly].erase(last, vertices[poly].end());

      std::vector<Point> pts(0);
      std::vector<int> global(0);
      for(int const &v : vertices[poly]) {
        pts.push_back(points_[v].first);
        global.push_back(v);
      }

      PairCells pc(pts, bound_max);
      pc.run();
      pc.enrichCascades(cascadeSet, ph, global);
    }
  }

  for(const Edge &e : cascadeSet)
    ph[3].emplace_back(e);
}

#endif