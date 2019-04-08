#include <QuadrangulationSubdivision.h>

#define MODULE_S "[QuadrangulationSubdivision] "

int ttk::QuadrangulationSubdivision::subdivise(
  std::vector<Quad> &currQuads, const std::vector<Quad> &prevQuads) {

  using edgeType = std::pair<long long, long long>;
  using vertexType = std::pair<long long, Point>;
  using std::make_pair;
  std::map<edgeType, vertexType> processedEdges;

  Timer t;

  // avoid reallocation in loop, causing invalid pointers
  outputPoints_->reserve(outputPoints_->size() * 5);

  for(auto &q : prevQuads) {
    assert(q.n == 4); // magic number...

    auto i = static_cast<size_t>(q.i);
    auto j = static_cast<size_t>(q.j);
    auto k = static_cast<size_t>(q.k);
    auto l = static_cast<size_t>(q.l);

    Point *pi = &(*outputPoints_)[i];
    Point *pj = &(*outputPoints_)[j];
    Point *pk = &(*outputPoints_)[k];
    Point *pl = &(*outputPoints_)[l];

    // middles of edges
    auto midij = (*pi + *pj) * 0.5;
    auto midjk = (*pj + *pk) * 0.5;
    auto midkl = (*pk + *pl) * 0.5;
    auto midli = (*pl + *pi) * 0.5;

    // quad barycenter
    auto bary = (*pi + *pj + *pk + *pl) * 0.25;

    // order edges to avoid duplicates (ij vs. ji)
    auto ij = make_pair(std::min(q.i, q.j), std::max(q.i, q.j));
    auto jk = make_pair(std::min(q.j, q.k), std::max(q.j, q.k));
    auto kl = make_pair(std::min(q.k, q.l), std::max(q.k, q.l));
    auto li = make_pair(std::min(q.l, q.i), std::max(q.l, q.i));

    // add to outputPoints_ after computing new point coordinates to
    // avoid invalidating pointers
    if(processedEdges.find(ij) == processedEdges.end()) {
      processedEdges.insert(
        make_pair(ij, make_pair(outputPoints_->size(), midij)));
      outputPoints_->emplace_back(midij);
    }

    if(processedEdges.find(jk) == processedEdges.end()) {
      processedEdges.insert(
        make_pair(jk, make_pair(outputPoints_->size(), midjk)));
      outputPoints_->emplace_back(midjk);
    }

    if(processedEdges.find(kl) == processedEdges.end()) {
      processedEdges.insert(
        make_pair(kl, make_pair(outputPoints_->size(), midkl)));
      outputPoints_->emplace_back(midkl);
    }

    if(processedEdges.find(li) == processedEdges.end()) {
      processedEdges.insert(
        make_pair(li, make_pair(outputPoints_->size(), midli)));
      outputPoints_->emplace_back(midli);
    }

    // barycenter index in outputPoints_
    auto baryIdx = static_cast<long long>(outputPoints_->size());
    outputPoints_->emplace_back(bary);

    // add the four new quads
    currQuads.emplace_back(Quad{
      4, q.i, processedEdges[ij].first, baryIdx, processedEdges[li].first});
    currQuads.emplace_back(Quad{
      4, q.j, processedEdges[jk].first, baryIdx, processedEdges[ij].first});
    currQuads.emplace_back(Quad{
      4, q.k, processedEdges[kl].first, baryIdx, processedEdges[jk].first});
    currQuads.emplace_back(Quad{
      4, q.l, processedEdges[li].first, baryIdx, processedEdges[kl].first});
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Subdivised " << prevQuads.size() << " quads into "
        << currQuads.size() << " new quads (" << outputPoints_->size()
        << " points) in " << t.getElapsedTime() << "s" << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}

ttk::SimplexId ttk::QuadrangulationSubdivision::findNearestInputVertex(
  const ttk::SimplexId i) const {

  // current point to project
  Point *vert = &(*outputPoints_)[i];

  // holds the distance to the nearest vertex
  std::pair<float, SimplexId> nearestVertex = std::make_pair(-1.0f, -1);

  // iterate over all vertices of the input mesh, find the nearest one
  for(SimplexId j = 0; j < triangulation_->getNumberOfVertices(); j++) {

    // get vertex coordinates
    Point inMesh{};
    triangulation_->getVertexPoint(j, inMesh.x, inMesh.y, inMesh.z);

    // get square distance to vertex
    float dist = (vert->x - inMesh.x) * (vert->x - inMesh.x)
                 + (vert->y - inMesh.y) * (vert->y - inMesh.y)
                 + (vert->z - inMesh.z) * (vert->z - inMesh.z);

    if(nearestVertex.first < 0.0f || dist < nearestVertex.first) {
      nearestVertex.first = dist;
      nearestVertex.second = j;
    }
  }

  return nearestVertex.second;
}

ttk::SimplexId
  ttk::QuadrangulationSubdivision::findVerticesMiddle(const ttk::SimplexId m,
                                                      const ttk::SimplexId n) {
  Point M{}, N{};
  triangulation_->getVertexPoint(m, M.x, M.y, M.z);
  triangulation_->getVertexPoint(n, N.x, N.y, N.z);

  if(m == n) {
    return m;
  }

  // find shortest path between M and N in graph using Dijkstra
  using pq_t = std::pair<float, SimplexId>;

  // priority queue
  std::priority_queue<pq_t, std::vector<pq_t>, std::greater<pq_t>> pq;
  // vector of distance to M
  std::vector<float> dists(
    vertexNumber_, std::numeric_limits<float>::infinity());
  // vector of parent vertices
  std::vector<SimplexId> parent(vertexNumber_, -1);

  // init pipeline
  pq.push(std::make_pair(0.0f, m));
  dists[m] = 0.0f;

  while(!pq.empty()) {
    auto elem = pq.top();
    pq.pop();
    auto a = elem.second;
    Point A{};
    triangulation_->getVertexPoint(a, A.x, A.y, A.z);

    auto nneigh = triangulation_->getVertexNeighborNumber(a);
    for(SimplexId j = 0; j < nneigh; j++) {
      SimplexId b{};
      triangulation_->getVertexNeighbor(a, j, b);
      Point B{};
      triangulation_->getVertexPoint(b, B.x, B.y, B.z);
      float distAB
        = powf(B.x - A.x, 2.0f) + powf(B.y - A.y, 2.0f) + powf(B.z - A.z, 2.0f);
      if(dists[b] > dists[a] + distAB) {
        dists[b] = dists[a] + distAB;
        parent[b] = a;
        if(b == n) {
          break;
        }
        pq.push(std::make_pair(dists[b], b));
      }
    }
  }

  // vertices in shortest path from M to N
  std::vector<SimplexId> pathMN;

  SimplexId j = n;
  while(j != m) {
    pathMN.emplace_back(parent[j]);
    j = parent[j];
  }

  // middle of shortest path between M and N
  return pathMN[pathMN.size() / 2];
}

ttk::QuadrangulationSubdivision::Point
  ttk::QuadrangulationSubdivision::findProjectionInTriangle(
    const ttk::SimplexId i) const {

  // current point to project
  Point *vert = &(*outputPoints_)[i];

  // projected point into triangle
  Point proj{};
  // found a projection in one triangle
  bool success = false;
  // list of triangle IDs to test to find a potential projection
  std::queue<SimplexId> trianglesToTest;
  // list of triangle IDs already tested
  // (takes more memory to reduce computation time)
  std::vector<bool> trianglesTested(
    triangulation_->getNumberOfTriangles(), false);

  // number of triangles around nearest vertex
  SimplexId triangleNumber
    = triangulation_->getVertexTriangleNumber(nearestVertexIdentifier_[i]);
  // init pipeline by checking in every triangle around selected vertex
  for(SimplexId j = 0; j < triangleNumber; j++) {
    SimplexId ntid;
    triangulation_->getVertexTriangle(nearestVertexIdentifier_[i], j, ntid);
    trianglesToTest.push(ntid);
  }

  while(!trianglesToTest.empty()) {
    SimplexId tid = trianglesToTest.front();
    trianglesToTest.pop();

    // skip if already tested
    if(trianglesTested[tid]) {
      continue;
    }

    // get triangle vertices
    SimplexId tverts[3];
    triangulation_->getTriangleVertex(tid, 0, tverts[0]);
    triangulation_->getTriangleVertex(tid, 1, tverts[1]);
    triangulation_->getTriangleVertex(tid, 2, tverts[2]);

    // get coordinates of triangle vertices
    Point pa{}, pb{}, pc{};
    triangulation_->getVertexPoint(tverts[0], pa.x, pa.y, pa.z);
    triangulation_->getVertexPoint(tverts[1], pb.x, pb.y, pb.z);
    triangulation_->getVertexPoint(tverts[2], pc.x, pc.y, pc.z);

    // triangle normal: cross product of two edges
    Point crossP{};
    // ab, ac vectors
    Point ab = pb - pa, ac = pc - pa;
    // compute ab ^ ac
    Geometry::crossProduct(&ab.x, &ac.x, &crossP.x);
    // unitary normal vector
    Point norm = crossP / Geometry::magnitude(&crossP.x);

    Point tmp = *vert - pa;
    // projected point into triangle
    proj = *vert - norm * Geometry::dotProduct(&norm.x, &tmp.x);

    // check if projection in triangle
    if(Geometry::isPointInTriangle(&pa.x, &pb.x, &pc.x, &proj.x)) {
      success = true;
      // should we check if we have the nearest triangle?
      break;
    }

    // mark triangle as tested
    trianglesTested[tid] = true;

    // (re-)compute barycentric coords of projection
    std::vector<float> baryCoords;
    Geometry::computeBarycentricCoordinates(
      &pa.x, &pb.x, &pc.x, &proj.x, baryCoords);

    // extrema values in baryCoords
    auto extrema = std::minmax_element(baryCoords.begin(), baryCoords.end());

    // find the nearest triangle vertices (with the highest/positive
    // values in baryCoords) from proj
    std::vector<SimplexId> vertices(2);
    vertices[0] = tverts[extrema.second - baryCoords.begin()];
    for(size_t j = 0; j < baryCoords.size(); j++) {
      if(j != static_cast<size_t>(extrema.first - baryCoords.begin())
         && j != static_cast<size_t>(extrema.second - baryCoords.begin())) {
        vertices[1] = tverts[j];
      }
    }
    vertices[1] = tverts[extrema.second - baryCoords.begin()];

    // triangles to test next
    std::set<SimplexId> common_triangles;

    // look for triangles sharing the two edges with max values in
    // baryCoords
    for(auto &vert : vertices) {
      SimplexId tnum = triangulation_->getVertexTriangleNumber(vert);
      for(SimplexId j = 0; j < tnum; j++) {
        SimplexId trid;
        triangulation_->getVertexTriangle(vert, j, trid);
        if(trid == tid) {
          continue;
        }
        common_triangles.insert(trid);
      }
    }

    for(auto &ntid : common_triangles) {
      if(!trianglesTested[ntid]) {
        trianglesToTest.push(ntid);
      }
    }
  }

  if(!success) {
    // replace proj by the nearest vertex?
    triangulation_->getVertexPoint(
      nearestVertexIdentifier_[i], proj.x, proj.y, proj.z);
  }

  return proj;
}

int ttk::QuadrangulationSubdivision::project(const size_t firstPointIdx) {
  Timer t;

  nearestVertexIdentifier_.resize(outputPoints_->size());

  // main loop
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = firstPointIdx; i < outputPoints_->size(); i++) {

    nearestVertexIdentifier_[i] = findNearestInputVertex(i);

    // replace curr in outputPoints_ by its projection
    (*outputPoints_)[i] = findProjectionInTriangle(i);
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Projected " << outputPoints_->size() - firstPointIdx
        << " points in " << t.getElapsedTime() << "s" << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}

int ttk::QuadrangulationSubdivision::getQuadNeighbors(
  const bool secondNeighbors) {
  Timer t;

  quadNeighbors_.clear();
  quadNeighbors_.resize(outputPoints_->size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t a = inputVertexNumber_; a < outputPoints_->size(); a++) {
    for(auto &q : *outputQuads_) {
      auto i = static_cast<size_t>(q.i);
      auto j = static_cast<size_t>(q.j);
      auto k = static_cast<size_t>(q.k);
      auto l = static_cast<size_t>(q.l);
      if(secondNeighbors) {
        if(i == a) {
          quadNeighbors_[a].insert(j);
          quadNeighbors_[a].insert(k);
          quadNeighbors_[a].insert(l);
        } else if(j == a) {
          quadNeighbors_[a].insert(i);
          quadNeighbors_[a].insert(k);
          quadNeighbors_[a].insert(l);
        } else if(k == a) {
          quadNeighbors_[a].insert(i);
          quadNeighbors_[a].insert(j);
          quadNeighbors_[a].insert(l);
        } else if(l == a) {
          quadNeighbors_[a].insert(i);
          quadNeighbors_[a].insert(j);
          quadNeighbors_[a].insert(k);
        }
      } else {
        if(i == a || k == a) {
          quadNeighbors_[a].insert(j);
          quadNeighbors_[a].insert(l);
        }
        if(j == a || l == a) {
          quadNeighbors_[a].insert(k);
          quadNeighbors_[a].insert(i);
        }
      }
    }
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Computed neighbors mapping of "
        << outputPoints_->size() - inputVertexNumber_ << " points in "
        << t.getElapsedTime() << "s" << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}

int ttk::QuadrangulationSubdivision::relax() {
  Timer t;

  // loop over output points, do not touch input MSC critical points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = inputVertexNumber_; i < outputPoints_->size(); i++) {
    Point *curr = &(*outputPoints_)[i];

    // barycenter of curr neighbors
    Point relax{};
    for(auto &neigh : quadNeighbors_[i]) {
      relax = relax + (*outputPoints_)[neigh];
    }
    relax = relax * (1.0f / static_cast<float>(quadNeighbors_[i].size()));

    *curr = relax;
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Relaxed " << outputPoints_->size() - inputVertexNumber_
        << " points in " << t.getElapsedTime() << "s" << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}

// main routine
int ttk::QuadrangulationSubdivision::execute() {

  using std::cout;
  using std::endl;

  Timer t;

  outputQuads_->clear();
  outputPoints_->clear();

  // store input points (MSC critical points)
  for(size_t i = 0; i < inputVertexNumber_; i++) {
    outputPoints_->emplace_back(inputVertices_[i]);
  }

  // vector of input quadrangles copied from the inputQuads_ array
  std::vector<Quad> inputQuadsVert;
  // loop variables: pointers to quadrangle vectors
  std::vector<Quad> *tmp0 = &inputQuadsVert, *tmp1 = outputQuads_;
  // holds the subdivision bounds in the outputPoints_ vector
  std::vector<size_t> newPointsRange;
  newPointsRange.emplace_back(outputPoints_->size());

  // copy input quads into vector
  for(size_t i = 0; i < inputQuadNumber_; i++) {
    inputQuadsVert.emplace_back(inputQuads_[i]);
  }

  // main loop
  for(size_t i = 0; i < subdivisionLevel_; i++) {
    // 1. we subdivise each quadrangle by creating five new points, at
    // the center of each edge (4) and at the barycenter of the four
    // vertices (1).

    // index of first point inserted in the outputPoints_ vector
    // during subdivision

    subdivise(*tmp1, *tmp0);
    tmp0->clear();

    auto swap = tmp0;
    tmp0 = tmp1;
    tmp1 = swap;

    // 2. we project every new point on the original 2D mesh, finding
    // the nearest triangle

    project(newPointsRange.back());
    newPointsRange.emplace_back(outputPoints_->size());
  }

  // remainder iteration: since we used outputQuads_ to store
  // temporary quadrangles data, we sometimes need to copy the true
  // output quadrangles into it
  if(subdivisionLevel_ % 2 == 0) {
    outputQuads_->reserve(tmp0->size());
    for(auto &q : *tmp0) {
      outputQuads_->emplace_back(q);
    }
  }

  // retrieve mapping between every vertex and its neighbors
  getQuadNeighbors();

  // 3. we "relax" the new points, i.e. we replace it by the
  // barycenter of its four neighbors
  for(size_t i = 0; i < relaxationIterations_; i++) {
    relax();

    // project all points except MSC critical points
    project(newPointsRange.front());
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Produced " << outputQuads_->size() << " quadrangles with "
        << outputPoints_->size() << " points in " << t.getElapsedTime() << "s ("
        << threadNumber_ << " thread(s))" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
