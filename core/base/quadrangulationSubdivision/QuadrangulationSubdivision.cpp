#include <Dijkstra.h>
#include <QuadrangulationSubdivision.h>

#define MODULE_S "[QuadrangulationSubdivision] "

ttk::SimplexId
  ttk::QuadrangulationSubdivision::findEdgeMiddle(const size_t a,
                                                  const size_t b) const {

  std::vector<float> sum(vertexDistance_[a].size());

  // euclidian barycenter of a and b
  Point edgeEuclBary = ((*outputPoints_)[a] + (*outputPoints_)[b]) * 0.5F;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sum.size(); ++i) {
    float m = vertexDistance_[a][i];
    float n = vertexDistance_[b][i];
    // stay on the shortest path between a and b
    sum[i] = m + n;
    if(m != std::numeric_limits<float>::infinity()
       && n != std::numeric_limits<float>::infinity()) {
      // try to get the middle of the shortest path
      sum[i] += std::abs(m - n);
    }

    // get the euclidian distance to AB
    Point curr{};
    triangulation_->getVertexPoint(i, curr.x, curr.y, curr.z);
    // try to minimize the euclidian distance to AB too
    sum[i] += Geometry::distance(&curr.x, &edgeEuclBary.x);
  }

  return std::min_element(sum.begin(), sum.end()) - sum.begin();
}

ttk::SimplexId ttk::QuadrangulationSubdivision::findQuadBary(
  const std::vector<size_t> &quadVertices) const {

  std::vector<float> sum(vertexDistance_[*quadVertices.begin()].size(),
                         std::numeric_limits<float>::infinity());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sum.size(); ++i) {

    // skip following computation if too far from any parent quad vertex
    bool skip = std::any_of(
      quadVertices.begin(), quadVertices.end(), [&](const size_t &a) {
        return vertexDistance_[a][i] == std::numeric_limits<float>::infinity();
      });

    if(skip) {
      continue;
    }

    float m = vertexDistance_[quadVertices[0]][i];
    float n = vertexDistance_[quadVertices[1]][i];
    float o = vertexDistance_[quadVertices[2]][i];
    float p = vertexDistance_[quadVertices[3]][i];

    // try to be "near" the four parent vertices
    sum[i] = m + n + o + p;

    // try to be on the diagonals intersection
    sum[i] += std::abs(m - o);
    sum[i] += std::abs(n - p);
  }

  return std::min_element(sum.begin(), sum.end()) - sum.begin();
}

int ttk::QuadrangulationSubdivision::subdivise() {

  using edgeType = std::pair<long long, long long>;
  using vertexType = std::pair<long long, Point>;
  using std::make_pair;
  std::map<edgeType, vertexType> processedEdges;

  // deep copy of coarse input quads
  auto prevQuads(*outputQuads_);

  // clear input quads buffer before re-writing it
  outputQuads_->clear();

  Timer t;

  // avoid reallocation in loop, causing invalid pointers
  outputPoints_->reserve(outputPoints_->size() * 5);

  vertexDistance_.resize(outputPoints_->size());

  // get all other vertices sharing a quad
  getQuadNeighbors(prevQuads, quadNeighbors_, true);

  // compute shortest distance from every vertex to all other that share a quad
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_->size(); ++i) {

    // skip if already computed on a coarser subdivision
    if(vertexDistance_[i].empty()) {

      // do not propagate on the whole mesh
      std::vector<SimplexId> bounds;
      for(auto &p : quadNeighbors_[i]) {
        bounds.emplace_back(nearestVertexIdentifier_[p]);
      }

      Dijkstra::shortestPath(nearestVertexIdentifier_[i], *triangulation_,
                             vertexDistance_[i], bounds);
    }
  }

  for(auto &q : prevQuads) {
    assert(q.n == 4); // magic number...

    auto i = static_cast<size_t>(q.i);
    auto j = static_cast<size_t>(q.j);
    auto k = static_cast<size_t>(q.k);
    auto l = static_cast<size_t>(q.l);

    // middles of edges
    auto ijid = findEdgeMiddle(i, j);
    auto jkid = findEdgeMiddle(j, k);
    auto klid = findEdgeMiddle(k, l);
    auto liid = findEdgeMiddle(l, i);

    Point midij{};
    triangulation_->getVertexPoint(ijid, midij.x, midij.y, midij.z);
    Point midjk{};
    triangulation_->getVertexPoint(jkid, midjk.x, midjk.y, midjk.z);
    Point midkl{};
    triangulation_->getVertexPoint(klid, midkl.x, midkl.y, midkl.z);
    Point midli{};
    triangulation_->getVertexPoint(liid, midli.x, midli.y, midli.z);

    std::vector<size_t> quadVertices{i, j, k, l};
    // barycenter TTK identifier
    auto baryid = findQuadBary(quadVertices);
    // barycenter 3D coordinates
    Point bary{};
    triangulation_->getVertexPoint(baryid, bary.x, bary.y, bary.z);

    // order edges to avoid duplicates (ij vs. ji)
    auto ij = make_pair(std::min(q.i, q.j), std::max(q.i, q.j));
    auto jk = make_pair(std::min(q.j, q.k), std::max(q.j, q.k));
    auto kl = make_pair(std::min(q.k, q.l), std::max(q.k, q.l));
    auto li = make_pair(std::min(q.l, q.i), std::max(q.l, q.i));

#define PROCESS_EDGE_MIDDLE(PARENT_PAIR, POINT_3D, TTK_ID)                 \
  /* check if edge already processed by a neighbor quad */                 \
  if(processedEdges.find(PARENT_PAIR) == processedEdges.end()) {           \
    processedEdges.insert(                                                 \
      make_pair(PARENT_PAIR, make_pair(outputPoints_->size(), POINT_3D))); \
    /* add new point 3d coordinates to vector of output points */          \
    outputPoints_->emplace_back(POINT_3D);                                 \
    /* new point is an edge middle */                                      \
    outputVertType_->emplace_back(1);                                      \
    /* store also TTK identifier of triangular mesh vertex */              \
    nearestVertexIdentifier_.emplace_back(TTK_ID);                         \
  }

    PROCESS_EDGE_MIDDLE(ij, midij, ijid);
    PROCESS_EDGE_MIDDLE(jk, midjk, jkid);
    PROCESS_EDGE_MIDDLE(kl, midkl, klid);
    PROCESS_EDGE_MIDDLE(li, midli, liid);

    // barycenter index in outputPoints_
    auto baryIdx = static_cast<long long>(outputPoints_->size());
    outputPoints_->emplace_back(bary);
    outputVertType_->emplace_back(2);
    nearestVertexIdentifier_.emplace_back(baryid);

    // add the four new quads
    outputQuads_->emplace_back(Quad{
      4, q.i, processedEdges[ij].first, baryIdx, processedEdges[li].first});
    outputQuads_->emplace_back(Quad{
      4, q.j, processedEdges[jk].first, baryIdx, processedEdges[ij].first});
    outputQuads_->emplace_back(Quad{
      4, q.k, processedEdges[kl].first, baryIdx, processedEdges[jk].first});
    outputQuads_->emplace_back(Quad{
      4, q.l, processedEdges[li].first, baryIdx, processedEdges[kl].first});
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Subdivised " << prevQuads.size() << " quads into "
        << outputQuads_->size() << " new quads (" << outputPoints_->size()
        << " points) in " << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}

ttk::QuadrangulationSubdivision::Point
  ttk::QuadrangulationSubdivision::findProjectionInTriangle(
    const ttk::SimplexId i, const bool lastIter) {

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
    std::array<SimplexId, 3> tverts;
    triangulation_->getTriangleVertex(tid, 0, tverts[0]);
    triangulation_->getTriangleVertex(tid, 1, tverts[1]);
    triangulation_->getTriangleVertex(tid, 2, tverts[2]);

    // get coordinates of triangle vertices
    Point pa{};
    triangulation_->getVertexPoint(tverts[0], pa.x, pa.y, pa.z);
    Point pb{};
    triangulation_->getVertexPoint(tverts[1], pb.x, pb.y, pb.z);
    Point pc{};
    triangulation_->getVertexPoint(tverts[2], pc.x, pc.y, pc.z);

    // triangle normal: cross product of two edges
    Point crossP{};
    // ab, ac vectors
    Point ab = pb - pa;
    Point ac = pc - pa;
    // compute ab ^ ac
    Geometry::crossProduct(&ab.x, &ac.x, &crossP.x);
    // unitary normal vector
    Point norm = crossP / Geometry::magnitude(&crossP.x);

    Point tmp = *vert - pa;
    // projected point into triangle
    proj = *vert - norm * Geometry::dotProduct(&norm.x, &tmp.x);

    // compute barycentric coords of projection
    std::vector<float> baryCoords;
    Geometry::computeBarycentricCoordinates(
      &pa.x, &pb.x, &pc.x, &proj.x, baryCoords);

    // check if projection in triangle
    bool inTriangle = true;
    for(auto &coord : baryCoords) {
      if(coord < powf(10, -FLT_DIG)) {
        inTriangle = false;
      }
      if(coord > 1 + powf(10, -FLT_DIG)) {
        inTriangle = false;
      }
    }

    if(inTriangle) {
      success = true;
      // should we check if we have the nearest triangle?
      break;
    }

    // mark triangle as tested
    trianglesTested[tid] = true;

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
    for(auto &avert : vertices) {
      SimplexId tnum = triangulation_->getVertexTriangleNumber(avert);
      for(SimplexId j = 0; j < tnum; j++) {
        SimplexId trid;
        triangulation_->getVertexTriangle(avert, j, trid);
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

  // fill in debug info
  if(lastIter) {
    (*trianglesChecked_)[i]
      = std::count(trianglesTested.begin(), trianglesTested.end(), true);
    (*projSucceeded_)[i] = success ? 1 : 0;
  }

  return proj;
}

int ttk::QuadrangulationSubdivision::project(const std::set<size_t> &filtered,
                                             const bool lastIter) {
  Timer t;

  if(lastIter) {
    trianglesChecked_->clear();
    projSucceeded_->clear();
    trianglesChecked_->resize(outputPoints_->size());
    projSucceeded_->resize(outputPoints_->size());
  }

  // compute nearest vertex in triangular mesh and stores it in
  // nearestVertexIdentifier_
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_->size(); i++) {

    // skip computation if i in filtered
    if(filtered.find(i) != filtered.end()) {
      continue;
    }

    const Point *vert = &(*outputPoints_)[i];
    // distance to triangular mesh vertex
    float min_dist = std::numeric_limits<float>::infinity();
    // iterate over the whole triangular mesh
    for(SimplexId j = 0; j < vertexNumber_; ++j) {
      Point p{};
      triangulation_->getVertexPoint(j, p.x, p.y, p.z);
      float curr_dist = Geometry::distance(&vert->x, &p.x);
      if(curr_dist < min_dist) {
        min_dist = curr_dist;
        nearestVertexIdentifier_[i] = j;
      }
    }
  }

  // main loop
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_->size(); i++) {

    // skip computation if i in filtered
    if(filtered.find(i) != filtered.end()) {
      continue;
    }

    // replace curr in outputPoints_ by its projection
    (*outputPoints_)[i] = findProjectionInTriangle(i, lastIter);
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Projected " << outputPoints_->size() - filtered.size()
        << " points in " << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}

int ttk::QuadrangulationSubdivision::getQuadNeighbors(
  const std::vector<Quad> &quads,
  std::vector<std::set<size_t>> &neighbors,
  const bool secondNeighbors) {
  Timer t;

  quadNeighbors_.clear();
  quadNeighbors_.resize(outputPoints_->size());

  for(auto &q : quads) {
    auto i = static_cast<size_t>(q.i);
    auto j = static_cast<size_t>(q.j);
    auto k = static_cast<size_t>(q.k);
    auto l = static_cast<size_t>(q.l);
    if(secondNeighbors) {
      neighbors[i].insert(j);
      neighbors[i].insert(k);
      neighbors[i].insert(l);
      neighbors[j].insert(i);
      neighbors[j].insert(k);
      neighbors[j].insert(l);
      neighbors[k].insert(i);
      neighbors[k].insert(j);
      neighbors[k].insert(l);
      neighbors[l].insert(i);
      neighbors[l].insert(j);
      neighbors[l].insert(k);
    } else {
      neighbors[i].insert(j);
      neighbors[i].insert(l);
      neighbors[k].insert(j);
      neighbors[k].insert(l);
      neighbors[j].insert(k);
      neighbors[j].insert(i);
      neighbors[l].insert(k);
      neighbors[l].insert(i);
    }
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Computed neighbors mapping of " << outputPoints_->size()
        << " points in " << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}

int ttk::QuadrangulationSubdivision::relax(const std::set<size_t> &filtered) {
  Timer t;

  // outputPoints_ deep copy
  auto tmp(*outputPoints_);

  // loop over output points, do not touch input MSC critical points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_->size(); i++) {

    // skip computation if i in filtered
    if(filtered.find(i) != filtered.end()) {
      continue;
    }

    // barycenter of curr neighbors
    Point relax{};
    for(auto &neigh : quadNeighbors_[i]) {
      relax = relax + tmp[neigh];
    }
    relax = relax * (1.0F / static_cast<float>(quadNeighbors_[i].size()));

    (*outputPoints_)[i] = relax;
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Relaxed " << outputPoints_->size() - inputVertexNumber_
        << " points in " << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}

int ttk::QuadrangulationSubdivision::findExtraordinaryVertices(
  std::set<size_t> &output) {

  // clear output
  output.clear();

  // hold input quads in a vector
  std::vector<Quad> inputQuads(inputQuadNumber_);

  std::vector<std::set<size_t>> neighbors(inputVertexNumber_);

  // use outputQuads_ here because it contains the input quadrangles before
  // subdivision wrt function call position in execute()
  getQuadNeighbors(*outputQuads_, neighbors);

  const size_t NORMAL_VALENCE = 4;

  for(size_t i = 0; i < neighbors.size(); ++i) {
    if(neighbors[i].size() != NORMAL_VALENCE) {
      output.insert(i);
    }
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
  outputVertType_->clear();

  // store input points (MSC critical points)
  for(size_t i = 0; i < inputVertexNumber_; i++) {
    outputPoints_->emplace_back(inputVertices_[i]);
  }

  // copy input quads into vector
  for(size_t i = 0; i < inputQuadNumber_; i++) {
    outputQuads_->emplace_back(inputQuads_[i]);
  }

  // fill outputInfos_ with input data (critical points)
  outputVertType_->resize(outputPoints_->size());
  std::fill(outputVertType_->begin(), outputVertType_->end(), 0);

  // vertices to filter from relaxation, projection
  std::set<size_t> filtered{};
  if(!lockAllInputVertices) {
    if(lockInputExtrema) {
      // get extraordinary vertices
      findExtraordinaryVertices(filtered);
    }
  } else {
    // fill vector with all input points indices from 0 to inputVertexNumber_
    for(size_t i = 0; i < inputVertexNumber_; ++i) {
      filtered.insert(i);
    }
  }

  // main loop
  for(size_t i = 0; i < subdivisionLevel_; i++) {
    // subdivise each quadrangle by creating five new points, at the
    // center of each edge (4) and at the barycenter of the four
    // vertices (1).
    subdivise();
  }

  // retrieve mapping between every vertex and its neighbors
  getQuadNeighbors(*outputQuads_, quadNeighbors_);

  // "relax" the new points, i.e. replace it by the barycenter of its
  // four neighbors
  for(size_t i = 0; i < relaxationIterations_; i++) {
    relax(filtered);

    // project all points on the nearest triangle (except MSC critical
    // points)
    project(filtered, (i == relaxationIterations_ - 1));
  }

  // compute valence of every quadrangle vertex
  outputValences_->resize(outputPoints_->size());
  std::transform(
    quadNeighbors_.begin(), quadNeighbors_.end(), outputValences_->begin(),
    [&](const std::set<size_t> &neighbors) { return neighbors.size(); });

  {
    std::stringstream msg;
    msg << MODULE_S "Produced " << outputQuads_->size() << " quadrangles with "
        << outputPoints_->size() << " points in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
