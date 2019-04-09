#include <Dijkstra.h>
#include <QuadrangulationSubdivision.h>

#define MODULE_S "[QuadrangulationSubdivision] "

static ttk::SimplexId minVecAdd(const std::vector<float> vec0,
                                const std::vector<float> vec1) {
  std::vector<float> sum(vec0.size());

  std::transform(
    vec0.begin(), vec0.end(), vec1.begin(), sum.begin(), std::plus<float>());

  return std::min_element(sum.begin(), sum.end()) - sum.begin();
}

static ttk::SimplexId minVecAdd4(const std::vector<float> vec0,
                                 const std::vector<float> vec1,
                                 const std::vector<float> vec2,
                                 const std::vector<float> vec3) {
  std::vector<float> sum(vec0.size());

  std::transform(
    vec0.begin(), vec0.end(), vec1.begin(), sum.begin(), std::plus<float>());
  std::transform(
    vec2.begin(), vec2.end(), sum.begin(), sum.begin(), std::plus<float>());
  std::transform(
    vec3.begin(), vec3.end(), sum.begin(), sum.begin(), std::plus<float>());

  return std::min_element(sum.begin(), sum.end()) - sum.begin();
}

int ttk::QuadrangulationSubdivision::subdivise(
  std::vector<Quad> &currQuads, const std::vector<Quad> &prevQuads) {

  using edgeType = std::pair<long long, long long>;
  using vertexType = std::pair<long long, Point>;
  using std::make_pair;
  std::map<edgeType, vertexType> processedEdges;

  Timer t;

  // avoid reallocation in loop, causing invalid pointers
  outputPoints_->reserve(outputPoints_->size() * 5);

  vertexDistance_.resize(outputPoints_->size());

  // get all other vertices sharing a quad
  getQuadNeighbors(prevQuads, true);

  // compute shortest distance from every vertex to all other that share a quad
  for(size_t i = 0; i < outputPoints_->size(); ++i) {

    // skip if already computed on a coarser subdivision
    if(vertexDistance_[i].empty()) {

      // do not propagate on the whole mesh
      std::vector<SimplexId> bounds;
      for(auto &p : quadNeighbors_[i]) {
        bounds.emplace_back(nearestVertexIdentifier_[p]);
      }

      Dijkstra::shortestPath(nearestVertexIdentifier_[i], *triangulation_,
                             vertexDistance_[i], &bounds);
    }
  }

  for(auto &q : prevQuads) {
    assert(q.n == 4); // magic number...

    auto i = static_cast<size_t>(q.i);
    auto j = static_cast<size_t>(q.j);
    auto k = static_cast<size_t>(q.k);
    auto l = static_cast<size_t>(q.l);

    // middles of edges
    auto ijid = minVecAdd(vertexDistance_[i], vertexDistance_[j]);
    auto jkid = minVecAdd(vertexDistance_[j], vertexDistance_[k]);
    auto klid = minVecAdd(vertexDistance_[k], vertexDistance_[l]);
    auto liid = minVecAdd(vertexDistance_[l], vertexDistance_[i]);

    Point midij;
    triangulation_->getVertexPoint(ijid, midij.x, midij.y, midij.z);
    Point midjk;
    triangulation_->getVertexPoint(jkid, midjk.x, midjk.y, midjk.z);
    Point midkl;
    triangulation_->getVertexPoint(klid, midkl.x, midkl.y, midkl.z);
    Point midli;
    triangulation_->getVertexPoint(liid, midli.x, midli.y, midli.z);

    // quad barycenter
    auto baryid = minVecAdd4(vertexDistance_[i], vertexDistance_[j],
                             vertexDistance_[k], vertexDistance_[l]);
    Point bary;
    triangulation_->getVertexPoint(baryid, bary.x, bary.y, bary.z);

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
      nearestVertexIdentifier_.emplace_back(ijid);
    }

    if(processedEdges.find(jk) == processedEdges.end()) {
      processedEdges.insert(
        make_pair(jk, make_pair(outputPoints_->size(), midjk)));
      outputPoints_->emplace_back(midjk);
      nearestVertexIdentifier_.emplace_back(jkid);
    }

    if(processedEdges.find(kl) == processedEdges.end()) {
      processedEdges.insert(
        make_pair(kl, make_pair(outputPoints_->size(), midkl)));
      outputPoints_->emplace_back(midkl);
      nearestVertexIdentifier_.emplace_back(klid);
    }

    if(processedEdges.find(li) == processedEdges.end()) {
      processedEdges.insert(
        make_pair(li, make_pair(outputPoints_->size(), midli)));
      outputPoints_->emplace_back(midli);
      nearestVertexIdentifier_.emplace_back(liid);
    }

    // barycenter index in outputPoints_
    auto baryIdx = static_cast<long long>(outputPoints_->size());
    outputPoints_->emplace_back(bary);
    nearestVertexIdentifier_.emplace_back(baryid);

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
  const std::vector<Quad> &quads, const bool secondNeighbors) {
  Timer t;

  quadNeighbors_.clear();
  quadNeighbors_.resize(outputPoints_->size());

  for(auto &q : quads) {
    auto i = static_cast<size_t>(q.i);
    auto j = static_cast<size_t>(q.j);
    auto k = static_cast<size_t>(q.k);
    auto l = static_cast<size_t>(q.l);
    if(secondNeighbors) {
      quadNeighbors_[i].insert(j);
      quadNeighbors_[i].insert(k);
      quadNeighbors_[i].insert(l);
      quadNeighbors_[j].insert(i);
      quadNeighbors_[j].insert(k);
      quadNeighbors_[j].insert(l);
      quadNeighbors_[k].insert(i);
      quadNeighbors_[k].insert(j);
      quadNeighbors_[k].insert(l);
      quadNeighbors_[l].insert(i);
      quadNeighbors_[l].insert(j);
      quadNeighbors_[l].insert(k);
    } else {
      quadNeighbors_[i].insert(j);
      quadNeighbors_[i].insert(l);
      quadNeighbors_[k].insert(j);
      quadNeighbors_[k].insert(l);
      quadNeighbors_[j].insert(k);
      quadNeighbors_[j].insert(i);
      quadNeighbors_[l].insert(k);
      quadNeighbors_[l].insert(i);
    }
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Computed neighbors mapping of " << outputPoints_->size()
        << " points in " << t.getElapsedTime() << "s" << std::endl;
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
  getQuadNeighbors(*outputQuads_);

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
