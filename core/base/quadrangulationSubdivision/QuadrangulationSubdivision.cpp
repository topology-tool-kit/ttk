#include <QuadrangulationSubdivision.h>

#define MODULE_S "[QuadrangulationSubdivision] "

int ttk::QuadrangulationSubdivision::subdivise(
  std::vector<Quad> &currQuads, const std::vector<Quad> &prevQuads) {

  using edgeType = std::pair<long long, long long>;
  using vertexType = std::pair<long long, Point>;
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
    auto ij = std::make_pair(std::min(q.i, q.j), std::max(q.i, q.j));
    auto jk = std::make_pair(std::min(q.j, q.k), std::max(q.j, q.k));
    auto kl = std::make_pair(std::min(q.k, q.l), std::max(q.k, q.l));
    auto li = std::make_pair(std::min(q.l, q.i), std::max(q.l, q.i));

    // add to outputPoints_ after computing new point coordinates to
    // avoid invalidating pointers
    if(processedEdges.find(ij) == processedEdges.end()) {
      processedEdges.insert(
        std::make_pair(ij, std::make_pair(outputPoints_->size(), midij)));
      outputPoints_->emplace_back(midij);
    }

    if(processedEdges.find(jk) == processedEdges.end()) {
      processedEdges.insert(
        std::make_pair(jk, std::make_pair(outputPoints_->size(), midjk)));
      outputPoints_->emplace_back(midjk);
    }

    if(processedEdges.find(kl) == processedEdges.end()) {
      processedEdges.insert(
        std::make_pair(kl, std::make_pair(outputPoints_->size(), midkl)));
      outputPoints_->emplace_back(midkl);
    }

    if(processedEdges.find(li) == processedEdges.end()) {
      processedEdges.insert(
        std::make_pair(li, std::make_pair(outputPoints_->size(), midli)));
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

int ttk::QuadrangulationSubdivision::project(const size_t firstPointIdx) {
  Timer t;

  // main loop
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = firstPointIdx; i < outputPoints_->size(); i++) {

    // current point to project
    Point *curr = &(*outputPoints_)[i];
    // holds the distance to the nearest vertex
    std::pair<float, SimplexId> nearestVertex = std::make_pair(-1.0f, -1);

    // iterate over all vertices of the input mesh, find the nearest one
    for(SimplexId j = 0; j < triangulation_->getNumberOfVertices(); j++) {

      // get vertex coordinates
      Point inMesh{};
      triangulation_->getVertexPoint(j, inMesh.x, inMesh.y, inMesh.z);

      // get square distance to vertex
      float dist = (curr->x - inMesh.x) * (curr->x - inMesh.x)
                   + (curr->y - inMesh.y) * (curr->y - inMesh.y)
                   + (curr->z - inMesh.z) * (curr->z - inMesh.z);

      if(nearestVertex.first < 0.0f || dist < nearestVertex.first) {
        nearestVertex.first = dist;
        nearestVertex.second = j;
      }
    }

    // projected point into triangle
    Point proj{};
    // found a projection in one triangle
    bool success = false;
    // number of triangles around nearest vertex
    SimplexId triangleNumber
      = triangulation_->getVertexTriangleNumber(nearestVertex.second);

    // iterate over nearest vertex triangles, find a projection
    for(SimplexId j = 0; j < triangleNumber; j++) {
      SimplexId tid;
      triangulation_->getVertexTriangle(nearestVertex.second, j, tid);

      // get triangle vertices
      SimplexId a, b, c;
      triangulation_->getTriangleVertex(tid, 0, a);
      triangulation_->getTriangleVertex(tid, 1, b);
      triangulation_->getTriangleVertex(tid, 2, c);

      // get coordinates of triangle vertices
      Point pa{}, pb{}, pc{};
      triangulation_->getVertexPoint(a, pa.x, pa.y, pa.z);
      triangulation_->getVertexPoint(b, pb.x, pb.y, pb.z);
      triangulation_->getVertexPoint(c, pc.x, pc.y, pc.z);

      // triangle normal: cross product of two edges
      Point crossP{};
      // ab, ac vectors
      Point ab = pb - pa, ac = pc - pa;
      // compute ab ^ ac
      Geometry::crossProduct(&ab.x, &ac.x, &crossP.x);
      // unitary normal vector
      Point norm = crossP / Geometry::magnitude(&crossP.x);

      Point tmp = *curr - pa;
      // projected point into triangle point
      proj = *curr - norm * Geometry::dotProduct(&norm.x, &tmp.x);

      if(Geometry::isPointInTriangle(&pa.x, &pb.x, &pc.x, &proj.x)) {
        success = true;
        break;
      }
    }

    if(!success) {
      // replace proj by the nearest vertex
      triangulation_->getVertexPoint(
        nearestVertex.second, proj.x, proj.y, proj.z);
    }

    // replace curr in outputPoints_ by its projection
    *curr = proj;
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Projected " << outputPoints_->size() - firstPointIdx
        << " points in " << t.getElapsedTime() << "s" << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}

int ttk::QuadrangulationSubdivision::relax() {
  Timer t;

  // maps every vertex to its quad neighbors
  std::vector<std::set<size_t>> quadNeighbors(outputPoints_->size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t a = inputVertexNumber_; a < outputPoints_->size(); a++) {
    for(auto &q : *outputQuads_) {
      auto i = static_cast<size_t>(q.i);
      auto j = static_cast<size_t>(q.j);
      auto k = static_cast<size_t>(q.k);
      auto l = static_cast<size_t>(q.l);
      if(i == a || k == a) {
        quadNeighbors[a].insert(j);
        quadNeighbors[a].insert(l);
      }
      if(j == a || l == a) {
        quadNeighbors[a].insert(k);
        quadNeighbors[a].insert(i);
      }
    }
  }

  // loop over output points, do not touch input MSC critical points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = inputVertexNumber_; i < outputPoints_->size(); i++) {
    Point *curr = &(*outputPoints_)[i];

    // barycenter of curr neighbors
    Point relax{};
    for(auto &neigh : quadNeighbors[i]) {
      relax = relax + (*outputPoints_)[neigh];
    }
    relax = relax * (1.0f / static_cast<float>(quadNeighbors[i].size()));

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
