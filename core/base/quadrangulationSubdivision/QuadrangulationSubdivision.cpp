#include <QuadrangulationSubdivision.h>

ttk::QuadrangulationSubdivision::QuadrangulationSubdivision()
  : vertexNumber_{}, subdivisionLevel_{3}, relaxationIterations_{100},
    inputQuads_{}, inputVertexIdentifiers_{}, triangulation_{}, outputQuads_{},
    outputPoints_{} {
}

int ttk::QuadrangulationSubdivision::subdivise(
  std::vector<Quad> &currQuads, const std::vector<Quad> &prevQuads) {

  using edgeType = std::pair<long long, long long>;
  using vertexType = std::pair<long long, Point>;
  std::map<edgeType, vertexType> processedEdges;

  // avoid reallocation in loop, causing invalid pointers
  outputPoints_->reserve(outputPoints_->size() * 5);

  for(auto &q : prevQuads) {
    assert(q.n == 4); // magic number...

    Point *pi = &(*outputPoints_)[q.i];
    Point *pj = &(*outputPoints_)[q.j];
    Point *pk = &(*outputPoints_)[q.k];
    Point *pl = &(*outputPoints_)[q.l];

    // middles of edges
    auto ij = std::make_pair(q.i, q.j);
    auto midij = (*pi + *pj) * 0.5;
    auto jk = std::make_pair(q.j, q.k);
    auto midjk = (*pj + *pk) * 0.5;
    auto kl = std::make_pair(q.k, q.l);
    auto midkl = (*pk + *pl) * 0.5;
    auto li = std::make_pair(q.l, q.i);
    auto midli = (*pl + *pi) * 0.5;

    // quad barycenter
    auto bary = (*pi + *pj + *pk + *pl) * 0.25;

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
    long long baryIdx = outputPoints_->size();
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

  return 0;
}

int ttk::QuadrangulationSubdivision::project(const size_t firstPointIdx) {

  // holds the barycenter coordinate of every triangle of the input mesh
  std::vector<Point> triangleBary(triangulation_->getNumberOfTriangles());

  // compute input triangles barycenters
  for(SimplexId j = 0; j < triangulation_->getNumberOfTriangles(); j++) {

    // get triangle vertices
    SimplexId a, b, c;
    triangulation_->getTriangleVertex(j, 0, a);
    triangulation_->getTriangleVertex(j, 1, b);
    triangulation_->getTriangleVertex(j, 2, c);

    // get coordinates of triangle vertices
    Point pa, pb, pc;
    triangulation_->getVertexPoint(a, pa.x, pa.y, pa.z);
    triangulation_->getVertexPoint(b, pb.x, pb.y, pb.z);
    triangulation_->getVertexPoint(c, pc.x, pc.y, pc.z);

    // get triangle barycenter
    Point bary = (pa + pb + pc) * 1.0 / 3.0;
    triangleBary[j] = bary;
  }

  // main loop
  for(size_t i = firstPointIdx; i < outputPoints_->size(); i++) {

    // current point to project
    Point *curr = &(*outputPoints_)[i];
    // holds the distance to the nearest triangle
    std::pair<float, SimplexId> nearestTriangleDist = std::make_pair(-1.0, -1);

    // iterate over all triangles of the input mesh, find the nearest triangle
    for(SimplexId j = 0; j < triangulation_->getNumberOfTriangles(); j++) {

      // get triangle barycenter
      Point *bary = &triangleBary[j];

      // get distance to triangle barycenter
      float dist = Geometry::distance(&curr->x, &bary->x);

      if(nearestTriangleDist.first < 0.0 || dist < nearestTriangleDist.first) {
        nearestTriangleDist.first = dist;
        nearestTriangleDist.second = j;
      }
    }

    // get triangle vertices
    SimplexId a, b, c;
    triangulation_->getTriangleVertex(nearestTriangleDist.second, 0, a);
    triangulation_->getTriangleVertex(nearestTriangleDist.second, 1, b);
    triangulation_->getTriangleVertex(nearestTriangleDist.second, 2, c);

    // get coordinates of triangle vertices
    Point pa, pb, pc;
    triangulation_->getVertexPoint(a, pa.x, pa.y, pa.z);
    triangulation_->getVertexPoint(b, pb.x, pb.y, pb.z);
    triangulation_->getVertexPoint(c, pc.x, pc.y, pc.z);

    // triangle normal: cross product of two edges
    Point crossP;
    // ab, ac vectors
    Point ab = pb - pa, ac = pc - pa;
    // compute ab ^ ac
    Geometry::crossProduct(&ab.x, &ac.x, &crossP.x);
    // unitary normal vector
    Point norm = crossP / Geometry::magnitude(&crossP.x);

    Point tmp = *curr - pa;
    // projected point into triangle point
    Point proj = *curr - norm * Geometry::dotProduct(&norm.x, &tmp.x);

    // replace curr in outputPoints_ by its projection in the nearest triangle
    *curr = proj;
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

    // 3. we "relax" the new points, i.e. we replace it by the
    // barycenter of its four neighbors

    // we must end by a projection!
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

  {
    std::stringstream msg;
    msg << "[QuadrangulationSubdivision] Produced " << outputQuads_->size()
        << " quadrangles with " << outputPoints_->size() << " points in "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
