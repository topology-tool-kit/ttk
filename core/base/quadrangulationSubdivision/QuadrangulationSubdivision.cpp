#include <QuadrangulationSubdivision.h>

ttk::QuadrangulationSubdivision::QuadrangulationSubdivision()
  : vertexNumber_{}, subdivisionLevel_{5}, relaxationIterations_{100},
    inputQuads_{}, inputVertexIdentifiers_{}, triangulation_{}, outputQuads_{},
    outputPoints_{} {
}

int ttk::QuadrangulationSubdivision::subdivise(
  std::vector<Quad> &currQuads, const std::vector<Quad> &prevQuads) {

  using edgeType = std::pair<long long, long long>;
  using vertexType = std::pair<long long, Point>;
  std::map<edgeType, vertexType> processedEdges;

  for(auto &q : prevQuads) {
    assert(q.n == 4); // magic number...

    Point *pi = &(*outputPoints_)[q.i];
    Point *pj = &(*outputPoints_)[q.j];
    Point *pk = &(*outputPoints_)[q.k];
    Point *pl = &(*outputPoints_)[q.l];

    // middles of edges
    auto ij = std::make_pair(q.i, q.j);
    auto midij = (*pi + *pj) * 0.5;
    if(processedEdges.find(ij) == processedEdges.end()) {
      processedEdges.insert(
        std::make_pair(ij, std::make_pair(outputPoints_->size(), midij)));
      outputPoints_->emplace_back(midij);
    }

    auto jk = std::make_pair(q.j, q.k);
    auto midjk = (*pj + *pk) * 0.5;
    if(processedEdges.find(jk) == processedEdges.end()) {
      processedEdges.insert(
        std::make_pair(jk, std::make_pair(outputPoints_->size(), midjk)));
      outputPoints_->emplace_back(midjk);
    }

    auto kl = std::make_pair(q.k, q.l);
    auto midkl = (*pk + *pl) * 0.5;
    if(processedEdges.find(kl) == processedEdges.end()) {
      processedEdges.insert(
        std::make_pair(kl, std::make_pair(outputPoints_->size(), midkl)));
      outputPoints_->emplace_back(midkl);
    }

    auto li = std::make_pair(q.l, q.i);
    auto midli = (*pl + *pi) * 0.5;
    if(processedEdges.find(li) == processedEdges.end()) {
      processedEdges.insert(
        std::make_pair(li, std::make_pair(outputPoints_->size(), midli)));
      outputPoints_->emplace_back(midli);
    }

    // quad barycenter
    auto bary = (*pi + *pj + *pk + *pl) * 0.25;
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

  std::vector<Quad> tmp0, tmp1;

  // copy input quads into vector
  for(size_t i = 0; i < inputQuadNumber_; i++) {
    tmp0.emplace_back(inputQuads_[i]);
  }

  // main loop
  for(size_t i = 0; i < subdivisionLevel_ / 2; i++) {
    // 1. we subdivise each quadrangle by creating five new points, at
    // the center of each edge (4) and at the barycenter of the four
    // vertices (1).

    subdivise(tmp1, tmp0);
    tmp0.clear();
    subdivise(tmp0, tmp1);
    tmp1.clear();

    // 2. we project every new point on the original 2D mesh, finding
    // the nearest triangle

    // 3. we "relax" the new points, i.e. we replace it by the
    // barycenter of its four neighbors

    // we must end by a projection!
  }

  switch(subdivisionLevel_ % 2) {
    case 0:
      outputQuads_->reserve(tmp0.size());
      for(auto &q : tmp0) {
        outputQuads_->emplace_back(q);
      }
      break;
    case 1:
      subdivise(*outputQuads_, tmp0);
      break;
  }

  {
    std::stringstream msg;
    msg << "[QuadrangulationSubdivision] Ending computation after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
