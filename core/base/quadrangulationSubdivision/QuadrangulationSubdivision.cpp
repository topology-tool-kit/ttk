#include <QuadrangulationSubdivision.h>

ttk::QuadrangulationSubdivision::QuadrangulationSubdivision()
  : vertexNumber_{}, subdivisionLevel_{5}, relaxationIterations_{100},
    inputQuads_{}, inputVertexIdentifiers_{}, triangulation_{}, outputQuads_{},
    outputPoints_{} {
}

int ttk::QuadrangulationSubdivision::subdivise() {

  using edgeType = std::pair<long long, long long>;
  std::map<edgeType, ttk::QuadrangulationSubdivision::Point> processedEdges;

  // five values per quad
  for(size_t i = 0; i < inputQuadNumber_; i++) {

    auto q = inputQuads_[i];
    assert(q.n == 4); // magic number...

    Point quadVertices[4];

    // get current quad 3D coordinates
    triangulation_->getVertexPoint(inputVertexIdentifiers_[q.i],
                                   quadVertices[0].x, quadVertices[0].y,
                                   quadVertices[0].z);
    triangulation_->getVertexPoint(inputVertexIdentifiers_[q.j],
                                   quadVertices[1].x, quadVertices[1].y,
                                   quadVertices[1].z);
    triangulation_->getVertexPoint(inputVertexIdentifiers_[q.k],
                                   quadVertices[2].x, quadVertices[2].y,
                                   quadVertices[2].z);
    triangulation_->getVertexPoint(inputVertexIdentifiers_[q.l],
                                   quadVertices[3].x, quadVertices[3].y,
                                   quadVertices[3].z);

    // middles of edges
    auto ij = std::make_pair(q.i, q.j);
    auto midij = (quadVertices[0] + quadVertices[1]) * 0.5;
    if(processedEdges.find(ij) == processedEdges.end()) {
      processedEdges.insert(std::make_pair(ij, midij));
      outputPoints_->emplace_back(midij);
    }

    auto jk = std::make_pair(q.j, q.k);
    auto midjk = (quadVertices[1] + quadVertices[2]) * 0.5;
    if(processedEdges.find(jk) == processedEdges.end()) {
      processedEdges.insert(std::make_pair(jk, midjk));
      outputPoints_->emplace_back(midjk);
    }

    auto kl = std::make_pair(q.k, q.l);
    auto midkl = (quadVertices[2] + quadVertices[3]) * 0.5;
    if(processedEdges.find(kl) == processedEdges.end()) {
      processedEdges.insert(std::make_pair(kl, midkl));
      outputPoints_->emplace_back(midkl);
    }

    auto li = std::make_pair(q.l, q.i);
    auto midli = (quadVertices[3] + quadVertices[0]) * 0.5;
    if(processedEdges.find(li) == processedEdges.end()) {
      processedEdges.insert(std::make_pair(li, midli));
      outputPoints_->emplace_back(midli);
    }

    // quad barycenter
    auto bary
      = (quadVertices[0] + quadVertices[1] + quadVertices[2] + quadVertices[3])
        * 0.25;
    outputPoints_->emplace_back(bary);

    // add the four new quads
    // TODO
  }

  return 0;
}

// main routine
int ttk::QuadrangulationSubdivision::execute() {

  using std::cout;
  using std::endl;

  Timer t;

  // main loop
  for(size_t i = 0; i < subdivisionLevel_; i++) {
    // 1. we subdivise each quadrangle by creating five new points, at
    // the center of each edge (4) and at the barycenter of the four
    // vertices (1).

    // 2. we project every new point on the original 2D mesh, finding
    // the nearest triangle

    // 3. we "relax" the new points, i.e. we replace it by the
    // barycenter of its four neighbors

    // we must end by a projection!
  }

  subdivise();

  {
    std::stringstream msg;
    msg << "[QuadrangulationSubdivision] Ending computation after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
