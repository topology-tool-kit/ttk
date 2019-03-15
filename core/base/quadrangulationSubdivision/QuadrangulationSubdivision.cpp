#include <QuadrangulationSubdivision.h>

ttk::QuadrangulationSubdivision::QuadrangulationSubdivision()
  : vertexNumber_{}, subdivisionLevel_{5}, relaxationIterations_{100},
    inputQuadVertexNumber_{}, inputQuadrangles_{}, triangulation_{},
    outputQuads_{} {
}

// main routine
int ttk::QuadrangulationSubdivision::execute() const {

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
