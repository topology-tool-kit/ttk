#include <SurfaceQuadrangulation.h>

ttk::SurfaceQuadrangulation::SurfaceQuadrangulation()
  : vertexNumber_{}, subdivisionLevel_{5}, relaxationIterations_{100},
    criticalPointsNumber_{}, criticalPoints_{}, criticalPointsIdentifiers_{},
    separatriceNumber_{}, sepId_{}, sepSourceId_{}, sepDestId_{} {
}

// main routine
int ttk::SurfaceQuadrangulation::execute() const {

  using std::cout;
  using std::endl;
  using std::stringstream;

  Timer t;
  {
    stringstream msg;
    msg << "[SurfaceQuadrangulation] Ending computation after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
