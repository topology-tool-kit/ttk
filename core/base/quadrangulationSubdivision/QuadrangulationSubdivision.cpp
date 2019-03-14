#include <QuadrangulationSubdivision.h>

ttk::QuadrangulationSubdivision::QuadrangulationSubdivision()
  : vertexNumber_{}, subdivisionLevel_{5}, relaxationIterations_{100},
    criticalPointsNumber_{}, criticalPoints_{}, criticalPointsIdentifiers_{},
    separatriceNumber_{}, sepId_{}, sepSourceId_{}, sepDestId_{} {
}

// main routine
int ttk::QuadrangulationSubdivision::execute() const {

  using std::cout;
  using std::endl;

  Timer t;

  {
    std::stringstream msg;
    msg << "[QuadrangulationSubdivision] Ending computation after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
