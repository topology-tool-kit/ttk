#include <SurfaceQuadrangulation.h>

ttk::SurfaceQuadrangulation::SurfaceQuadrangulation()
  : vertexNumber_{}, triangulation_{}, subdivisionLevel_{5},
    relaxationIterations_{100}, inputScalarFieldPointer_{},
    inputOffsetIdentifiersFieldPointer_{} {
}

// main routine
int ttk::SurfaceQuadrangulation::execute() const {

  using std::cout;
  using std::endl;
  using std::stringstream;

  Timer t;

  MorseSmaleComplex msc{};
  SimplexId criticalPointsNumber{};
  std::vector<float> criticalPoints;
  std::vector<char> cpCellDims;
  std::vector<SimplexId> cpCellIds;
  std::vector<float> cpCellScalars;
  std::vector<char> cpIsOnBoundary;
  std::vector<SimplexId> cpPLVertexIdentifiers;
  std::vector<SimplexId> cpManifoldSize;

  msc.setupTriangulation(triangulation_);
  msc.setInputScalarField(inputScalarFieldPointer_);
  msc.setInputOffsets(inputOffsetIdentifiersFieldPointer_);
  msc.setOutputCriticalPoints(
    &criticalPointsNumber, &criticalPoints, &cpCellDims, &cpCellIds,
    &cpCellScalars, &cpIsOnBoundary, &cpPLVertexIdentifiers, &cpManifoldSize);

  {
    stringstream msg;
    msg << "[SurfaceQuadrangulation] Ending computation after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
