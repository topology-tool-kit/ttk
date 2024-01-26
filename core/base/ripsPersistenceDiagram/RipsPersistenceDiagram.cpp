#include <RipsPersistenceDiagram.h>

ttk::RipsPersistenceDiagram::RipsPersistenceDiagram() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("RipsPersistenceDiagram");
}

int ttk::RipsPersistenceDiagram::execute(
  const std::vector<std::vector<double>> &points,
  std::vector<std::vector<pers_pair_t> >& ph) const {

  Ripser(points, Threshold, MaxDim, InputIsDistanceMatrix, ph);

  return 0;
}