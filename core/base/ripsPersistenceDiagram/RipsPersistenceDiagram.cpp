#include <RipsPersistenceDiagram.h>

static bool isPrime(int n) {
  if(n <= 1)
    return false;
  for(int d = 2; d * d <= n; ++d) {
    if(n % d == 0)
      return false;
  }
  return true;
}

ttk::RipsPersistenceDiagram::RipsPersistenceDiagram() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("RipsPersistenceDiagram");
}

int ttk::RipsPersistenceDiagram::execute(
  const rpd::PointCloud &points,
  rpd::MultidimensionalDiagram &ph,
  std::vector<rpd::Generator> &generators) const {

  if (BackEnd == BACKEND::RIPSER) {
    if(isPrime(FieldOfCoefficients))
      ripser::ripser(points, ph, SimplexMaximumDiameter, SimplexMaximumDimension,
                     InputIsDistanceMatrix, false, FieldOfCoefficients);
    else
      printErr("The chosen p=" + std::to_string(FieldOfCoefficients)
               + " is not prime");
  }

  else if (BackEnd == BACKEND::GEOMETRY) {
#ifdef TTK_ENABLE_CGAL
    if (points[0].size() == 2) {
      rpd::FastRipsPersistenceDiagram2 FRPD(points);
      FRPD.setDebugLevel(debugLevel_);
      if (DelaunayRips)
        FRPD.computeDelaunayRips0And1Persistence(ph);
      else
        FRPD.computeRips0And1Persistence(ph, false, false);
      if (OutputGenerators)
        FRPD.exportRips1Generators(generators);
    }
    else
      printErr("Geometric method only implemented for dim 2.");
#else
    TTK_FORCE_USE(generators);
    this->printErr("TTK was not compiled with CGAL support enabled.");
#endif
  }

  return 0;
}