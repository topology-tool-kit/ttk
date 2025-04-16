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
  const std::vector<std::vector<double>> &points,
  rpd::MultidimensionalDiagram &ph) const {

  if(isPrime(FieldOfCoefficients))
    ripser::ripser(points, ph, SimplexMaximumDiameter, SimplexMaximumDimension,
                   InputIsDistanceMatrix, false, FieldOfCoefficients);
  else
    printErr("The chosen p=" + std::to_string(FieldOfCoefficients)
             + " is not prime");

  return 0;
}