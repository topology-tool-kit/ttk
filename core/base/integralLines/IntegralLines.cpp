#include <IntegralLines.h>

using namespace std;
using namespace ttk;

IntegralLines::IntegralLines()
  : vertexNumber_{}, seedNumber_{}, triangulation_{}, inputScalarField_{},
    inputOffsets_{}, vertexIdentifierScalarField_{}, outputTrajectories_{} {
}

IntegralLines::~IntegralLines() {
}
