#include <IntegralLines.h>

using namespace std;
using namespace ttk;

IntegralLines::IntegralLines()
  : vertexNumber_{}, seedNumber_{}, inputScalarField_{}, inputOffsets_{},
    vertexIdentifierScalarField_{}, outputTrajectories_{} {
  this->setDebugMsgPrefix("IntegralLines");
}

IntegralLines::~IntegralLines() = default;
