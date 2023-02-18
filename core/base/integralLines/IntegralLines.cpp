#include <IntegralLines.h>

using namespace std;
using namespace ttk;

IntegralLines::IntegralLines()
  : vertexNumber_{}, seedNumber_{}, inputScalarField_{}, inputOffsets_{},
    vertexIdentifierScalarField_{}, outputIntegralLines_{} {
  this->setDebugMsgPrefix("IntegralLines");
#ifdef TTK_ENABLE_MPI
  hasMPISupport_ = true;
#endif
}

IntegralLines::~IntegralLines() = default;
