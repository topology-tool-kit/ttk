#include <ContinuousScatterPlot.h>

using namespace std;
using namespace ttk;

ContinuousScatterPlot::ContinuousScatterPlot()
  : vertexNumber_{}, withDummyValue_{}, dummyValue_{}, resolutions_{},
    scalarMin_{}, scalarMax_{}, density_{} {
  this->setDebugMsgPrefix("ContinuousScatterPlot");
}

ContinuousScatterPlot::~ContinuousScatterPlot() = default;