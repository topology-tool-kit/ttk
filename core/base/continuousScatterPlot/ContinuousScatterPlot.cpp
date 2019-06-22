#include <ContinuousScatterPlot.h>

using namespace std;
using namespace ttk;

ContinuousScatterPlot::ContinuousScatterPlot()
  : vertexNumber_{}, triangulation_{}, withDummyValue_{}, dummyValue_{},
    resolutions_{}, inputScalarField1_{}, inputScalarField2_{}, scalarMin_{},
    scalarMax_{}, density_{} {
}

ContinuousScatterPlot::~ContinuousScatterPlot() {
}
