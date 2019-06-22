#include <PersistenceCurve.h>

using namespace std;
using namespace ttk;

PersistenceCurve::PersistenceCurve()
  : ComputeSaddleConnectors{},

    triangulation_{}, inputScalars_{},
    inputOffsets_{}, JTPlot_{}, MSCPlot_{}, STPlot_{}, CTPlot_{} {
}

PersistenceCurve::~PersistenceCurve() {
}
