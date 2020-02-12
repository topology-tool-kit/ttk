#include <ExplicitTriangulation.h>

using namespace std;
using namespace ttk;

ExplicitTriangulation::ExplicitTriangulation() {

  setDebugMsgPrefix("ExplicitTriangulation");

  clear();
}

ExplicitTriangulation::~ExplicitTriangulation() {
}

int ExplicitTriangulation::clear() {

  AbstractTriangulation::clear();

  vertexNumber_ = 0;
  cellNumber_ = 0;
  doublePrecision_ = false;

  printMsg(
    "[ExplicitTriangulation] Triangulation cleared.", debug::Priority::DETAIL);

  return AbstractTriangulation::clear();
}
