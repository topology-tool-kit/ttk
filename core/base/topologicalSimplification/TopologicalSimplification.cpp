#include <TopologicalSimplification.h>

using namespace std;
using namespace ttk;

TopologicalSimplification::TopologicalSimplification()
  : triangulation_{}, vertexNumber_{}, constraintNumber_{},
    inputScalarFieldPointer_{}, vertexIdentifierScalarFieldPointer_{},
    inputOffsetScalarFieldPointer_{}, considerIdentifierAsBlackList_{},
    addPerturbation_{}, outputScalarFieldPointer_{},
    outputOffsetScalarFieldPointer_{} {
  considerIdentifierAsBlackList_ = false;
  addPerturbation_ = false;
}

TopologicalSimplification::~TopologicalSimplification() {
}
