#include <DistanceField.h>

using namespace std;
using namespace ttk;

DistanceField::DistanceField()
  : vertexNumber_{}, sourceNumber_{}, vertexIdentifierScalarFieldPointer_{},
    outputScalarFieldPointer_{}, outputIdentifiers_{}, outputSegmentation_{} {
  this->setDebugMsgPrefix("DistanceField");
}

DistanceField::~DistanceField() {
}
