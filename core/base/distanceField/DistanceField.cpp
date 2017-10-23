#include<DistanceField.h>

DistanceField::DistanceField():
  vertexNumber_{},
  sourceNumber_{},
  triangulation_{},
  vertexIdentifierScalarFieldPointer_{},
  outputScalarFieldPointer_{},
  outputIdentifiers_{},
  outputSegmentation_{}
{}

DistanceField::~DistanceField(){
}

