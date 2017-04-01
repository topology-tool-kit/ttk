#include                  <Triangulation.h>

Triangulation::Triangulation(){

  debugLevel_ = 0;
  gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;
  abstractTriangulation_ = NULL;
}

Triangulation::~Triangulation(){
  
}

