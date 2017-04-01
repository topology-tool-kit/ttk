#include                  <ExplicitTriangulation.h>

ExplicitTriangulation::ExplicitTriangulation(){

  clear();
}

ExplicitTriangulation::~ExplicitTriangulation(){

}

int ExplicitTriangulation::clear(){
  
  vertexNumber_ = 0;
  cellNumber_ = 0;
  
  {
    stringstream msg;
    msg << "[ExplicitTriangulation] Triangulation cleared." << endl;
    dMsg(cout, msg.str(), detailedInfoMsg);
  }
  
  return AbstractTriangulation::clear();
}
