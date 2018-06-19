#include                  <PersistenceDiagramsBarycenter.h>

using namespace std;
using namespace ttk;

PersistenceDiagramsBarycenter::PersistenceDiagramsBarycenter(){

  inputData_ = NULL;
  vertexNumber_ = 0;
  numberOfInputs_ = 0;
}

PersistenceDiagramsBarycenter::~PersistenceDiagramsBarycenter(){
  if(inputData_)
    free(inputData_);
}

