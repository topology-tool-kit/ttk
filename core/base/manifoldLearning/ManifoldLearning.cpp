#include<ManifoldLearning.h>
#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include<Python.h>
#include<numpy/arrayobject.h>

using namespace std;
using namespace ttk;

ManifoldLearning::ManifoldLearning():
  matrixDimension_{0},
  numberOfComponents_{0},
  numberOfNeighbors_{0},
  matrix_{nullptr},
  components_{nullptr},
  majorVersion_{'0'}
{
  auto finalize_callback=[](){
    Py_Finalize();
  };

  if(!Py_IsInitialized()){
    Py_Initialize();
    atexit(finalize_callback);
  }

  const char* version=Py_GetVersion();
  if(version[0]>='3'){
    stringstream msg;
    msg << "[ManifoldLearning] Initializing Python: " <<
      version[0] << version[1] << version[2] << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
  else{
    cerr << "[ManifoldLearning] Error: Python 3+ is required:\n" <<
      version << " is provided." << endl;
  }

  majorVersion_=version[0];
}

ManifoldLearning::~ManifoldLearning(){
}

int ManifoldLearning::execute() const{
  return 0;
}
