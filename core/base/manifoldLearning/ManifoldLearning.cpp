#include<ManifoldLearning.h>
#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include<Python.h>
#include<numpy/arrayobject.h>

using namespace std;
using namespace ttk;

ManifoldLearning::ManifoldLearning():
  numberOfRows_{0},
  numberOfColumns_{0},
  numberOfComponents_{0},
  numberOfNeighbors_{0},
  matrix_{nullptr},
  embedding_{nullptr},
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
#ifndef TTK_ENABLE_KAMIKAZE
  if(majorVersion_<'3') return -1;
  if(modulePath_.length()<=0) return -1;
  if(moduleName_.length()<=0) return -1;
  if(functionName_.length()<=0) return -1;
  if(!matrix_) return -1;
  // if(!components_) return -1;
#endif

  const int numberOfComponents=std::max(2,numberOfComponents_);
  const int numberOfNeighbors=std::max(1,numberOfNeighbors_);

  // declared here to avoid crossing initialization with goto
  vector<PyObject*> gc;
  PyObject* pArray;
  PyObject* pPath;
  PyObject* pSys;
  PyObject* pName;
  PyObject* pModule;
  PyObject* pFunc;
  PyObject* pMethod;
  PyObject* pNumberOfComponents;
  PyObject* pNumberOfNeighbors;
  PyObject* pJobs;
  PyObject* pReturn;
  PyObject* pNRows;
  PyObject* pNColumns;
  PyObject* pEmbedding;
  PyArrayObject* npArr;
  PyArrayObject* npEmbedding;

   string modulePath;

  if(PyArray_API==NULL){
    import_array();
  }

  // convert the input matrix into a NumPy array.
  const int numberOfDimensions=2;
  npy_intp dimensions[2]{numberOfRows_, numberOfColumns_};

  pArray=PyArray_SimpleNewFromData(numberOfDimensions,
      dimensions, NPY_DOUBLE, matrix_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pArray){
    cerr << "[ManifoldLearning] Python error: failed to convert the array." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pArray);

  npArr=reinterpret_cast<PyArrayObject*>(pArray);

  pSys=PyImport_ImportModule("sys");
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pSys){
    cerr << "[ManifoldLearning] Python error: failed to load the sys module." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pSys);

  pPath=PyObject_GetAttrString(pSys, "path");
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pPath){
    cerr << "[ManifoldLearning] Python error: failed to get the path variable." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pPath);

  if(modulePath_=="default")
    modulePath=VALUE(TTK_SCRIPTS_PATH);
  else
    modulePath=modulePath_;

  {
    stringstream msg;
    msg << "[ManifoldLearning] Loading Python script from: "
      << modulePath << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
  PyList_Append(pPath, PyUnicode_FromString(modulePath.data()));

  // set other parameters
  pNumberOfComponents=PyLong_FromLong(numberOfComponents);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pNumberOfComponents){
    cerr << "[ManifoldLearning] Python error: cannot convert pNumberOfComponents." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pNumberOfComponents);

  pNumberOfNeighbors=PyLong_FromLong(numberOfNeighbors);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pNumberOfNeighbors){
    cerr << "[ManifoldLearning] Python error: cannot convert pNumberOfNeighbors." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pNumberOfNeighbors);

  pMethod=PyLong_FromLong(method_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pMethod){
    cerr << "[ManifoldLearning] Python error: cannot convert pMethod." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pMethod);

  pJobs=PyLong_FromLong(threadNumber_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pJobs){
    cerr << "[ManifoldLearning] Python error: cannot convert pJobs." << endl;
    goto collect_garbage;
  }
#endif

  // load module
  pName=PyUnicode_FromString(moduleName_.data());
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pName){
    cerr << "[ManifoldLearning] Python error: moduleName parsing failed." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pName);

  pModule=PyImport_Import(pName);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pModule){
    cerr << "[ManifoldLearning] Python error: module import failed." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pModule);

  // configure function
  pFunc=PyObject_GetAttrString(pModule, functionName_.data());
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pFunc){
    cerr << "[ManifoldLearning] Python error: functionName parsing failed." << endl;
    goto collect_garbage;
  }

  if(!PyCallable_Check(pFunc)){
    cerr << "[ManifoldLearning] Python error: function call failed." << endl;
    goto collect_garbage;
  }
#endif

   pReturn=PyObject_CallFunctionObjArgs(pFunc, npArr, pMethod, pNumberOfComponents,
       pNumberOfNeighbors, pJobs, NULL);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pReturn){
    cerr << "[ManifoldLearning] Python error: function returned invalid object." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pReturn);

  pNRows=PyList_GetItem(pReturn, 0);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pNRows){
    cerr << "[SpectralEmbedding] Python error: function returned invalid number of rows." << endl;
    goto collect_garbage;
  }
#endif

  pNColumns=PyList_GetItem(pReturn, 1);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pNColumns){
    cerr << "[SpectralEmbedding] Python error: function returned invalid number of columns." << endl;
    goto collect_garbage;
  }
#endif

  pEmbedding=PyList_GetItem(pReturn, 2);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pEmbedding){
    cerr << "[SpectralEmbedding] Python error: function returned invalid embedding data." << endl;
    goto collect_garbage;
  }
#endif

  if(PyLong_AsLong(pNRows)==numberOfRows_ and PyLong_AsLong(pNColumns)==numberOfComponents){
    npEmbedding=reinterpret_cast<PyArrayObject*>(pEmbedding);

    embedding_->resize(numberOfComponents);
    for(int i=0; i<numberOfComponents; ++i){
      if(PyArray_TYPE(npEmbedding)==NPY_FLOAT){
        float* c_out=reinterpret_cast<float*>(PyArray_DATA(npEmbedding));
        for(int j=0; j<numberOfRows_; ++j)
          (*embedding_)[i].push_back(c_out[i*numberOfRows_+j]);
      }
      else if(PyArray_TYPE(npEmbedding)==NPY_DOUBLE){
        double* c_out=reinterpret_cast<double*>(PyArray_DATA(npEmbedding));
        for(int j=0; j<numberOfRows_; ++j)
          (*embedding_)[i].push_back(c_out[i*numberOfRows_+j]);
      }
    }
  }

#ifndef TTK_ENABLE_KAMIKAZE
collect_garbage:
#endif
  for(auto i : gc)
    Py_DECREF(i);

  return 0;
}
