#include <DimensionReduction.h>
#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)

#ifdef TTK_ENABLE_SCIKIT_LEARN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif

using namespace std;
using namespace ttk;

DimensionReduction::DimensionReduction()
  : numberOfRows_{0}, numberOfColumns_{0}, numberOfComponents_{0},
    numberOfNeighbors_{0}, randomState_{0}, matrix_{nullptr},
    embedding_{nullptr}, majorVersion_{'0'} {
#ifdef TTK_ENABLE_SCIKIT_LEARN
  auto finalize_callback = []() { Py_Finalize(); };

  if(!Py_IsInitialized()) {
    Py_Initialize();
    atexit(finalize_callback);
  }

  const char *version = Py_GetVersion();
  if(version[0] >= '3') {
    stringstream msg;
    msg << "[DimensionReduction] Initializing Python: " << version[0]
        << version[1] << version[2] << endl;
    dMsg(cout, msg.str(), infoMsg);
  } else {
    cerr << "[DimensionReduction] Error: Python 3+ is required:\n"
         << version << " is provided." << endl;
  }

  majorVersion_ = version[0];
#endif
}

DimensionReduction::~DimensionReduction() {
}

bool DimensionReduction::isPythonFound() const {
#ifdef TTK_ENABLE_SCIKIT_LEARN
  return true;
#else
  stringstream msg;
  msg << "[DimensionReduction] "
      << "Warning: scikit-learn support disabled :(" << endl;
  msg << "[DimensionReduction] "
      << "Python/Numpy may not be installed properly." << endl;
  msg << "[DimensionReduction] Features disabled..." << endl;
  dMsg(cerr, msg.str(), fatalMsg);
  return false;
#endif
}

int DimensionReduction::execute() const {
#ifdef TTK_ENABLE_SCIKIT_LEARN
#ifndef TTK_ENABLE_KAMIKAZE
  if(majorVersion_ < '3')
    return -1;
  if(modulePath_.length() <= 0)
    return -1;
  if(moduleName_.length() <= 0)
    return -1;
  if(functionName_.length() <= 0)
    return -1;
  if(!matrix_)
    return -1;
#endif

  Timer t;

  const int numberOfComponents = std::max(2, numberOfComponents_);
  const int numberOfNeighbors = std::max(1, numberOfNeighbors_);

  // declared here to avoid crossing initialization with goto
  vector<PyObject *> gc;
  PyObject *pArray;
  PyObject *pPath;
  PyObject *pSys;
  PyObject *pName;
  PyObject *pModule;
  PyObject *pFunc;
  PyObject *pMethod;
  PyObject *pNumberOfComponents;
  PyObject *pNumberOfNeighbors;
  PyObject *pJobs;
  PyObject *pIsDeterministic;
  PyObject *pReturn;
  PyObject *pNRows;
  PyObject *pNColumns;
  PyObject *pEmbedding;
  PyObject *pSEParams;
  PyObject *pLLEParams;
  PyObject *pMDSParams;
  PyObject *pTSNEParams;
  PyObject *pISOParams;
  PyObject *pPCAParams;
  PyObject *pParams;
  PyArrayObject *npArr;
  PyArrayObject *npEmbedding;

  string modulePath;

  if(PyArray_API == NULL) {
    import_array1(-1);
  }

  // convert the input matrix into a NumPy array.
  const int numberOfDimensions = 2;
  npy_intp dimensions[2]{numberOfRows_, numberOfColumns_};

  pArray = PyArray_SimpleNewFromData(
    numberOfDimensions, dimensions, NPY_DOUBLE, matrix_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pArray) {
    cerr << "[DimensionReduction] Python error: failed to convert the array."
         << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pArray);

  npArr = reinterpret_cast<PyArrayObject *>(pArray);

  pSys = PyImport_ImportModule("sys");
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pSys) {
    cerr << "[DimensionReduction] Python error: failed to load the sys module."
         << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pSys);

  pPath = PyObject_GetAttrString(pSys, "path");
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pPath) {
    cerr
      << "[DimensionReduction] Python error: failed to get the path variable."
      << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pPath);

  if(modulePath_ == "default")
    modulePath = VALUE(TTK_SCRIPTS_PATH);
  else
    modulePath = modulePath_;

  {
    stringstream msg;
    msg << "[DimensionReduction] Loading Python script from: " << modulePath
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
  PyList_Append(pPath, PyUnicode_FromString(modulePath.data()));

  // set other parameters
  pNumberOfComponents = PyLong_FromLong(numberOfComponents);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pNumberOfComponents) {
    cerr << "[DimensionReduction] Python error: cannot convert "
            "pNumberOfComponents."
         << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pNumberOfComponents);

  pNumberOfNeighbors = PyLong_FromLong(numberOfNeighbors);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pNumberOfNeighbors) {
    cerr
      << "[DimensionReduction] Python error: cannot convert pNumberOfNeighbors."
      << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pNumberOfNeighbors);

  pMethod = PyLong_FromLong(method_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pMethod) {
    cerr << "[DimensionReduction] Python error: cannot convert pMethod."
         << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pMethod);

  pJobs = PyLong_FromLong(threadNumber_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pJobs) {
    cerr << "[DimensionReduction] Python error: cannot convert pJobs." << endl;
    goto collect_garbage;
  }
#endif

  pIsDeterministic = PyLong_FromLong(randomState_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pIsDeterministic) {
    cerr
      << "[DimensionReduction] Python error: cannot convert pIsDeterministic."
      << endl;
    goto collect_garbage;
  }
#endif

  // load module
  pName = PyUnicode_FromString(moduleName_.data());
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pName) {
    cerr << "[DimensionReduction] Python error: moduleName parsing failed."
         << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pName);

  pModule = PyImport_Import(pName);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pModule) {
    cerr << "[DimensionReduction] Python error: module import failed." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pModule);

  // configure function
  pFunc = PyObject_GetAttrString(pModule, functionName_.data());
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pFunc) {
    cerr << "[DimensionReduction] Python error: functionName parsing failed."
         << endl;
    goto collect_garbage;
  }

  if(!PyCallable_Check(pFunc)) {
    cerr << "[DimensionReduction] Python error: function call failed." << endl;
    goto collect_garbage;
  }
#endif

  pSEParams = PyList_New(0);
  PyList_Append(pSEParams, PyUnicode_FromString(se_Affinity.data()));
  PyList_Append(pSEParams, PyFloat_FromDouble(se_Gamma));
  PyList_Append(pSEParams, PyUnicode_FromString(se_EigenSolver.data()));

  pLLEParams = PyList_New(0);
  PyList_Append(pLLEParams, PyFloat_FromDouble(lle_Regularization));
  PyList_Append(pLLEParams, PyUnicode_FromString(lle_EigenSolver.data()));
  PyList_Append(pLLEParams, PyFloat_FromDouble(lle_Tolerance));
  PyList_Append(pLLEParams, PyLong_FromLong(lle_MaxIteration));
  PyList_Append(pLLEParams, PyUnicode_FromString(lle_Method.data()));
  PyList_Append(pLLEParams, PyFloat_FromDouble(lle_HessianTolerance));
  PyList_Append(pLLEParams, PyFloat_FromDouble(lle_ModifiedTolerance));
  PyList_Append(
    pLLEParams, PyUnicode_FromString(lle_NeighborsAlgorithm.data()));

  pMDSParams = PyList_New(0);
  PyList_Append(pMDSParams, PyBool_FromLong(mds_Metric));
  PyList_Append(pMDSParams, PyLong_FromLong(mds_Init));
  PyList_Append(pMDSParams, PyLong_FromLong(mds_MaxIteration));
  PyList_Append(pMDSParams, PyLong_FromLong(mds_Verbose));
  PyList_Append(pMDSParams, PyFloat_FromDouble(mds_Epsilon));
  PyList_Append(pMDSParams, PyUnicode_FromString(mds_Dissimilarity.data()));

  pTSNEParams = PyList_New(0);
  PyList_Append(pTSNEParams, PyFloat_FromDouble(tsne_Perplexity));
  PyList_Append(pTSNEParams, PyFloat_FromDouble(tsne_Exaggeration));
  PyList_Append(pTSNEParams, PyFloat_FromDouble(tsne_LearningRate));
  PyList_Append(pTSNEParams, PyLong_FromLong(tsne_MaxIteration));
  PyList_Append(pTSNEParams, PyLong_FromLong(tsne_MaxIterationProgress));
  PyList_Append(pTSNEParams, PyFloat_FromDouble(tsne_GradientThreshold));
  PyList_Append(pTSNEParams, PyUnicode_FromString(tsne_Metric.data()));
  PyList_Append(pTSNEParams, PyUnicode_FromString(tsne_Init.data()));
  PyList_Append(pTSNEParams, PyLong_FromLong(tsne_Verbose));
  PyList_Append(pTSNEParams, PyUnicode_FromString(tsne_Method.data()));
  PyList_Append(pTSNEParams, PyFloat_FromDouble(tsne_Angle));

  pISOParams = PyList_New(0);
  PyList_Append(pISOParams, PyUnicode_FromString(iso_EigenSolver.data()));
  PyList_Append(pISOParams, PyFloat_FromDouble(iso_Tolerance));
  PyList_Append(pISOParams, PyLong_FromLong(iso_MaxIteration));
  PyList_Append(pISOParams, PyUnicode_FromString(iso_PathMethod.data()));
  PyList_Append(
    pISOParams, PyUnicode_FromString(iso_NeighborsAlgorithm.data()));

  pPCAParams = PyList_New(0);
  PyList_Append(pPCAParams, PyBool_FromLong(pca_Copy));
  PyList_Append(pPCAParams, PyBool_FromLong(pca_Whiten));
  PyList_Append(pPCAParams, PyUnicode_FromString(pca_SVDSolver.data()));
  PyList_Append(pPCAParams, PyFloat_FromDouble(pca_Tolerance));
  PyList_Append(pPCAParams, PyUnicode_FromString(pca_MaxIteration.data()));

  pParams = PyList_New(0);
  gc.push_back(pParams);

  PyList_Append(pParams, pSEParams);
  PyList_Append(pParams, pLLEParams);
  PyList_Append(pParams, pMDSParams);
  PyList_Append(pParams, pTSNEParams);
  PyList_Append(pParams, pISOParams);
  PyList_Append(pParams, pPCAParams);

  pReturn = PyObject_CallFunctionObjArgs(
    pFunc, npArr, pMethod, pNumberOfComponents, pNumberOfNeighbors, pJobs,
    pIsDeterministic, pParams, NULL);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pReturn) {
    cerr
      << "[DimensionReduction] Python error: function returned invalid object."
      << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pReturn);

  pNRows = PyList_GetItem(pReturn, 0);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pNRows) {
    cerr
      << "[DimensionReduction] Python error: function returned invalid number "
      << "of rows." << endl;
    goto collect_garbage;
  }
#endif

  pNColumns = PyList_GetItem(pReturn, 1);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pNColumns) {
    cerr << "[DimensionReduction] Python error: function returned invalid "
         << "number of columns." << endl;
    goto collect_garbage;
  }
#endif

  pEmbedding = PyList_GetItem(pReturn, 2);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pEmbedding) {
    cerr << "[DimensionReduction] Python error: function returned invalid"
         << " embedding data." << endl;
    goto collect_garbage;
  }
#endif

  if(PyLong_AsLong(pNRows) == numberOfRows_
     and PyLong_AsLong(pNColumns) == numberOfComponents) {
    npEmbedding = reinterpret_cast<PyArrayObject *>(pEmbedding);

    embedding_->resize(numberOfComponents);
    for(int i = 0; i < numberOfComponents; ++i) {
      if(PyArray_TYPE(npEmbedding) == NPY_FLOAT) {
        float *c_out = reinterpret_cast<float *>(PyArray_DATA(npEmbedding));
        for(int j = 0; j < numberOfRows_; ++j)
          (*embedding_)[i].push_back(c_out[i * numberOfRows_ + j]);
      } else if(PyArray_TYPE(npEmbedding) == NPY_DOUBLE) {
        double *c_out = reinterpret_cast<double *>(PyArray_DATA(npEmbedding));
        for(int j = 0; j < numberOfRows_; ++j)
          (*embedding_)[i].push_back(c_out[i * numberOfRows_ + j]);
      }
    }
  }

  // normal control-flow
  for(auto i : gc)
    Py_DECREF(i);

  {
    stringstream msg;
    msg << "[DimensionReduction] ";
    switch(method_) {
      case 0:
        msg << "SE";
        break;
      case 1:
        msg << "LLE";
        break;
      case 2:
        msg << "MDS";
        break;
      case 3:
        msg << "t-SNE";
        break;
      case 4:
        msg << "IsoMap";
        break;
      case 5:
        msg << "PCA";
        break;
    }
    msg << " computed in " << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;

  // error control-flow
#ifndef TTK_ENABLE_KAMIKAZE
collect_garbage:
#endif
  for(auto i : gc)
    Py_DECREF(i);
  return -1;

#endif

  return 0;
}
