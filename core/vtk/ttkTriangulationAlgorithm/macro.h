#pragma once

#include <vtkVersion.h>
#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#define TTK_COMMA ,

#define ttkVtkTemplateMacroCase(                         \
  dataType, triangulationType, triangulationClass, call) \
  case triangulationType: {                              \
    typedef triangulationClass TTK_TT;                   \
    switch(dataType) { vtkTemplateMacro((call)); };      \
  }; break;

#define ttkVtkTemplateMacro(dataType, triangulationType, call)            \
  switch(triangulationType) {                                             \
    ttkVtkTemplateMacroCase(dataType, ttk::Triangulation::Type::EXPLICIT, \
                            ttk::ExplicitTriangulation, call);            \
    ttkVtkTemplateMacroCase(dataType, ttk::Triangulation::Type::IMPLICIT, \
                            ttk::ImplicitTriangulation, call);            \
    ttkVtkTemplateMacroCase(dataType, ttk::Triangulation::Type::PERIODIC, \
                            ttk::PeriodicImplicitTriangulation, call);    \
  }

#if VTK_MAJOR_VERSION <= 8
#define vtkTemplate2Macro(call)                                             \
  vtkTemplate2MacroCase1(VTK_DOUBLE, double, call);                         \
  vtkTemplate2MacroCase1(VTK_FLOAT, float, call);                           \
  vtkTemplate2MacroCase1(VTK_LONG_LONG, long long, call);                   \
  vtkTemplate2MacroCase1(VTK_UNSIGNED_LONG_LONG, unsigned long long, call); \
  vtkTemplate2MacroCase1(VTK_ID_TYPE, vtkIdType, call);                     \
  vtkTemplate2MacroCase1(VTK_LONG, long, call);                             \
  vtkTemplate2MacroCase1(VTK_UNSIGNED_LONG, unsigned long, call);           \
  vtkTemplate2MacroCase1(VTK_INT, int, call);                               \
  vtkTemplate2MacroCase1(VTK_UNSIGNED_INT, unsigned int, call);             \
  vtkTemplate2MacroCase1(VTK_SHORT, short, call);                           \
  vtkTemplate2MacroCase1(VTK_UNSIGNED_SHORT, unsigned short, call);         \
  vtkTemplate2MacroCase1(VTK_CHAR, char, call);                             \
  vtkTemplate2MacroCase1(VTK_SIGNED_CHAR, signed char, call);               \
  vtkTemplate2MacroCase1(VTK_UNSIGNED_CHAR, unsigned char, call)
#ifndef vtkTemplate2MacroCase1
#define vtkTemplate2MacroCase1(type1N, type1, call)                            \
  vtkTemplate2MacroCase2(type1N, type1, VTK_DOUBLE, double, call);             \
  vtkTemplate2MacroCase2(type1N, type1, VTK_FLOAT, float, call);               \
  vtkTemplate2MacroCase2(type1N, type1, VTK_LONG_LONG, long long, call);       \
  vtkTemplate2MacroCase2(                                                      \
    type1N, type1, VTK_UNSIGNED_LONG_LONG, unsigned long long, call);          \
  vtkTemplate2MacroCase2(type1N, type1, VTK_ID_TYPE, vtkIdType, call);         \
  vtkTemplate2MacroCase2(type1N, type1, VTK_LONG, long, call);                 \
  vtkTemplate2MacroCase2(                                                      \
    type1N, type1, VTK_UNSIGNED_LONG, unsigned long, call);                    \
  vtkTemplate2MacroCase2(type1N, type1, VTK_INT, int, call);                   \
  vtkTemplate2MacroCase2(type1N, type1, VTK_UNSIGNED_INT, unsigned int, call); \
  vtkTemplate2MacroCase2(type1N, type1, VTK_SHORT, short, call);               \
  vtkTemplate2MacroCase2(                                                      \
    type1N, type1, VTK_UNSIGNED_SHORT, unsigned short, call);                  \
  vtkTemplate2MacroCase2(type1N, type1, VTK_CHAR, char, call);                 \
  vtkTemplate2MacroCase2(type1N, type1, VTK_SIGNED_CHAR, signed char, call);   \
  vtkTemplate2MacroCase2(type1N, type1, VTK_UNSIGNED_CHAR, unsigned char, call)
#endif
#define vtkTemplate2MacroCase2(type1N, type1, type2N, type2, call) \
  case vtkTemplate2PackMacro(type1N, type2N): {                    \
    typedef type1 VTK_T1;                                          \
    typedef type2 VTK_T2;                                          \
    call;                                                          \
  }; break
#ifndef vtkTemplate2PackMacro
#define vtkTemplate2PackMacro(type1N, type2N) \
  ((((type1N)&0xFF) << 8) | ((type2N)&0xFF))
#endif
#endif

#ifdef TTK_ENABLE_64BIT_IDS
using ttkSimplexIdTypeArray = vtkIdTypeArray;
#else
using ttkSimplexIdTypeArray = vtkIntArray;
#endif

//------------------------------------------------------------------

// Macros for vtkWrappers
#define TTK_IMAGE_DATA_NEW(i, ouputInformation, dataType)          \
  if(dataType == "vtkImageData") {                                 \
    ttkImageData *data = ttkImageData::SafeDownCast(               \
      outputInformation->Get(vtkDataObject::DATA_OBJECT()));       \
                                                                   \
    if(!data) {                                                    \
      data = ttkImageData::New();                                  \
      outputInformation->Set(vtkDataObject::DATA_OBJECT(), data);  \
      data->FastDelete();                                          \
      data->CopyInformationFromPipeline(outputInformation);        \
      GetOutputPortInformation(i)->Set(                            \
        vtkDataObject::DATA_EXTENT_TYPE(), data->GetExtentType()); \
    }                                                              \
                                                                   \
    if((data) && (!data->getTriangulation())) {                    \
      data->allocate();                                            \
    }                                                              \
  }

//------------------------------------------------------------------

#define TTK_POLY_DATA_NEW(i, ouputInformation, dataType)           \
  if(dataType == "vtkPolyData") {                                  \
    ttkPolyData *data = ttkPolyData::SafeDownCast(                 \
      outputInformation->Get(vtkDataObject::DATA_OBJECT()));       \
                                                                   \
    if(!data) {                                                    \
      data = ttkPolyData::New();                                   \
      outputInformation->Set(vtkDataObject::DATA_OBJECT(), data);  \
      data->FastDelete();                                          \
      data->CopyInformationFromPipeline(outputInformation);        \
      GetOutputPortInformation(i)->Set(                            \
        vtkDataObject::DATA_EXTENT_TYPE(), data->GetExtentType()); \
    }                                                              \
                                                                   \
    if((data) && (!data->getTriangulation())) {                    \
      data->allocate();                                            \
    }                                                              \
  }

//------------------------------------------------------------------

#define TTK_UNSTRUCTURED_GRID_NEW(i, ouputInformation, dataType)   \
  if(dataType == "vtkUnstructuredGrid") {                          \
                                                                   \
    ttkUnstructuredGrid *data = ttkUnstructuredGrid::SafeDownCast( \
      outputInformation->Get(vtkDataObject::DATA_OBJECT()));       \
                                                                   \
    if(!data) {                                                    \
      data = ttkUnstructuredGrid::New();                           \
      outputInformation->Set(vtkDataObject::DATA_OBJECT(), data);  \
      data->FastDelete();                                          \
      data->CopyInformationFromPipeline(outputInformation);        \
      GetOutputPortInformation(i)->Set(                            \
        vtkDataObject::DATA_EXTENT_TYPE(), data->GetExtentType()); \
    }                                                              \
                                                                   \
    if((data) && (!data->getTriangulation())) {                    \
      data->allocate();                                            \
    }                                                              \
  }

//----------------------------------------------------------------------

#define TTK_OUTPUT_MANAGEMENT()                                        \
protected:                                                             \
  int RequestDataObject(vtkInformation *request,                       \
                        vtkInformationVector **inputVector,            \
                        vtkInformationVector *outputVector) override { \
                                                                       \
    int ret = vtkDataSetAlgorithm::RequestDataObject(                  \
      request, inputVector, outputVector);                             \
                                                                       \
    for(int i = 0; i < GetNumberOfOutputPorts(); i++) {                \
                                                                       \
      vtkInformation *outputInformation                                \
        = outputVector->GetInformationObject(i);                       \
                                                                       \
      vtkDataSet *output = vtkDataSet::SafeDownCast(                   \
        outputInformation->Get(vtkDataObject::DATA_OBJECT()));         \
                                                                       \
      if(output) {                                                     \
                                                                       \
        std::string dataType = output->GetClassName();                 \
                                                                       \
        TTK_UNSTRUCTURED_GRID_NEW(i, outputInformation, dataType);     \
        TTK_POLY_DATA_NEW(i, outputInformation, dataType);             \
      }                                                                \
    }                                                                  \
                                                                       \
    return ret;                                                        \
  }                                                                    \
                                                                       \
public:

//-------------------------------------------------------------------------------

#define TTK_PIPELINE_REQUEST()                                                 \
protected:                                                                     \
  std::vector<vtkSmartPointer<ttkTriangulationAlgorithm>>                      \
    inputTriangulations_;                                                      \
  int RequestData(vtkInformation *request, vtkInformationVector **inputVector, \
                  vtkInformationVector *outputVector) override {               \
                                                                               \
    if((int)inputTriangulations_.size() != GetNumberOfInputPorts()) {          \
      inputTriangulations_.resize(GetNumberOfInputPorts());                    \
      for(int i = 0; i < (int)inputTriangulations_.size(); i++) {              \
        inputTriangulations_[i]                                                \
          = vtkSmartPointer<ttkTriangulationAlgorithm>::New();                 \
      }                                                                        \
    }                                                                          \
                                                                               \
    std::vector<vtkDataSet *> inputs(GetNumberOfInputPorts(), nullptr);        \
    std::vector<vtkDataSet *> outputs(GetNumberOfOutputPorts(), nullptr);      \
                                                                               \
    for(int i = 0; i < GetNumberOfInputPorts(); i++) {                         \
      vtkDataSet *input = vtkDataSet::GetData(inputVector[i]);                 \
      if(input) {                                                              \
        inputTriangulations_[i]->SetInputData(input);                          \
        inputTriangulations_[i]->Update();                                     \
        inputs[i] = inputTriangulations_[i]->GetOutput();                      \
      }                                                                        \
    }                                                                          \
    for(int i = 0; i < GetNumberOfOutputPorts(); i++) {                        \
      outputs[i] = vtkDataSet::SafeDownCast(                                   \
                                                                               \
        outputVector->GetInformationObject(i)->Get(                            \
          vtkDataObject::DATA_OBJECT()));                                      \
    }                                                                          \
                                                                               \
    int ret = 1;                                                               \
    if(((int)inputs.size() == GetNumberOfInputPorts())                         \
       && ((int)outputs.size() == GetNumberOfOutputPorts()))                   \
      ret = doIt(inputs, outputs);                                             \
                                                                               \
    return !ret;                                                               \
  }

//-------------------------------------------------------------------------------

#define TTK_SETUP()                                                         \
public:                                                                     \
  TTK_PIPELINE_REQUEST();                                                   \
  TTK_OUTPUT_MANAGEMENT();                                                  \
  void SetThreads() {                                                       \
    if(!UseAllCores)                                                        \
      threadNumber_ = ThreadNumber;                                         \
    else {                                                                  \
      threadNumber_ = ttk::OsCall::getNumberOfCores();                      \
    }                                                                       \
    Modified();                                                             \
  }                                                                         \
                                                                            \
private:                                                                    \
  int doIt(                                                                 \
    std::vector<vtkDataSet *> &inputs, std::vector<vtkDataSet *> &outputs); \
                                                                            \
  bool needsToAbort() override {                                            \
    return GetAbortExecute();                                               \
  };                                                                        \
                                                                            \
  int updateProgress(const float &progress) override {                      \
    UpdateProgress(progress);                                               \
    return 0;                                                               \
  };                                                                        \
                                                                            \
protected:                                                                  \
  bool UseAllCores;                                                         \
  int ThreadNumber;

//----------------------------------------------------------------------------

// Internal things to allow the ttkTriangulation to travel through a VTK
// pipeline.
#define ttkTypeMacro(thisClass, superClass)                      \
protected:                                                       \
  const char *GetClassNameInternal() const override {            \
    return #superClass;                                          \
  }                                                              \
                                                                 \
public:                                                          \
  typedef superClass Superclass;                                 \
  static bool IsTypeOf(const char *type) {                       \
    if(!strcmp("superClass", type)) {                            \
      return 1;                                                  \
    }                                                            \
    return superClass::IsTypeOf(type);                           \
  }                                                              \
  int IsA(const char *type) override {                           \
    return this->thisClass::IsTypeOf(type);                      \
  }                                                              \
  static thisClass *SafeDownCast(vtkObjectBase *o) {             \
    if((o) && (o->IsA("thisClass"))) {                           \
      return static_cast<thisClass *>(o);                        \
    }                                                            \
    return NULL;                                                 \
  }                                                              \
  thisClass *NewInstance() const {                               \
    return thisClass::SafeDownCast(this->NewInstanceInternal()); \
  }                                                              \
                                                                 \
protected:                                                       \
  vtkObjectBase *NewInstanceInternal() const override {          \
    return thisClass::New();                                     \
  }                                                              \
                                                                 \
public:
