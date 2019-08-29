/// \defgroup vtk vtk
/// \brief The Topology ToolKit - VTK wrapping code for the processing
/// packages.
/// @{
/// \ingroup vtk
/// \class ttkWrapper
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2016.
///

#pragma once

#include <vtkDataSetAlgorithm.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#include <Wrapper.h>
#include <ttkTriangulation.h>

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

// Macros for vtkWrappers
#define TTK_POLY_DATA_NEW(i, ouputInformation, dataTYpe)           \
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

#define TTK_UNSTRUCTURED_GRID_NEW(i, ouputInformation, dataTYpe)   \
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

#define TTK_OUTPUT_MANAGEMENT()                                        \
protected:                                                             \
  int RequestDataObject(vtkInformation *request,                       \
                        vtkInformationVector **inputVector,            \
                        vtkInformationVector *outputVector) override { \
                                                                       \
    vtkDataSetAlgorithm::RequestDataObject(                            \
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
                                                                       \
        TTK_POLY_DATA_NEW(i, outputInformation, dataType);             \
      }                                                                \
    }                                                                  \
                                                                       \
    return 1;                                                          \
  }                                                                    \
                                                                       \
public:

#define TTK_PIPELINE_REQUEST()                                                 \
protected:                                                                     \
  std::vector<vtkSmartPointer<ttkTriangulationFilter>> inputTriangulations_;   \
  int RequestData(vtkInformation *request, vtkInformationVector **inputVector, \
                  vtkInformationVector *outputVector) override {               \
                                                                               \
    if((int)inputTriangulations_.size() != GetNumberOfInputPorts()) {          \
      inputTriangulations_.resize(GetNumberOfInputPorts());                    \
      for(int i = 0; i < (int)inputTriangulations_.size(); i++) {              \
        inputTriangulations_[i]                                                \
          = vtkSmartPointer<ttkTriangulationFilter>::New();                    \
      }                                                                        \
    }                                                                          \
                                                                               \
    std::vector<vtkDataSet *> inputs(GetNumberOfInputPorts(), NULL);           \
    std::vector<vtkDataSet *> outputs(GetNumberOfOutputPorts(), NULL);         \
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
    if(((int)inputs.size() == GetNumberOfInputPorts())                         \
       && ((int)outputs.size() == GetNumberOfOutputPorts()))                   \
      doIt(inputs, outputs);                                                   \
                                                                               \
    return 1;                                                                  \
  }

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

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTriangulationFilter
#else
class ttkTriangulationFilter
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkTriangulationFilter *New();

  vtkTypeMacro(ttkTriangulationFilter, vtkDataSetAlgorithm);

protected:
  ttkTriangulationFilter();
  ~ttkTriangulationFilter(){};

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int RequestDataObject(vtkInformation *request,
                        vtkInformationVector **inputVector,
                        vtkInformationVector *outputVector) override;

private:
  bool needsToAbort() override;

  int updateProgress(const float &progress) override;
};

/// @}
