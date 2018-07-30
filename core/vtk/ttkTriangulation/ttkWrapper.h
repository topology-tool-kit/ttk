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

#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkImageData.h>
#include                  <vtkPolyData.h>
#include                  <vtkUnstructuredGrid.h>

#include                  <Wrapper.h>
#include                  <ttkTriangulation.h>

#ifdef TTK_USE_64BIT_IDS
using ttkSimplexIdTypeArray = vtkIdTypeArray;
#else
using ttkSimplexIdTypeArray = vtkIntArray;
#endif

#define TTK_COMMA ,

// Macros for vtkWrappers
#define TTK_POLY_DATA_NEW(i, ouputInformation, dataTYpe)\
  if(dataType == "vtkPolyData"){\
    ttkPolyData *data = ttkPolyData::SafeDownCast(\
      outputInformation->Get(vtkDataObject::DATA_OBJECT()));\
    \
    if(!data){\
      data = ttkPolyData::New();\
      outputInformation->Set(vtkDataObject::DATA_OBJECT(), data);\
      data->FastDelete();\
      data->CopyInformationFromPipeline(outputInformation);\
      GetOutputPortInformation(i)->Set(\
        vtkDataObject::DATA_EXTENT_TYPE(), data->GetExtentType());\
    }\
    \
    if((data)&&(!data->getTriangulation())){\
      data->allocate();\
    }\
  }

#define TTK_UNSTRUCTURED_GRID_NEW(i, ouputInformation, dataTYpe)\
  if(dataType == "vtkUnstructuredGrid"){\
    \
    ttkUnstructuredGrid *data = ttkUnstructuredGrid::SafeDownCast(\
      outputInformation->Get(vtkDataObject::DATA_OBJECT()));\
    \
    if(!data){\
      data = ttkUnstructuredGrid::New();\
      outputInformation->Set(vtkDataObject::DATA_OBJECT(), data);\
      data->FastDelete();\
      data->CopyInformationFromPipeline(outputInformation);\
      GetOutputPortInformation(i)->Set(\
        vtkDataObject::DATA_EXTENT_TYPE(), data->GetExtentType());\
    }\
    \
    if((data)&&(!data->getTriangulation())){\
      data->allocate();\
    }\
  }

#define TTK_OUTPUT_MANAGEMENT() \
  protected: \
    int RequestDataObject(vtkInformation *request, \
      vtkInformationVector **inputVector, vtkInformationVector *outputVector) override {\
      \
      vtkDataSetAlgorithm::RequestDataObject(\
        request, inputVector, outputVector);\
        \
      for(int i = 0; i < GetNumberOfOutputPorts(); i++){\
        \
        vtkInformation *outputInformation =\
          outputVector->GetInformationObject(i);\
        \
        vtkDataSet *output = vtkDataSet::SafeDownCast(\
          outputInformation->Get(vtkDataObject::DATA_OBJECT()));\
        \
        if(output){\
          \
          std::string dataType = output->GetClassName();\
          \
          TTK_UNSTRUCTURED_GRID_NEW(i, outputInformation, dataType);\
          \
          TTK_POLY_DATA_NEW(i, outputInformation, dataType);\
        }\
      }\
      \
      return 1;\
    }\
    public:
      
#define TTK_PIPELINE_REQUEST() \
    protected:\
      std::vector<vtkSmartPointer<ttkTriangulationFilter> > inputTriangulations_;\
      int RequestData(vtkInformation *request, \
        vtkInformationVector **inputVector, \
        vtkInformationVector *outputVector) override {\
        \
        if((int) inputTriangulations_.size() != GetNumberOfInputPorts()){\
          inputTriangulations_.resize(GetNumberOfInputPorts());\
          for(int i = 0; i < (int) inputTriangulations_.size(); i++){\
            inputTriangulations_[i] = \
              vtkSmartPointer<ttkTriangulationFilter>::New();\
          }\
        }\
        \
        std::vector<vtkDataSet *> inputs(GetNumberOfInputPorts(), NULL);\
        std::vector<vtkDataSet *> outputs(GetNumberOfOutputPorts(), NULL);\
        \
        for(int i = 0; i < GetNumberOfInputPorts(); i++){\
          vtkDataSet *input = vtkDataSet::GetData(inputVector[i]);\
          if(input){\
            inputTriangulations_[i]->SetInputData(input);\
            inputTriangulations_[i]->Update();\
            inputs[i] = inputTriangulations_[i]->GetOutput();\
          }\
        }\
        for(int i = 0; i < GetNumberOfOutputPorts(); i++){\
          outputs[i] = vtkDataSet::SafeDownCast(\
          \
          outputVector->GetInformationObject(i)->Get(\
            vtkDataObject::DATA_OBJECT()));\
        }\
        \
        if(((int) inputs.size() == GetNumberOfInputPorts())\
          &&((int) outputs.size() == GetNumberOfOutputPorts()))\
          doIt(inputs, outputs);\
        \
        return 1;\
      }

#define TTK_SETUP()\
    public:\
      TTK_PIPELINE_REQUEST();\
      TTK_OUTPUT_MANAGEMENT();\
      void SetThreads(){\
        if(!UseAllCores)\
          threadNumber_ = ThreadNumber;\
        else{\
          threadNumber_ = ttk::OsCall::getNumberOfCores();\
        }\
        Modified();\
      }\
    private:\
      int doIt(std::vector<vtkDataSet *> &inputs, std::vector<vtkDataSet *> &outputs);\
      \
      bool needsToAbort() override { return GetAbortExecute();};\
      \
      int updateProgress(const float &progress) override { \
        UpdateProgress(progress);\
        return 0;\
      };\
    protected:\
      bool                  UseAllCores;\
      int                   ThreadNumber;\
      
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTriangulationFilter
#else
class ttkTriangulationFilter
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
      
    static ttkTriangulationFilter* New();
    
    vtkTypeMacro(ttkTriangulationFilter, vtkDataSetAlgorithm);
    
  protected:
   
    ttkTriangulationFilter();
    ~ttkTriangulationFilter(){};
    
    int RequestData(vtkInformation *request, 
      vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;
   
    int RequestDataObject(vtkInformation *request,
      vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;
    
  private:
    
    bool needsToAbort() override;
    
    int updateProgress(const float &progress) override;
   
};

/// @}
