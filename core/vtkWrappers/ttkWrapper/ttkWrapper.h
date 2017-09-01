/// \defgroup ttkWrappers ttkWrappers
/// \brief The Topology ToolKit - VTK wrapping code for the processing 
/// packages.
/// @{
/// \ingroup ttkWrappers
/// \class ttkWrapper
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2016.
/// 

#pragma once

#include                  <Wrapper.h>
#include                  <ttkTriangulation.h>

#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkImageData.h>
#include                  <vtkPolyData.h>
#include                  <ttkTriangulation.h>
#include                  <vtkUnstructuredGrid.h>

// Macros for ttkWrappers
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
      vtkInformationVector **inputVector, vtkInformationVector *outputVector){\
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
          string dataType = output->GetClassName();\
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
      vector<vtkSmartPointer<ttkTriangulationFilter> > inputTriangulations_;\
      int RequestData(vtkInformation *request, \
        vtkInformationVector **inputVector, \
        vtkInformationVector *outputVector){\
        \
        if((int) inputTriangulations_.size() != GetNumberOfInputPorts()){\
          inputTriangulations_.resize(GetNumberOfInputPorts());\
          for(int i = 0; i < (int) inputTriangulations_.size(); i++){\
            inputTriangulations_[i] = \
              vtkSmartPointer<ttkTriangulationFilter>::New();\
          }\
        }\
        \
        vector<vtkDataSet *> inputs(GetNumberOfInputPorts(), NULL);\
        vector<vtkDataSet *> outputs(GetNumberOfOutputPorts(), NULL);\
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
          threadNumber_ = OsCall::getNumberOfCores();\
        }\
        Modified();\
      }\
    private:\
      int doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs);\
      \
      bool needsToAbort(){ return GetAbortExecute();};\
      \
      int updateProgress(const float &progress) { \
        UpdateProgress(progress);\
        return 0;\
      };\
    protected:\
      bool                  UseAllCores;\
      int                   ThreadNumber;\
      

class VTKFILTERSCORE_EXPORT ttkTriangulationFilter 
  : public vtkDataSetAlgorithm, public Wrapper{

  public:
      
    static ttkTriangulationFilter* New();
    
    vtkTypeMacro(ttkTriangulationFilter, vtkDataSetAlgorithm);
    
  protected:
   
    ttkTriangulationFilter();
    ~ttkTriangulationFilter(){};
    
    int RequestData(vtkInformation *request, 
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);
   
    int RequestDataObject(vtkInformation *request,
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);
    
  private:
    
    bool needsToAbort();
    
    int updateProgress(const float &progress);
   
};

/// @}
