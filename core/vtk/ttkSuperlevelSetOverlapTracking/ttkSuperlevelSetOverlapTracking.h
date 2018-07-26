/// ingroup vtk
/// class ttkSuperlevelSetOverlapTracking
/// author Your Name Here <Your Email Address Here>
/// date The Date Here.
///
/// brief TTK VTK-filter that wraps the superlevelSetOverlapTracking processing package.
///
/// VTK wrapping code for the @SuperlevelSetOverlapTracking package.
///
/// param Input Input scalar field (vtkDataSet)
/// param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// sa ttk::SuperlevelSetOverlapTracking
#pragma once

// VTK includes -- to adapt
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkTable.h>

// ttk code includes
#include                  <SuperlevelSetOverlapTracking.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkSuperlevelSetOverlapTracking
#else
class ttkSuperlevelSetOverlapTracking
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkSuperlevelSetOverlapTracking* New();
    vtkTypeMacro(ttkSuperlevelSetOverlapTracking, vtkDataSetAlgorithm)

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = ttk::OsCall::getNumberOfCores();
      }
      Modified();
    }

    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    vtkSetMacro(Level, double);
    vtkGetMacro(Level, double);

    int FillInputPortInformation(int port, vtkInformation *info) override {

      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
          break;
        default:
          break;
      }

      return 1;
    }

    int FillOutputPortInformation(int port, vtkInformation *info) override {
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
          break;
        default:
          break;
      }

      return 1;
    }


  protected:
    ttkSuperlevelSetOverlapTracking(){
      Level = 0;
      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }
    ~ttkSuperlevelSetOverlapTracking(){};

    int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

    int RequestDataObject(vtkInformation *request,
      vtkInformationVector **inputVector, vtkInformationVector *outputVector) override {

      vtkDataSetAlgorithm::RequestDataObject(request, inputVector, outputVector);

      for(int i = 0; i < GetNumberOfOutputPorts(); i++){

        vtkInformation *outputInformation =
          outputVector->GetInformationObject(i);

        vtkDataSet *output = vtkDataSet::SafeDownCast(
          outputInformation->Get(vtkDataObject::DATA_OBJECT()));

        if(output){

          std::string dataType = output->GetClassName();

          TTK_UNSTRUCTURED_GRID_NEW(i, outputInformation, dataType);

          TTK_POLY_DATA_NEW(i, outputInformation, dataType);
        }
      }

      return 1;
    }

    bool UseAllCores;
    int ThreadNumber;
    // TTK_SETUP();
    // TTK_OUTPUT_MANAGEMENT();

private:

    double                Level;
    ttk::SuperlevelSetOverlapTracking superlevelSetOverlapTracking_;

    // base code features
    int doIt(vtkTable* inputs, vtkUnstructuredGrid* outputs);

    bool needsToAbort() override { return GetAbortExecute();};
    int updateProgress(const float &progress) override {
        UpdateProgress(progress);
        return 0;
    };

};
