/// ingroup vtk
/// class ttkCinemaQuery
/// author Your Name Here <Your Email Address Here>
/// date The Date Here.
///
/// brief TTK VTK-filter that wraps the cinemaQuery processing package.
///
/// VTK wrapping code for the @CinemaQuery package.
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
/// sa ttk::CinemaQuery
#pragma once

// VTK includes -- to adapt
#include                <vtkCharArray.h>
#include                <vtkDataArray.h>
#include                <vtkDataSet.h>
#include                <vtkTableAlgorithm.h>
#include                <vtkDoubleArray.h>
#include                <vtkFiltersCoreModule.h>
#include                <vtkFloatArray.h>
#include                <vtkInformation.h>
#include                <vtkIntArray.h>
#include                <vtkObjectFactory.h>
#include                <vtkPointData.h>
#include                <vtkSmartPointer.h>
#include                <vtkTable.h>

// ttk code includes
#include                  <CinemaQuery.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCinemaQuery
#else
class ttkCinemaQuery
#endif
  : public vtkTableAlgorithm, public ttk::Wrapper{

  public:

    static ttkCinemaQuery* New();
    vtkTypeMacro(ttkCinemaQuery, vtkTableAlgorithm)

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

    vtkSetMacro(QueryString, std::string);
    vtkGetMacro(QueryString, std::string);

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
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
          break;
        default:
          break;
      }

      return 1;
    }

  protected:

    ttkCinemaQuery(){
      QueryString = "";
      UseAllCores = false;

      // Specify the number of input and output ports.
      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }
    ~ttkCinemaQuery(){};

    int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

    bool UseAllCores;
    int ThreadNumber;
    // TTK_SETUP();

  private:

    std::string      QueryString;
    ttk::CinemaQuery cinemaQuery_;

    // base code features
    int doIt(vtkTable* input, vtkTable* output);

    bool needsToAbort() override { return GetAbortExecute();};
    int updateProgress(const float &progress) override {
        UpdateProgress(progress);
        return 0;
    };
};
