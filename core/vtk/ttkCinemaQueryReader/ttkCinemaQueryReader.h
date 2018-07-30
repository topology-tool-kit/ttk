/// ingroup vtk
/// class ttkCinemaQueryLoader
/// author Your Name Here <Your Email Address Here>
/// date The Date Here.
///
/// brief TTK VTK-filter that wraps the cinemaQueryLoader processing package.
///
/// VTK wrapping code for the @CinemaQueryLoader package.
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
/// sa ttk::CinemaQueryLoader
#pragma once

#include                  <vtkTable.h>
#include                  <vtkImageData.h>

#include                  <vtkMultiBlockDataSetAlgorithm.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkInformation.h>
#include                  <vtkMultiBlockDataSet.h>

// ttk code includes
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkImageAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCinemaQueryLoader
#else
class ttkCinemaQueryLoader
#endif
  : public vtkMultiBlockDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkCinemaQueryLoader* New();
    vtkTypeMacro(ttkCinemaQueryLoader, vtkMultiBlockDataSetAlgorithm)

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

    vtkSetMacro(Index, int);
    vtkGetMacro(Index, int);

    vtkSetMacro(NumberOfRows, int);
    vtkGetMacro(NumberOfRows, int);

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
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
          break;
        default:
          break;
      }

      return 1;
    }


  protected:

    ttkCinemaQueryLoader(){
      Index = 0;
      NumberOfRows = 1;

      UseAllCores = false;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }

    ~ttkCinemaQueryLoader(){};

    // TTK_SETUP();
    int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;
    // int RequestInformation (vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

    bool UseAllCores;
    int ThreadNumber;

  private:

    int                    Index;
    int                    NumberOfRows;

    int doIt(vtkTable* inputTable, vtkImageData* output);

    bool needsToAbort() override { return GetAbortExecute();};
    int updateProgress(const float &progress) override {
        UpdateProgress(progress);
        return 0;
    };

};
