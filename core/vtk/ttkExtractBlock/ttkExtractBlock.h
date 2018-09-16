/// \ingroup vtk
/// \class ttkExtractBlock
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the extractBlock processing package.
///
/// VTK wrapping code for the @ExtractBlock package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::ExtractBlock
#pragma once

// VTK includes -- to adapt
#include                  <vtkMultiBlockDataSetAlgorithm.h>
#include                  <vtkInformation.h>

// ttk code includes
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkExtractBlock
#else
class ttkExtractBlock
#endif
  : public vtkMultiBlockDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkExtractBlock* New();
    vtkTypeMacro(ttkExtractBlock, vtkMultiBlockDataSetAlgorithm)

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

    int FillInputPortInformation(int port, vtkInformation *info) override {

      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
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

    ttkExtractBlock(){

        // init
      Index = 0;

      UseAllCores = true;
      ThreadNumber = 1;
      debugLevel_ = 3;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }
    ~ttkExtractBlock(){};

    int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

        bool UseAllCores;
    int ThreadNumber;

  private:

    int Index;

    bool needsToAbort() override { return GetAbortExecute();};
    int updateProgress(const float &progress) override {
        UpdateProgress(progress);
        return 0;
    };
};
