/// \ingroup vtk
/// \class ttkA3Refinement
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the a3Refinement processing package.
///
/// VTK wrapping code for the @A3Refinement package.
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
/// \sa ttk::A3Refinement
#pragma once

// VTK includes -- to adapt
#include                  <vtkHyperTreeGridAlgorithm.h>
#include                  <vtkInformation.h>
#include                  <vtkSmartPointer.h>

// ttk code includes
#include                  <A3Refinement.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkA3Refinement
#else
class ttkA3Refinement
#endif
  : public vtkHyperTreeGridAlgorithm, public ttk::Wrapper{

  public:

    static ttkA3Refinement* New();
    vtkTypeMacro(ttkA3Refinement, vtkHyperTreeGridAlgorithm)

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else
        threadNumber_ = ttk::OsCall::getNumberOfCores();
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

    vtkSetMacro(SomeIntegerArgument, int);
    vtkGetMacro(SomeIntegerArgument, int);

    vtkSetMacro(SomeDoubleArgument, double);
    vtkGetMacro(SomeDoubleArgument, double);

    vtkSetMacro(SomeOption, bool);
    vtkGetMacro(SomeOption, bool);

    vtkSetMacro(ScalarField, std::string);
    vtkGetMacro(ScalarField, std::string);


    int FillInputPortInformation(int port, vtkInformation *info) override {
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
          break;
        default:
          break;
      }

      return 1;
    }

    int FillOutputPortInformation(int port, vtkInformation *info) override {
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkHyperTreeGrid");
          break;
        default:
          break;
      }

      return 1;
    }


  protected:

    ttkA3Refinement(){
        // init
      SomeIntegerArgument = 1;
      SomeDoubleArgument = 1;
      SomeOption = true;
      outputScalarField_ = NULL;

      UseAllCores = true;
      ThreadNumber = 1;
      debugLevel_ = 3;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }

    ~ttkA3Refinement(){};

    int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;
    int ProcessTrees(vtkHyperTreeGrid*, vtkDataObject*) override;

    bool UseAllCores;
    int ThreadNumber;


  private:

    int                   SomeIntegerArgument;
    double                SomeDoubleArgument;
    bool                  SomeOption;
    std::string           ScalarField;
    vtkDataArray          *outputScalarField_;
    ttk::A3Refinement            a3Refinement_;

    bool needsToAbort() override { return GetAbortExecute();};
    int updateProgress(const float &progress) override {
        UpdateProgress(progress);
        return 0;
    };
};
