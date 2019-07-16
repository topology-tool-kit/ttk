/// \ingroup vtk
/// \class ttkBlank
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the blank processing package.
///
/// VTK wrapping code for the @Blank package.
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
/// \sa ttk::Blank
#pragma once

// VTK includes -- to adapt
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>

// ttk code includes
#include <Blank.h>
#include <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkBlank
#else
class ttkBlank
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkBlank *New();
  vtkTypeMacro(ttkBlank, vtkDataSetAlgorithm)

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }
  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  // TODO-4
  // set-getters macros to define from each variable you want to access from
  // the outside (in particular from paraview) - to adapt.
  // Note that the XML file for the ParaView plug-in specification needs to be
  // edited accordingly.
  vtkSetMacro(SomeIntegerArgument, int);
  vtkGetMacro(SomeIntegerArgument, int);

  vtkSetMacro(SomeDoubleArgument, double);
  vtkGetMacro(SomeDoubleArgument, double);

  vtkSetMacro(SomeOption, bool);
  vtkGetMacro(SomeOption, bool);

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);
  // end of TODO-4

  // TODO-2
  // Over-ride the input types.
  // By default, this filter has one input and one output, of the same type.
  // Here, you can re-define the input types, on a per input basis.
  // In this example, the first input type is forced to vtkUnstructuredGrid.
  // The second input type is forced to vtkImageData.
  //     int FillInputPortInformation(int port, vtkInformation *info) override {
  //
  //       switch(port){
  //         case 0:
  //           info->Set(vtkDataObject::DATA_TYPE_NAME(),
  //           "vtkUnstructuredGrid"); break;
  //         case 1:
  //           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  //           break;
  //         default:
  //           break;
  //       }
  //
  //       return 1;
  //     }
  // end of TODO-2

  // TODO-3
  // Over-ride the output types.
  // By default, this filter has one input and one output, of the same type.
  // Here, you can re-define the output types, on a per output basis.
  // In this example, the first output type is forced to vtkUnstructuredGrid.
  // The second output type is forced to vtkImageData.
  //     int FillOutputPortInformation(int port, vtkInformation *info) override
  //     {
  //
  //       switch(port){
  //         case 0:
  //           info->Set(vtkDataObject::DATA_TYPE_NAME(),
  //           "vtkUnstructuredGrid"); break;
  //         case 1:
  //           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  //           break;
  //         default:
  //           break;
  //       }
  //
  //       return 1;
  //     }
  // end of TODO-3

protected:
  ttkBlank() {

    // init
    SomeIntegerArgument = 1;
    SomeDoubleArgument = 1;
    SomeOption = true;
    outputScalarField_ = NULL;

    UseAllCores = true;

    // TODO-1
    // Specify the number of input and output ports.
    // By default, this filter has one input and one output.
    // In this example, we define 2 inputs and 2 outputs.
    //       SetNumberOfInputPorts(2);
    //       SetNumberOfOutputPorts(2);
    // end of TODO-1
  }

  ~ttkBlank(){};

  TTK_SETUP();

private:
  int SomeIntegerArgument;
  double SomeDoubleArgument;
  bool SomeOption;
  std::string ScalarField;
  vtkDataArray *outputScalarField_;
  ttk::blank::Blank blank_;
};
