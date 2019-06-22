/// \ingroup vtk
/// \class ttkIdentifyByScalarField
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the identifyByScalarField processing
/// package.
///
/// VTK wrapping code for the @IdentifyByScalarField package.
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
/// \sa ttk::IdentifyByScalarField
#pragma once

#include <vtkCellData.h>
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
#include <vtkSmartPointer.h>

// ttk code includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkIdentifyByScalarField
#else
class ttkIdentifyByScalarField
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkIdentifyByScalarField *New();
  vtkTypeMacro(ttkIdentifyByScalarField, vtkDataSetAlgorithm)

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

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(IncreasingOrder, bool);
  vtkGetMacro(IncreasingOrder, bool);

  vtkSetMacro(StartByOne, bool);
  vtkGetMacro(StartByOne, bool);

  int FillInputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
        break;
      default:
        break;
    }

    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
        break;
      default:
        break;
    }

    return 1;
  }

  int getScalars(vtkDataSet *input);

protected:
  ttkIdentifyByScalarField() {
    UseAllCores = true;
    IncreasingOrder = false;
    StartByOne = false;
    ScalarFieldId = 0;
    inputScalars_ = nullptr;
    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }

  ~ttkIdentifyByScalarField(){};

  TTK_SETUP();

private:
  int ScalarFieldId;
  bool IncreasingOrder;
  bool StartByOne;
  std::string ScalarField;

  vtkDataArray *inputScalars_;
};
