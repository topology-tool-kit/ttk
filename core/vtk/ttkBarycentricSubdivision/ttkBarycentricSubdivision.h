/// \ingroup vtk
/// \class ttkBarycentricSubdivision
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the barycentricSubdivision processing
/// package.
///
/// VTK wrapping code for the @BarycentricSubdivision package.
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
/// \sa ttk::BarycentricSubdivision
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
#include <vtkCellData.h>
#include <vtkSmartPointer.h>

// TTK code includes
#include <BarycentricSubdivision.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkBarycentricSubdivision
#else
class ttkBarycentricSubdivision
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkBarycentricSubdivision *New();
  vtkTypeMacro(ttkBarycentricSubdivision, vtkDataSetAlgorithm);

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

  int FillInputPortInformation(int /*port*/, vtkInformation *info) override {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }

protected:
  ttkBarycentricSubdivision() {
    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }

  TTK_SETUP();

private:
  // base worker
  ttk::BarycentricSubdivision baseWorker_{};
};
