/// \ingroup vtk
/// \class ttkDataSetInterpolator
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the dataSetInterpolator processing package.
///
/// VTK wrapping code for the @DataSetInterpolator package.
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
/// \sa ttk::DataSetInterpolator
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
#include <vtkProbeFilter.h>
#include <vtkSmartPointer.h>

#include <Wrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkDataSetInterpolator
#else
class ttkDataSetInterpolator
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkDataSetInterpolator *New();
  vtkTypeMacro(ttkDataSetInterpolator, vtkDataSetAlgorithm)

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

  void SetThreads() {
    if(!UseAllCores)
      threadNumber_ = ThreadNumber;
    else {
      threadNumber_ = ttk::OsCall::getNumberOfCores();
    }
    Modified();
  }

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

protected:
  ttkDataSetInterpolator() {
    UseAllCores = true;

    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(1);
  }

  ~ttkDataSetInterpolator(){};

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool UseAllCores;
  int ThreadNumber;

  int doIt(vtkDataSet *source, vtkDataSet *target, vtkDataSet *output);
  bool needsToAbort() override;
  int updateProgress(const float &progress) override;
};
