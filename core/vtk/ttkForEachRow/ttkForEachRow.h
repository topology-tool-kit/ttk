/// \ingroup vtk
/// \class ttkForEachRow
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2018
///
/// \brief TTK VTK-filter that iterates over rows of a vtkTable.
///
/// This filter works in conjunction with the ttkEndFor filter to iterate over
/// rows of a vtkTable.
///
/// \param Input vtkTable table that will be iterated over
/// \param Output vtkTable table that contains only one row of the input

#pragma once

// VTK includes
#include <vtkMultiBlockDataSetAlgorithm.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkForEachRow
#else
class ttkForEachRow
#endif
  : public vtkMultiBlockDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkForEachRow *New();
  vtkTypeMacro(ttkForEachRow, vtkMultiBlockDataSetAlgorithm)

    // default ttk setters
    vtkSetMacro(debugLevel_, int);
  void SetThreads() {
    threadNumber_
      = !UseAllCores ? ThreadNumber : ttk::OsCall::getNumberOfCores();
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

  int FillInputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
        break;
      default:
        return 0;
    }
    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
        break;
      default:
        return 0;
    }
    return 1;
  }

protected:
  ttkForEachRow() {
    UseAllCores = false;

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }
  ~ttkForEachRow(){};

  bool UseAllCores;
  int ThreadNumber;

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool needsToAbort() override {
    return GetAbortExecute();
  };
  int updateProgress(const float &progress) override {
    UpdateProgress(progress);
    return 0;
  };
};