/// \ingroup vtk
/// \class ttkCinemaLayout
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.10.2018
///
/// \brief TTK VTK-filter that arranges vtkDataSets on a grid.
///
/// This filter arranges vtkDataSets which are stored in a vtkMultiBlockDataSet
/// on a grid.
///
/// \param Input vtkMultiBlockDataSet
/// \param Output vtkMultiBlockDataSet

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkXMLPMultiBlockDataWriter.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCinemaLayout
#else
class ttkCinemaLayout
#endif
  : public vtkXMLPMultiBlockDataWriter,
    public ttk::Wrapper {

public:
  static ttkCinemaLayout *New();
  vtkTypeMacro(ttkCinemaLayout, vtkXMLPMultiBlockDataWriter)

    vtkSetMacro(RowAxis, int);
  vtkGetMacro(RowAxis, int);

  vtkSetMacro(RowGap, double);
  vtkGetMacro(RowGap, double);

  vtkSetMacro(ColumnAxis, int);
  vtkGetMacro(ColumnAxis, int);

  vtkSetMacro(ColumnGap, double);
  vtkGetMacro(ColumnGap, double);

  vtkSetMacro(NumberOfRows, int);
  vtkGetMacro(NumberOfRows, int);

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
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
        break;
      default:
        return 0;
    }
    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
        break;
      default:
        return 0;
    }
    return 1;
  }

protected:
  ttkCinemaLayout() {
    SetRowAxis(0);
    SetRowGap(0);
    SetColumnAxis(0);
    SetColumnGap(0);
    SetNumberOfRows(0);

    UseAllCores = false;

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }
  ~ttkCinemaLayout(){};

  int RowAxis;
  double RowGap;

  int ColumnAxis;
  double ColumnGap;

  int NumberOfRows;

  bool UseAllCores;
  int ThreadNumber;

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
