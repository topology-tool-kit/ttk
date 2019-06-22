/// \ingroup vtk
/// \class ttkCinemaProductReader
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that reads the data products that are referenced in a
/// vtkTable.
///
/// This filter reads the products that are referenced in a vtkTable. The
/// results are stored in a vtkMultiBlockDataSet where each block corresponds to
/// a row of the table with consistent ordering.
///
/// \param Input vtkTable that contains data product references (vtkTable)
/// \param Output vtkMultiBlockDataSet where each block is a referenced product
/// of an input table row (vtkMultiBlockDataSet)

#pragma once

// VTK includes
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiBlockDataSetAlgorithm.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCinemaProductReader
#else
class ttkCinemaProductReader
#endif
  : public vtkMultiBlockDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkCinemaProductReader *New();
  vtkTypeMacro(ttkCinemaProductReader, vtkMultiBlockDataSetAlgorithm)

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

  void SetFilepathColumnName(
    int idx, int port, int connection, int fieldAssociation, const char *name) {
    this->FilepathColumnName = std::string(name);
    this->Modified();
  };

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
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
        break;
      default:
        return 0;
    }
    return 1;
  }

protected:
  ttkCinemaProductReader() {
    UseAllCores = false;

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }
  ~ttkCinemaProductReader(){};

  bool UseAllCores;
  int ThreadNumber;

  std::string FilepathColumnName;

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