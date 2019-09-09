/// \ingroup vtk
/// \class ttkCinemaQuery
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that uses a SQL statement to select a subset of a
/// vtkTable.
///
/// This filter creates a temporary SQLite3 database from the input table,
/// performs a SQL query, and then returns the result as a vtkTable.
///
/// VTK wrapping code for the @CinemaQuery package.
///
/// \param Input Input table (vtkTable)
/// \param Output Output table (vtkTable)
///
/// sa ttk::CinemaQuery

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkTableAlgorithm.h>

// TTK includes
#include <CinemaQuery.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCinemaQuery
#else
class ttkCinemaQuery
#endif
  : public vtkTableAlgorithm,
    public ttk::Wrapper {

public:
  static ttkCinemaQuery *New();
  vtkTypeMacro(ttkCinemaQuery, vtkTableAlgorithm)

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

  vtkSetMacro(QueryString, std::string);
  vtkGetMacro(QueryString, std::string);

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
  ttkCinemaQuery() {
    QueryString = "";
    UseAllCores = false;

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }
  ~ttkCinemaQuery(){};

  bool UseAllCores;
  int ThreadNumber;

  std::string QueryString;
  ttk::CinemaQuery cinemaQuery;

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