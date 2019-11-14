/// \ingroup vtk
/// \class ttkCinemaWriter
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that writes input to disk.
///
/// This filter stores the input as a VTK dataset to disk and updates the
/// data.csv file of a Cinema Spec D database.
///
/// \param Input vtkDataSet to be stored (vtkDataSet)

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkXMLPMultiBlockDataWriter.h>

// TTK includes
#include <ttkTopologicalCompressionWriter.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCinemaWriter
#else
class ttkCinemaWriter
#endif
  : public vtkXMLPMultiBlockDataWriter,
    public ttk::Wrapper {

public:
  static ttkCinemaWriter *New();
  vtkTypeMacro(ttkCinemaWriter, vtkXMLPMultiBlockDataWriter)

    vtkSetMacro(DatabasePath, std::string);
  vtkGetMacro(DatabasePath, std::string);

  vtkSetMacro(OverrideDatabase, bool);
  vtkGetMacro(OverrideDatabase, bool);

  vtkSetMacro(CompressLevel, int);
  vtkGetMacro(CompressLevel, int);

  vtkSetMacro(UseTopologicalCompression, bool);
  vtkGetMacro(UseTopologicalCompression, bool);

#define TopoCompWriterGetSetMacro(NAME, TYPE) \
  void Set##NAME(const TYPE _arg) {           \
    this->ttkCompWriter_->Set##NAME(_arg);    \
    this->Modified();                         \
  }                                           \
  TYPE Get##NAME() {                          \
    return this->ttkCompWriter_->Get##NAME(); \
  }

  TopoCompWriterGetSetMacro(ScalarField, std::string);
  TopoCompWriterGetSetMacro(Tolerance, double);
  TopoCompWriterGetSetMacro(MaximumError, double);
  TopoCompWriterGetSetMacro(ZFPBitBudget, double);
  TopoCompWriterGetSetMacro(ZFPOnly, bool);
  TopoCompWriterGetSetMacro(CompressionType, int);
  TopoCompWriterGetSetMacro(Subdivide, bool);
  TopoCompWriterGetSetMacro(UseTopologicalSimplification, bool);

  void SetSQMethodPV(const int arg) {
    this->ttkCompWriter_->SetSQMethodPV(arg);
  }

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
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
        break;
      default:
        return 0;
    }
    return 1;
  }

protected:
  ttkCinemaWriter() {
    SetDatabasePath("");
    SetOverrideDatabase(true);
    SetCompressLevel(9);
    SetUseTopologicalCompression(false);

    UseAllCores = false;

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }
  ~ttkCinemaWriter(){};

  bool UseAllCores;
  int ThreadNumber;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string DatabasePath;
  bool OverrideDatabase;
  int CompressLevel;
  bool UseTopologicalCompression;
  vtkNew<ttkTopologicalCompressionWriter> ttkCompWriter_;

  bool needsToAbort() override {
    return GetAbortExecute();
  };
  int updateProgress(const float &progress) override {
    UpdateProgress(progress);
    return 0;
  };
};
