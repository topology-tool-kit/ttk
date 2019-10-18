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

  void SetTolerance(double arg) {
    this->ttkCompWriter_->SetTolerance(arg);
  }
  double GetTolerance() {
    return this->ttkCompWriter_->GetTolerance();
  }
  void SetMaximumError(double arg) {
    this->ttkCompWriter_->SetMaximumError(arg);
  }
  double GetMaximumError() {
    return this->ttkCompWriter_->GetMaximumError();
  }
  void SetZFPBitBudget(double arg) {
    this->ttkCompWriter_->SetZFPBitBudget(arg);
  }
  double GetZFPBitBudget() {
    return this->ttkCompWriter_->GetZFPBitBudget();
  }
  void SetZFPOnly(bool arg) {
    this->ttkCompWriter_->SetZFPOnly(arg);
  }
  bool GetZFPOnly() {
    return this->ttkCompWriter_->GetZFPOnly();
  }
  void SetCompressionType(int arg) {
    this->ttkCompWriter_->SetCompressionType(arg);
  }
  int GetCompressionType() {
    return this->ttkCompWriter_->GetCompressionType();
  }
  void SetNbSegments(int arg) {
    this->ttkCompWriter_->SetNbSegments(arg);
  }
  int GetNbSegments() {
    return this->ttkCompWriter_->GetNbSegments();
  }
  void SetNbVertices(int arg) {
    this->ttkCompWriter_->SetNbVertices(arg);
  }
  int GetNbVertices() {
    return this->ttkCompWriter_->GetNbVertices();
  }
  void SetSQMethod(std::string arg) {
    this->ttkCompWriter_->SetSQMethod(arg);
  }
  std::string GetSQMethod() {
    return this->ttkCompWriter_->GetSQMethod();
  }
  void SetSubdivide(bool arg) {
    this->ttkCompWriter_->SetSubdivide(arg);
  }
  bool GetSubdivide() {
    return this->ttkCompWriter_->GetSubdivide();
  }
  void SetUseTopologicalSimplification(bool arg) {
    this->ttkCompWriter_->SetUseTopologicalSimplification(arg);
  }
  bool GetUseTopologicalSimplification() {
    return this->ttkCompWriter_->GetUseTopologicalSimplification();
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
