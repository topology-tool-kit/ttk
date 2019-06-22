/// \ingroup vtk
/// \class ttkAddFieldData
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.10.2018
///
/// \brief TTK VTK-filter that adds field data to a vtkDataObject.
///
/// This filter adds field data to a vtkDataObject based on a string and/or
/// point/cell/field data of an optional vtkDataObject.
///
/// VTK wrapping code for the @AddFieldData package.
///
/// \param Output vtkTable

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkPassInputTypeAlgorithm.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkAddFieldData
#else
class ttkAddFieldData
#endif
  : public vtkPassInputTypeAlgorithm,
    public ttk::Wrapper {

public:
  static ttkAddFieldData *New();
  vtkTypeMacro(ttkAddFieldData, vtkPassInputTypeAlgorithm)

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

  vtkSetMacro(FieldDataString, std::string);
  vtkGetMacro(FieldDataString, std::string);

  /**
   * Adds an array to pass through.
   * fieldType where the array that should be passed (point data, cell data,
   * etc.). It should be one of the constants defined in the
   * vtkDataObject::AttributeTypes enumeration.
   */
  virtual void AddArray(int fieldType, const char *name);
  virtual void AddPointDataArray(const char *name);
  virtual void AddCellDataArray(const char *name);
  virtual void AddFieldDataArray(const char *name);

  virtual void RemoveArray(int fieldType, const char *name, bool deleteType);
  virtual void RemovePointDataArray(const char *name);
  virtual void RemoveCellDataArray(const char *name);
  virtual void RemoveFieldDataArray(const char *name);

  /**
   * Clear all arrays to pass through.
   */
  virtual void ClearArrays();
  virtual void ClearPointDataArrays();
  virtual void ClearCellDataArrays();
  virtual void ClearFieldDataArrays();

  int FillInputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
        break;
      case 1:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
        info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
        break;
      default:
        return 0;
    }
    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
        break;
      default:
        return 0;
    }
    return 1;
  }

protected:
  ttkAddFieldData() {
    FieldDataString = "";
    UseAllCores = false;

    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(1);
  }
  ~ttkAddFieldData(){};

  bool UseAllCores;
  int ThreadNumber;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string FieldDataString;
  std::vector<std::pair<int, std::string>> ArraySelection;

  bool needsToAbort() override {
    return GetAbortExecute();
  };
  int updateProgress(const float &progress) override {
    UpdateProgress(progress);
    return 0;
  };
};