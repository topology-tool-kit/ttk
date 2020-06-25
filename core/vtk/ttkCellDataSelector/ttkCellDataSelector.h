/// \ingroup vtk
/// \class ttkCellDataSelector
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date December 2017
///
/// \brief TTK VTK-filter that selects scalar fields on input with shallow copy.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
#pragma once

#include <limits>

// VTK includes
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataArraySelection.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedShortArray.h>

// VTK Module
#include <ttkCellDataSelectorModule.h>

// ttk code includes
#include <Wrapper.h>

class TTKCELLDATASELECTOR_EXPORT ttkCellDataSelector
  : public vtkDataSetAlgorithm,
    protected ttk::Wrapper {

public:
  static ttkCellDataSelector *New();
  vtkTypeMacro(ttkCellDataSelector, vtkDataSetAlgorithm);

  // default ttk setters
  void SetDebugLevel(int debugLevel) {
    setDebugLevel(debugLevel);
    Modified();
  }

  vtkSetMacro(RegexpString, std::string);

  vtkSetVector2Macro(RangeId, int);
  vtkGetVector2Macro(RangeId, int);

  vtkSetMacro(RenameSelected, bool);
  vtkGetMacro(RenameSelected, bool);

  vtkSetMacro(SelectedFieldName, std::string);
  vtkGetMacro(SelectedFieldName, std::string);

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

  void AddScalarField(std::string s) {
    SelectedFields.push_back(s);
    Modified();
  }

  void ClearScalarFields() {
    SelectedFields.clear();
    Modified();
  }

  vtkDataArraySelection *GetRangeIds() {
    vtkDataArraySelection *arr = vtkDataArraySelection::New();
    arr->SetArraySetting("0", true);
    arr->SetArraySetting(
      std::to_string(AvailableFields.size() - 1).c_str(), true);
    return arr;
  }

protected:
  ttkCellDataSelector() {
    UseAllCores = false;
    RenameSelected = false;

    RegexpString = ".*";
    SelectedFieldName = "SelectedField";

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);

    localFieldCopy_ = NULL;

    RangeId[0] = 0;
    RangeId[1] = std::numeric_limits<int>::max();
  }

  ~ttkCellDataSelector() override {
    if(localFieldCopy_)
      localFieldCopy_->Delete();
  };

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  void FillAvailableFields(vtkDataSet *input);

private:
  bool UseAllCores;
  int ThreadNumber;
  bool RenameSelected;
  std::string SelectedFieldName;
  std::vector<std::string> AvailableFields;
  std::vector<std::string> SelectedFields;
  std::string RegexpString;
  vtkDataArray *localFieldCopy_;
  int RangeId[2];

  int doIt(vtkDataSet *input, vtkDataSet *output);
  bool needsToAbort() override;
  int updateProgress(const float &progress) override;
};
