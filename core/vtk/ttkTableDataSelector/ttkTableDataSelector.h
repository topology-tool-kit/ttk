/// \ingroup vtk
/// \class ttkTableDataSelector
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date September 2018
///
/// \brief TTK VTK-filter that selects scalar fields on input with shallow copy.
///
/// \param Input Input scalar field (vtkTable)
/// \param Output Output scalar field (vtkTable)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
#pragma once

#include <limits>

// VTK includes
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataArraySelection.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkTableAlgorithm.h>

// VTK Module
#include <ttkTableDataSelectorModule.h>

// ttk code includes
#include <Wrapper.h>

class TTKTABLEDATASELECTOR_EXPORT ttkTableDataSelector
  : public vtkTableAlgorithm,
    protected ttk::Wrapper {

public:
  static ttkTableDataSelector *New();
  vtkTypeMacro(ttkTableDataSelector, vtkTableAlgorithm);

  // default ttk setters
  void SetDebugLevel(int debugLevel) {
    setDebugLevel(debugLevel);
    Modified();
  }

  vtkSetMacro(RegexpString, std::string);

  vtkGetVector2Macro(RangeId, int);
  vtkSetVector2Macro(RangeId, int);

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

  void AddCol(std::string s) {
    SelectedCols.push_back(s);
    Modified();
  }

  void ClearCols() {
    SelectedCols.clear();
    Modified();
  }

  vtkDataArraySelection *GetRangeIds() {
    vtkDataArraySelection *arr = vtkDataArraySelection::New();
    arr->SetArraySetting("0", true);
    arr->SetArraySetting(
      std::to_string(AvailableCols.size() - 1).c_str(), true);
    return arr;
  }

protected:
  ttkTableDataSelector() {
    UseAllCores = true;

    RegexpString = ".*";

    RangeId[0] = 0;
    RangeId[1] = std::numeric_limits<int>::max();
  }

  ~ttkTableDataSelector() override{};

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  void FillAvailableCols(vtkTable *input);

private:
  bool UseAllCores;
  int ThreadNumber;
  std::vector<std::string> SelectedCols;
  std::vector<std::string> AvailableCols;
  std::string RegexpString;
  int RangeId[2];

  int doIt(vtkTable *input, vtkTable *output);
  bool needsToAbort() override;
  int updateProgress(const float &progress) override;
};
