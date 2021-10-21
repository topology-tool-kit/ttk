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
#include <vtkDataArray.h>
#include <vtkDataArraySelection.h>
#include <vtkDataSetAttributes.h>
#include <vtkInformation.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

// VTK Module
#include <ttkTableDataSelectorModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class TTKTABLEDATASELECTOR_EXPORT ttkTableDataSelector : public ttkAlgorithm {

public:
  static ttkTableDataSelector *New();
  vtkTypeMacro(ttkTableDataSelector, ttkAlgorithm);

  // default ttk setters

  vtkSetMacro(RegexpString, const std::string &);

  vtkGetVector2Macro(RangeId, int);
  vtkSetVector2Macro(RangeId, int);

  // end of default ttk setters

  void AddCol(const std::string &s) {
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
    this->setDebugMsgPrefix("TableDataSelector");

    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);

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

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  void FillAvailableCols(vtkTable *input);

private:
  std::vector<std::string> SelectedCols;
  std::vector<std::string> AvailableCols;
  std::string RegexpString;
  int RangeId[2];
};
