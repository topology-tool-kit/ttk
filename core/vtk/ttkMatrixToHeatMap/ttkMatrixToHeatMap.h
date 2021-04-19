/// \class ttkMatrixToHeatMap
/// \ingroup vtk
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \brief Generates a Heat Map from a Distance Matrix stored into a vtkTable
///
/// \param Input vtkTable containing a distance matrix
/// \param Output Heat map as a vtkUnstructuredGrid
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
#pragma once

#include <ttkMatrixToHeatMapModule.h>

#include <ttkAlgorithm.h>

#include <string>
#include <vector>

class TTKMATRIXTOHEATMAP_EXPORT ttkMatrixToHeatMap : public ttkAlgorithm {
public:
  static ttkMatrixToHeatMap *New();
  vtkTypeMacro(ttkMatrixToHeatMap, ttkAlgorithm);

  void SetScalarFields(const std::string &s) {
    ScalarFields.emplace_back(s);
    Modified();
  }
  void ClearScalarFields() {
    ScalarFields.clear();
    Modified();
  }

  vtkSetMacro(SelectFieldsWithRegexp, bool);
  vtkGetMacro(SelectFieldsWithRegexp, bool);

  vtkSetMacro(RegexpString, const std::string &);
  vtkGetMacro(RegexpString, std::string);

protected:
  ttkMatrixToHeatMap();
  ~ttkMatrixToHeatMap() = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool SelectFieldsWithRegexp{false};
  std::string RegexpString{".*"};
  std::vector<std::string> ScalarFields;
};
