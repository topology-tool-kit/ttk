#include <ttkMatrixToHeatMap.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkFiltersCoreModule.h>
#include <vtkObjectFactory.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

#include <regex>

vtkStandardNewMacro(ttkMatrixToHeatMap);

ttkMatrixToHeatMap::ttkMatrixToHeatMap() {
  this->setDebugMsgPrefix("MatrixToHeatMap");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkMatrixToHeatMap::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkMatrixToHeatMap::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkMatrixToHeatMap::RequestData(vtkInformation * /*request*/,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {

  const auto input = vtkTable::GetData(inputVector[0], 0);
  const auto output = vtkUnstructuredGrid::GetData(outputVector, 0);

  if(SelectFieldsWithRegexp) {
    // select all input columns whose name is matching the regexp
    ScalarFields.clear();
    const auto n = input->GetNumberOfColumns();
    for(int i = 0; i < n; ++i) {
      const auto &name = input->GetColumnName(i);
      if(std::regex_match(name, std::regex(RegexpString))) {
        ScalarFields.emplace_back(name);
      }
    }
  }

  return 1;
}
