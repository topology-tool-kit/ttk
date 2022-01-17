#include <ttkRipsComplex.h>
#include <ttkUtils.h>

#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

#include <regex>

vtkStandardNewMacro(ttkRipsComplex);

ttkRipsComplex::ttkRipsComplex() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkRipsComplex::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkRipsComplex::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkRipsComplex::RequestData(vtkInformation *ttkNotUsed(request),
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {

  using ttk::SimplexId;

  auto *input = vtkTable::GetData(inputVector[0]);
  auto *output = vtkUnstructuredGrid::GetData(outputVector);

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
