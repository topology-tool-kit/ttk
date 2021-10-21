#include <ttkPointSetToCurve.h>

#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkUnstructuredGrid.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <ttkUtils.h>

#include <array>
#include <map>
#include <set>

vtkStandardNewMacro(ttkPointSetToCurve);

ttkPointSetToCurve::ttkPointSetToCurve() {
  this->setDebugMsgPrefix("PointSetToCurve");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkPointSetToCurve::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkPointSetToCurve::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

template <typename VTK_TT>
void ttkPointSetToCurve::dispatch(
  std::vector<std::pair<vtkIdType, double>> &storage,
  const VTK_TT *const values,
  const size_t nvalues) {

  for(size_t i = 0; i < nvalues; ++i) {
    storage.emplace_back(i, static_cast<double>(values[i]));
  }
}

int ttkPointSetToCurve::RequestData(vtkInformation *ttkNotUsed(request),
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  const auto input = vtkPointSet::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  if(input == nullptr || output == nullptr) {
    this->printErr("Null input data, aborting");
    return 0;
  }

  // ordering array
  const auto oa = this->GetInputArrayToProcess(0, inputVector);

  if(oa == nullptr) {
    this->printErr("Cannot find the required data array");
    return 0;
  }

  const auto nvalues = oa->GetNumberOfTuples();

  // store point index <-> ordering value in vector
  std::vector<std::pair<vtkIdType, double>> orderedValues{};

  switch(oa->GetDataType()) {
    vtkTemplateMacro(
      dispatch(orderedValues,
               static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(oa)), nvalues));
  }

  // compare two pairs of index/value according to their values
  const auto cmp
    = [](const std::pair<vtkIdType, double> &a,
         const std::pair<vtkIdType, double> &b) { return a.second < b.second; };

  // sort the vector of indices/values in ascending order
  std::sort(orderedValues.begin(), orderedValues.end(), cmp);

  // shallow-copy input into output
  output->ShallowCopy(input);

  // output cell array
  vtkNew<vtkCellArray> cells{};

  for(size_t i = 1; i < orderedValues.size(); ++i) {
    std::array<vtkIdType, 2> linePoints{
      orderedValues[i - 1].first, orderedValues[i].first};
    cells->InsertNextCell(2, linePoints.data());
  }

  if(CloseCurve) {
    std::array<vtkIdType, 2> linePoints{
      orderedValues.back().first, orderedValues.front().first};
    cells->InsertNextCell(2, linePoints.data());
  }

  output->SetCells(VTK_LINE, cells);

  return 1;
}
