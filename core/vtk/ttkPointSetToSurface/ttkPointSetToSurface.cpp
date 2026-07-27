#include <ttkPointSetToSurface.h>

#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkPolyData.h>

#include <ttkUtils.h>

#include <array>

vtkStandardNewMacro(ttkPointSetToSurface);

ttkPointSetToSurface::ttkPointSetToSurface() {
  this->setDebugMsgPrefix("PointSetToSurface");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkPointSetToSurface::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkPointSetToSurface::FillOutputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

template <typename VTK_T1, typename VTK_T2>
void ttkPointSetToSurface::dispatch(
  std::vector<std::tuple<vtkIdType, double, double>> &storage,
  const VTK_T1 *const values,
  const VTK_T2 *const values2,
  const size_t nvalues) {

  for(size_t i = 0; i < nvalues; ++i) {
    storage.emplace_back(
      i, static_cast<double>(values[i]), static_cast<double>(values2[i]));
  }
}

int ttkPointSetToSurface::RequestData(vtkInformation *ttkNotUsed(request),
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {
  const auto input = vtkPointSet::GetData(inputVector[0]);
  auto output = vtkPolyData::GetData(outputVector);

  if(input == nullptr || output == nullptr) {
    this->printErr("Null input data, aborting");
    return 0;
  }
  const auto oa = this->GetInputArrayToProcess(0, inputVector);
  const auto oa2 = this->GetInputArrayToProcess(1, inputVector);
  if(oa == nullptr) {
    this->printErr("Cannot find the required data X array");
    return 0;
  }
  if(oa2 == nullptr) {
    this->printErr("Cannot find the required data Y array");
    return 0;
  }

  const auto nvalues = oa->GetNumberOfTuples();

  // store point index <-> ordering value in vector
  std::vector<std::tuple<vtkIdType, double, double>> orderedValues;

#ifndef TTK_ENABLE_DOUBLE_TEMPLATING
  switch(oa->GetDataType()) {
    vtkTemplateMacro(dispatch(
      orderedValues, static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(oa)),
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(oa2)), nvalues));
  }
#else
  switch(vtkTemplate2PackMacro(oa->GetDataType(), oa2->GetDataType())) {
    vtkTemplate2Macro(dispatch(
      orderedValues, static_cast<VTK_T1 *>(ttkUtils::GetVoidPointer(oa)),
      static_cast<VTK_T2 *>(ttkUtils::GetVoidPointer(oa2)), nvalues));
  }
#endif // TTK_ENABLE_DOUBLE_TEMPLATING

  // Get number of unique values of each array
  std::vector<double> xValues(orderedValues.size()),
    yValues(orderedValues.size());
  for(unsigned int i = 0; i < orderedValues.size(); ++i) {
    auto tup = orderedValues[i];
    xValues[i] = std::get<1>(tup);
    yValues[i] = std::get<2>(tup);
  }
  TTK_PSORT(this->threadNumber_, xValues.begin(), xValues.end());
  const auto nUniqueXValues
    = std::unique(xValues.begin(), xValues.end()) - xValues.begin();
  TTK_PSORT(this->threadNumber_, yValues.begin(), yValues.end());
  const auto nUniqueYValues
    = std::unique(yValues.begin(), yValues.end()) - yValues.begin();

  if(nUniqueXValues * nUniqueYValues != input->GetNumberOfPoints()) {
    printErr(
      "Number of unique values in first array times the number of unique "
      "values in the second one does not equal the number of points");
    return 0;
  }

  // Compare two pairs of index/value according to their values
  double yRange[2] = {*std::min_element(yValues.begin(), yValues.end()),
                      *std::max_element(yValues.begin(), yValues.end())};
  double minXInterval = std::numeric_limits<double>::max();
  for(unsigned int i = 1; i < nUniqueXValues; ++i)
    minXInterval = std::min(minXInterval, xValues[i] - xValues[i - 1]);
  const auto normValue = [&](double v) {
    return (v - yRange[0]) / (yRange[1] - yRange[0]) * minXInterval * 0.99;
  };
  const auto cmp = [&](const std::tuple<int, double, double> &a,
                       const std::tuple<int, double, double> &b) {
    return std::get<1>(a) + normValue(std::get<2>(a))
           < std::get<1>(b) + normValue(std::get<2>(b));
  };

  // sort the vector of indices/values in ascending order
  TTK_PSORT(
    this->threadNumber_, orderedValues.begin(), orderedValues.end(), cmp);

  // Create point ids matrix
  std::vector<std::vector<vtkIdType>> orderedIds(
    nUniqueXValues, std::vector<vtkIdType>(nUniqueYValues));
  for(unsigned int i = 0; i < nUniqueXValues; ++i) {
    for(unsigned int j = 0; j < nUniqueYValues; ++j) {
      auto index = i * nUniqueYValues + j;
      orderedIds[i][j] = std::get<0>(orderedValues[index]);
    }
  }

  // Create new grid
  vtkNew<vtkPolyData> vtkOutput{};
  vtkOutput->DeepCopy(input);
  auto numberOfCells = (nUniqueXValues - 1) * nUniqueYValues
                       + nUniqueXValues * (nUniqueYValues - 1)
                       + (nUniqueXValues - 1) * (nUniqueYValues - 1);
  vtkOutput->Allocate(numberOfCells);

  for(unsigned int i = 0; i < nUniqueXValues; ++i) {
    for(unsigned int j = 0; j < nUniqueYValues; ++j) {
      if(j != 0) {
        std::array<vtkIdType, 2> linePoints{
          orderedIds[i][j - 1], orderedIds[i][j]};
        vtkOutput->InsertNextCell(VTK_LINE, 2, linePoints.data());
      }

      if(i != 0) {
        std::array<vtkIdType, 2> linePoints{
          orderedIds[i - 1][j], orderedIds[i][j]};
        vtkOutput->InsertNextCell(VTK_LINE, 2, linePoints.data());
      }

      if(i != 0 and j != 0) {
        std::array<vtkIdType, 4> cellPoints{
          orderedIds[i - 1][j - 1], orderedIds[i - 1][j], orderedIds[i][j],
          orderedIds[i][j - 1]};
        vtkOutput->InsertNextCell(VTK_QUAD, 4, cellPoints.data());
      }
    }
  }

  auto noCells = vtkOutput->GetNumberOfCells();
  vtkNew<vtkIntArray> cellTypeArray{};
  cellTypeArray->SetName("CellType");
  cellTypeArray->SetNumberOfTuples(noCells);
  for(int i = 0; i < noCells; ++i) {
    cellTypeArray->SetTuple1(i, vtkOutput->GetCellType(i));
  }
  vtkOutput->GetCellData()->AddArray(cellTypeArray);

  output->ShallowCopy(vtkOutput);

  return 1;
}
