#include <ttkProjectionFromTable.h>

#include <vtkInformation.h>

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <set>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkProjectionFromTable);

/**
 * Implement the filter constructor and destructor in the cpp file.
 */
ttkProjectionFromTable::ttkProjectionFromTable() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

/**
 * Specify the required input data type of each input port
 */
int ttkProjectionFromTable::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
  } else
    return 0;
  return 1;
}

/**
 * Specify the data object type of each output port
 */
int ttkProjectionFromTable::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  } else
    return 0;
  return 1;
}

/**
 * Pass VTK data to the base code and convert base code output to VTK
 */
template <typename VTK_T1, typename VTK_T2>
void ttkProjectionFromTable::dispatch(
  std::vector<std::tuple<int, double, double>> &storage,
  const VTK_T1 *const values,
  const VTK_T2 *const values2,
  const size_t nvalues) {

  for(size_t i = 0; i < nvalues; ++i) {
    storage.emplace_back(
      i, static_cast<double>(values[i]), static_cast<double>(values2[i]));
  }
}

int ttkProjectionFromTable::RequestData(vtkInformation *ttkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
  // --------------------------------------------------------------------------
  // --- Get input object from input vector
  // --------------------------------------------------------------------------
  // Get surface
  auto surface = vtkPolyData::GetData(inputVector[0], 0);
  auto triangulation = ttkAlgorithm::GetTriangulation(surface);
  if(!triangulation) {
    printErr("Unable to load triangulation.");
    return 0;
  }

  // Get surface arrays
  const auto surfaceXArray = this->GetInputArrayToProcess(0, inputVector);
  const auto surfaceYArray = this->GetInputArrayToProcess(1, inputVector);
  if(surfaceXArray == nullptr) {
    this->printErr("Cannot find the required surface data X array");
    return 0;
  }
  if(surfaceYArray == nullptr) {
    this->printErr("Cannot find the required surface data Y array");
    return 0;
  }

  // Get table arrays
  auto coefficients = vtkTable::GetData(inputVector[1], 0);
  auto xArrayName = surfaceXArray->GetName();
  auto tableXArray = coefficients->GetColumnByName(xArrayName);
  if(!tableXArray) {
    printErr("Can not find " + std::string{xArrayName} + " X array in table");
    return 0;
  }
  auto tableDataXArray = vtkDataArray::SafeDownCast(tableXArray);
  if(!tableDataXArray) {
    printErr("Can not load " + std::string{xArrayName}
             + " X data array in table");
    return 0;
  }
  auto yArrayName = surfaceYArray->GetName();
  auto tableYArray = coefficients->GetColumnByName(yArrayName);
  if(!tableYArray) {
    printErr("Can not find " + std::string{yArrayName} + " Y array in table");
    return 0;
  }
  auto tableDataYArray = vtkDataArray::SafeDownCast(tableYArray);
  if(!tableDataYArray) {
    printErr("Can not load " + std::string{yArrayName}
             + " Y data array in table");
    return 0;
  }

  const auto nSurfaceValues = surfaceXArray->GetNumberOfTuples();
  const auto nTableValues = tableDataXArray->GetNumberOfTuples();

  // --------------------------------------------------------------------------
  // --- Get ordered values
  // --------------------------------------------------------------------------
  // store point index <-> ordering value in vector
  std::vector<std::tuple<int, double, double>> orderedValues;

#ifndef TTK_ENABLE_DOUBLE_TEMPLATING
  switch(surfaceXArray->GetDataType()) {
    vtkTemplateMacro(
      dispatch(orderedValues,
               static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(surfaceXArray)),
               static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(surfaceYArray)),
               nSurfaceValues));
  }
#else
  switch(vtkTemplate2PackMacro(
    surfaceXArray->GetDataType(), surfaceYArray->GetDataType())) {
    vtkTemplate2Macro(
      dispatch(orderedValues,
               static_cast<VTK_T1 *>(ttkUtils::GetVoidPointer(surfaceXArray)),
               static_cast<VTK_T2 *>(ttkUtils::GetVoidPointer(surfaceYArray)),
               nSurfaceValues));
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
  const long nUniqueXValues
    = std::unique(xValues.begin(), xValues.end()) - xValues.begin();
  TTK_PSORT(this->threadNumber_, yValues.begin(), yValues.end());
  const long nUniqueYValues
    = std::unique(yValues.begin(), yValues.end()) - yValues.begin();
  std::array<const long, 2> surfaceDim{nUniqueXValues, nUniqueYValues};

  if(nUniqueXValues * nUniqueYValues != surface->GetNumberOfPoints()) {
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

  //---------------------------------------------------------------------------
  // --- Call base
  //---------------------------------------------------------------------------
  std::vector<std::vector<double>> inputPoints;
#ifndef TTK_ENABLE_DOUBLE_TEMPLATING
  switch(tableDataXArray->GetDataType()) {
    vtkTemplateMacro(computeInputPoints(
      triangulation, orderedValues, surfaceDim,
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(tableDataXArray)),
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(tableDataYArray)),
      nTableValues, inputPoints));
  }
#else
  switch(vtkTemplate2PackMacro(
    tableDataXArray->GetDataType(), tableDataYArray->GetDataType())) {
    vtkTemplate2Macro(computeInputPoints(
      triangulation, orderedValues, surfaceDim,
      static_cast<VTK_T1 *>(ttkUtils::GetVoidPointer(tableDataXArray)),
      static_cast<VTK_T2 *>(ttkUtils::GetVoidPointer(tableDataYArray)),
      nTableValues, inputPoints));
  }
#endif // TTK_ENABLE_DOUBLE_TEMPLATING

  //---------------------------------------------------------------------------
  // --- Create output
  //---------------------------------------------------------------------------
  auto output = vtkPolyData::GetData(outputVector, 0);
  output->DeepCopy(surface);

  std::set<std::string> toRemove;
  for(unsigned int i = 0; i < inputPoints.size(); ++i) {
    double point[3];
    for(unsigned int j = 0; j < 3; ++j) {
      point[j] = inputPoints[i][j];
    }
    output->GetPoints()->InsertNextPoint(point);

    // Merge data of arrays in surface also in table or add empty data
    for(int j = 0; j < output->GetPointData()->GetNumberOfArrays(); ++j) {
      auto outputArray = output->GetPointData()->GetAbstractArray(j);
      std::string const name = outputArray->GetName();
      auto array = coefficients->GetColumnByName(name.c_str());
      if(array) {
        outputArray->InsertNextTuple(i, array);
      } else {
        auto dataArray = vtkDataArray::SafeDownCast(outputArray);
        if(dataArray) {
          const double val = std::nan("");
          dataArray->InsertNextTuple(&val);
        } else {
          auto stringArray = vtkStringArray::SafeDownCast(outputArray);
          if(stringArray)
            stringArray->InsertNextValue("");
          else
            toRemove.insert(name);
        }
      }
    }
  }

  int const noPointsOri = output->GetNumberOfPoints() - inputPoints.size();

  // Merge data of arrays in table also in surface or add empty data
  std::set<std::string> toGet;
  for(int j = 0; j < coefficients->GetNumberOfColumns(); ++j) {
    auto inputArray = coefficients->GetColumn(j);
    auto dataArray = vtkDataArray::SafeDownCast(inputArray);
    auto stringArray = vtkStringArray::SafeDownCast(inputArray);
    std::string const name = inputArray->GetName();
    auto array = output->GetPointData()->GetAbstractArray(name.c_str());
    if(not array and (dataArray or stringArray))
      toGet.insert(name);
  }
  vtkNew<vtkTable> tableCopy{};
  if(toGet.size() != 0)
    tableCopy->DeepCopy(coefficients);
  for(auto &name : toGet) {
    auto inputArray = tableCopy->GetColumnByName(name.c_str());
    auto dataArray = vtkDataArray::SafeDownCast(inputArray);
    auto stringArray = vtkStringArray::SafeDownCast(inputArray);
    inputArray->SetNumberOfTuples(output->GetNumberOfPoints());
    for(int i = 0; i < output->GetNumberOfPoints(); ++i) {
      // if is surface
      if(i < noPointsOri) {
        if(dataArray) {
          double const val = std::nan("");
          dataArray->SetTuple1(i, val);
        } else if(stringArray) {
          stringArray->SetValue(i, "");
        } else
          printWrn("Can not convert " + name + " to dataArray or stringArray.");
      } else {
        inputArray->SetTuple(
          i, i - noPointsOri, coefficients->GetColumnByName(name.c_str()));
      }
    }
    output->GetPointData()->AddArray(inputArray);
  }

  for(auto &name : toRemove)
    output->GetPointData()->RemoveArray(name.c_str());

  // Some points data
  vtkNew<vtkIntArray> isSurfaceArray{};
  isSurfaceArray->SetName("isSurface");
  isSurfaceArray->SetNumberOfTuples(output->GetNumberOfPoints());
  for(int i = 0; i < output->GetNumberOfPoints(); ++i) {
    bool const isSurface = (i < noPointsOri);
    isSurfaceArray->SetTuple1(i, isSurface);
  }
  output->GetPointData()->AddArray(isSurfaceArray);

  // Some cells data
  auto noCells = output->GetNumberOfCells();
  vtkNew<vtkIntArray> cellTypeArray{};
  cellTypeArray->SetName("CellType");
  cellTypeArray->SetNumberOfTuples(noCells);
  for(int i = 0; i < noCells; ++i) {
    cellTypeArray->SetTuple1(i, output->GetCellType(i));
  }
  output->GetCellData()->AddArray(cellTypeArray);

  return 1;
}
