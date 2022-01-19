#include <ttkRipsComplex.h>
#include <ttkUtils.h>

#include <vtkAbstractArray.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkTable.h>
#include <vtkType.h>
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
  ttk::Timer tm{};

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

  const SimplexId numberOfRows = input->GetNumberOfRows();
  const SimplexId numberOfColumns = ScalarFields.size();

  if(numberOfRows <= 0 || numberOfColumns <= 0) {
    this->printErr("Input matrix has invalid dimensions (rows: "
                   + std::to_string(numberOfRows)
                   + ", columns: " + std::to_string(numberOfColumns) + ")");
    return 0;
  }

  std::array<vtkAbstractArray *, 3> components{
    input->GetColumnByName(this->XColumn.data()),
    input->GetColumnByName(this->YColumn.data()),
    input->GetColumnByName(this->ZColumn.data()),
  };

  // points from components arrays
  vtkNew<vtkPoints> points{};
  const auto nPoints = components[0]->GetNumberOfTuples();
  points->SetNumberOfPoints(nPoints);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(vtkIdType i = 0; i < nPoints; ++i) {
    std::array<float, 3> coords{};
    if(components[0] != nullptr) {
      coords[0] = components[0]->GetVariantValue(i).ToFloat();
    }
    if(components[1] != nullptr && components[1] != components[0]) {
      coords[1] = components[1]->GetVariantValue(i).ToFloat();
    }
    if(components[2] != nullptr && components[2] != components[0]
       && components[2] != components[1]) {
      coords[2] = components[2]->GetVariantValue(i).ToFloat();
    }
    points->SetPoint(i, coords.data());
  }

  std::vector<vtkAbstractArray *> arrays{};
  for(const auto &s : ScalarFields) {
    arrays.push_back(input->GetColumnByName(s.data()));
  }

  std::vector<std::vector<double>> inputMatrix(arrays.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < numberOfRows; ++i) {
    for(size_t j = 0; j < arrays.size(); ++j) {
      inputMatrix[i].emplace_back(arrays[j]->GetVariantValue(i).ToDouble());
    }
  }

  std::vector<ttk::LongSimplexId> vec_connectivity{};

  this->execute(vec_connectivity, inputMatrix);

  const auto nCells = vec_connectivity.size() / (this->OutputDimension + 1);
  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(nCells + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(vec_connectivity.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nCells + 1; ++i) {
    offsets->SetTuple1(i, i * (this->OutputDimension + 1));
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < vec_connectivity.size(); ++i) {
    connectivity->SetTuple1(i, vec_connectivity[i]);
  }

  // gather arrays to make the PolyData
  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);

  const auto getCellType = [](const SimplexId dim) -> int {
    if(dim == 3) {
      return VTK_TETRA;
    } else if(dim == 2) {
      return VTK_TRIANGLE;
    } else if(dim == 1) {
      return VTK_LINE;
    }
    return VTK_VERTEX;
  };

  output->SetPoints(points);
  output->SetCells(getCellType(this->OutputDimension), cells);

  if(this->KeepAllDataArrays) {
    auto pd = output->GetPointData();
    if(pd != nullptr) {
      for(vtkIdType i = 0; i < input->GetNumberOfColumns(); ++i) {
        pd->AddArray(input->GetColumn(i));
      }
    }
  }

  this->printMsg("Produced " + std::to_string(nCells) + " cells of dimension "
                   + std::to_string(this->OutputDimension),
                 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 1;
}
