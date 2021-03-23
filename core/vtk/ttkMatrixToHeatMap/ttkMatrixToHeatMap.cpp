#include <ttkMatrixToHeatMap.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkTable.h>

#include <array>
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
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkMatrixToHeatMap::RequestData(vtkInformation * /*request*/,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  ttk::Timer tm{};

  const auto input = vtkTable::GetData(inputVector[0]);
  const auto output = vtkPolyData::GetData(outputVector);

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

  const auto nInputs = ScalarFields.size();
  const auto nRows = static_cast<size_t>(input->GetNumberOfRows());
  if(nInputs != nRows) {
    this->printWrn("Distance matrix is not square");
  }

  // grid points
  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints((nInputs + 1) * (nRows + 1));
  // grid cells connectivity arrays
  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(nInputs * nRows + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(4 * nInputs * nRows);
  // distance and proximity arrays
  vtkNew<vtkDoubleArray> dist{}, prox{};
  dist->SetNumberOfComponents(1);
  dist->SetName("Distance");
  dist->SetNumberOfTuples(nInputs * nRows);
  prox->SetNumberOfComponents(1);
  prox->SetName("Proximity");
  prox->SetNumberOfTuples(nInputs * nRows);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs + 1; ++i) {
    for(size_t j = 0; j < nRows + 1; ++j) {
      points->SetPoint(i * (nRows + 1) + j, i, j, 0.0);
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs; ++i) {
    for(size_t j = 0; j < nRows; ++j) {
      // build the cell connectivity and offsets array
      const auto nptline{static_cast<vtkIdType>(nRows + 1)};
      const auto curr{static_cast<vtkIdType>(i * nptline + j)};
      const auto o = i * nRows + j;
      connectivity->SetTuple1(4 * o + 0, curr);
      connectivity->SetTuple1(4 * o + 1, curr + nptline);
      connectivity->SetTuple1(4 * o + 2, curr + nptline + 1);
      connectivity->SetTuple1(4 * o + 3, curr + 1);
      offsets->SetTuple1(o, 4 * o);

      // copy distance matrix to heat map cell data
      const auto val = input->GetColumnByName(ScalarFields[i].data())
                         ->GetVariantValue(j)
                         .ToDouble();
      dist->SetTuple1(i * nRows + j, val);
      prox->SetTuple1(i * nRows + j, std::exp(-val));
    }
  }
  offsets->SetTuple1(nInputs * nRows, connectivity->GetNumberOfTuples());

  // gather arrays to make the PolyData
  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  output->SetPoints(points);
  output->SetPolys(cells);
  output->GetCellData()->AddArray(dist);
  output->GetCellData()->AddArray(prox);

  this->printMsg("Complete (#inputs: " + std::to_string(nInputs) + ")", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);

  return 1;
}
