#include <ttkMatrixToHeatMap.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
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
  ttk::Timer tm{};

  const auto input = vtkTable::GetData(inputVector[0]);
  const auto output = vtkUnstructuredGrid::GetData(outputVector);

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

  // generate heat map points
  vtkNew<vtkPoints> points{};
  for(size_t i = 0; i < nInputs + 1; ++i) {
    for(size_t j = 0; j < nRows + 1; ++j) {
      points->InsertNextPoint(i, j, 0.0);
    }
  }

  // generate heat map cells
  vtkNew<vtkCellArray> cells{};
  for(size_t i = 0; i < nInputs; ++i) {
    for(size_t j = 0; j < nRows; ++j) {
      const auto nptline{static_cast<vtkIdType>(nRows + 1)};
      const auto curr{static_cast<vtkIdType>(i * nptline + j)};
      std::array<vtkIdType, 4> ptIds{
        curr, curr + nptline, curr + nptline + 1, curr + 1};
      cells->InsertNextCell(4, ptIds.data());
    }
  }

  // copy distance matrix to heat map cell data
  vtkNew<vtkDoubleArray> dist{};
  vtkNew<vtkDoubleArray> prox{};
  dist->SetNumberOfComponents(1);
  dist->SetName("Distance");
  prox->SetNumberOfComponents(1);
  prox->SetName("Proximity");
  for(size_t i = 0; i < nInputs; ++i) {
    for(size_t j = 0; j < nRows; ++j) {
      const auto val = input->GetColumnByName(ScalarFields[i].data())
                         ->GetVariantValue(j)
                         .ToDouble();
      dist->InsertNextValue(val);
      prox->InsertNextValue(std::exp(-val));
    }
  }

  output->SetPoints(points);
  output->SetCells(VTK_QUAD, cells);
  output->GetCellData()->AddArray(dist);
  output->GetCellData()->AddArray(prox);

  this->printMsg("Complete (#inputs: " + std::to_string(nInputs) + ")", 1.0,
                 tm.getElapsedTime(), 1);

  return 1;
}
