#include <ttkMacros.h>
#include <ttkMarchingTetrahedra.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSignedCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedLongLongArray.h>

namespace {
  template <typename vtkArrayType, typename vectorType>
  void setArray(vtkArrayType &vtkArray, vectorType &vector) {
    ttkUtils::SetVoidArray(vtkArray, vector.data(), vector.size(), 1);
  }
} // namespace

vtkStandardNewMacro(ttkMarchingTetrahedra);

ttkMarchingTetrahedra::ttkMarchingTetrahedra() {
  this->setDebugMsgPrefix("MarchingTetrahedra");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkMarchingTetrahedra::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkMarchingTetrahedra::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

template <typename scalarType, typename triangulationType>
int ttkMarchingTetrahedra::dispatch(vtkDataArray *const inputScalars,
                                    vtkPolyData *const outputSeparators,
                                    const triangulationType &triangulation) {

  const auto scalars = ttkUtils::GetPointer<scalarType>(inputScalars);
  const int dim = triangulation.getDimensionality();

  output_points_.clear();
  output_cells_labels_.clear();
  output_cells_connectivity_.clear();

  const int status
    = this->execute<scalarType, triangulationType>(scalars, triangulation);

  if(status != 0)
    return !this->printErr("MarchingTetrahedra.execute() error");

  vtkNew<vtkFloatArray> pointsCoords{};
  pointsCoords->SetNumberOfComponents(3);
  setArray(pointsCoords, output_points_);

  vtkNew<ttkSimplexIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(output_numberOfCells_ + 1);
  connectivity->SetNumberOfComponents(1);
  setArray(connectivity, output_cells_connectivity_);

  vtkNew<vtkUnsignedLongLongArray> hashArr{};
  hashArr->SetNumberOfComponents(1);
  hashArr->SetName("Hash");
  setArray(hashArr, output_cells_labels_);

  if(dim == 2 || dim == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < output_numberOfCells_ + 1; ++i) {
      offsets->SetTuple1(i, dim * i);
    }
  }

  vtkNew<vtkPoints> points{};
  points->SetData(pointsCoords);
  outputSeparators->SetPoints(points);

  vtkNew<vtkCellArray> cells{};
#ifndef TTK_ENABLE_64BIT_IDS
  cells->Use32BitStorage();
#endif // TTK_ENABLE_64BIT_IDS
  cells->SetData(offsets, connectivity);
  if(dim == 3) {
    outputSeparators->SetPolys(cells);
  } else {
    outputSeparators->SetLines(cells);
  }

  auto cellData = outputSeparators->GetCellData();
  cellData->AddArray(hashArr);

  return 1;
}

int ttkMarchingTetrahedra::RequestData(vtkInformation *ttkNotUsed(request),
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputSeparators = vtkPolyData::GetData(outputVector, 0);

  if(!input)
    return !this->printErr("Input pointer is NULL.");

  if(input->GetNumberOfPoints() == 0)
    return !this->printErr("Input has no point.");

  if(!outputSeparators)
    return !this->printErr("Output pointers are NULL.");

  const auto triangulation = ttkAlgorithm::GetTriangulation(input);

  if(triangulation == nullptr)
    return !this->printErr("Triangulation is null");

  const auto inputScalars = this->GetInputArrayToProcess(0, inputVector);

  if(inputScalars == nullptr)
    return !this->printErr("wrong scalars.");

  this->printMsg("Launching computation on field `"
                 + std::string(inputScalars->GetName()) + "'...");

  const SimplexId numberOfVertices = triangulation->getNumberOfVertices();

  if(!numberOfVertices)
    return !this->printErr("Input has no vertices.");

  int status{};

  ttkVtkTemplateMacro(inputScalars->GetDataType(), triangulation->getType(),
                      (status = dispatch<VTK_TT, TTK_TT>(
                         inputScalars, outputSeparators,
                         *static_cast<TTK_TT *>(triangulation->getData()))));

  return status;
}
