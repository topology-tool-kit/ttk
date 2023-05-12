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

vtkStandardNewMacro(ttkMarchingTetrahedra);

ttkMarchingTetrahedra::ttkMarchingTetrahedra() {
  this->setDebugMsgPrefix("MarchingTetrahedra");
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
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

template <typename vtkArrayType, typename vectorType>
void setArray(vtkArrayType &vtkArray, vectorType &vector) {
  ttkUtils::SetVoidArray(vtkArray, vector.data(), vector.size(), 1);
}

template <typename scalarType, typename triangulationType>
int ttkMarchingTetrahedra::dispatch(vtkDataArray *const inputScalars,
                                    vtkPolyData *const outputSeparators,
                                    const triangulationType &triangulation) {

  const auto scalars = ttkUtils::GetPointer<scalarType>(inputScalars);

  SimplexId numberOfPoints{};
  SimplexId numberOfCells{};
  output_points.clear();
  output_cells_mscIds.clear();

  this->setOutput(&numberOfPoints, &output_points, &numberOfCells,
                  &output_cells_connectivity, &output_cells_mscIds);

  const int ret
    = this->execute<scalarType, triangulationType>(scalars, triangulation);

  if(ret != 0)
    return !this->printErr("MarchingTetrahedra.execute() error");

  vtkNew<vtkFloatArray> pointsCoords{};
  pointsCoords->SetNumberOfComponents(3);
  setArray(pointsCoords, output_points);

  vtkNew<ttkSimplexIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(numberOfCells + 1);
  connectivity->SetNumberOfComponents(1);
  setArray(connectivity, output_cells_connectivity);

  vtkNew<vtkUnsignedLongLongArray> mscIds{};
  mscIds->SetNumberOfComponents(1);
  mscIds->SetName("MSCIds");
  setArray(mscIds, output_cells_mscIds);

  if(triangulation.getDimensionality() == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < numberOfCells + 1; ++i) {
      offsets->SetTuple1(i, 3 * i);
    }
  } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < numberOfCells + 1; ++i) {
      offsets->SetTuple1(i, 2 * i);
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
  if(triangulation.getDimensionality() == 3) {
    outputSeparators->SetPolys(cells);
  } else {
    outputSeparators->SetLines(cells);
  }

  auto cellData = outputSeparators->GetCellData();
  cellData->AddArray(mscIds);

  return ret;
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

  int ret{};

  ttkVtkTemplateMacro(inputScalars->GetDataType(), triangulation->getType(),
                      (ret = dispatch<VTK_TT, TTK_TT>(
                         inputScalars, outputSeparators,
                         *static_cast<TTK_TT *>(triangulation->getData()))));

  if(ret != 0) {
    return -1;
  }

  return !ret;
}
