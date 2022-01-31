#include <ttkMetricDistortion.h>

#include <vtkInformation.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMetricDistortion);

/**
 * Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkMetricDistortion::ttkMetricDistortion() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkMetricDistortion::~ttkMetricDistortion() {
}

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMetricDistortion::FillInputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  } else
    return 0;
  return 1;
}

/**
 * Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkMetricDistortion::FillOutputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

/**
 * Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
int ttkMetricDistortion::RequestData(vtkInformation *ttkNotUsed(request),
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {

  // --------------------------------------------------------------------------
  // Get input object from input vector
  // --------------------------------------------------------------------------
  // Load Surface
  vtkPolyData *inputSurface = vtkPolyData::GetData(inputVector[0]);
  if(!inputSurface)
    return 0;

  auto triangulation = ttkAlgorithm::GetTriangulation(inputSurface);
  if(!triangulation)
    return -1;
  this->preconditionTriangulation(triangulation);

  // Load Table
  vtkTable *distanceMatrixVTK = vtkTable::GetData(inputVector[1]);
  std::vector<std::vector<double>> distanceMatrix;
  if(distanceMatrixVTK) {
    auto noRows = distanceMatrixVTK->GetNumberOfRows();
    auto noCols = distanceMatrixVTK->GetNumberOfColumns();
    distanceMatrix = std::vector<std::vector<double>>(
      noRows, std::vector<double>(noCols, 0.0));
    for(unsigned int i = 0; i < noRows; ++i)
      for(unsigned int j = 0; j < noCols; ++j)
        distanceMatrix[i][j]
          = distanceMatrixVTK->GetColumn(j)->GetVariantValue(i).ToDouble();
  }

  // --------------------------------------------------------------------------
  // Call base
  // --------------------------------------------------------------------------
  std::vector<double> surfaceArea, metricArea, ratioArea;
  computeSurfaceArea(triangulation->getData(), distanceMatrix, surfaceArea,
                     metricArea, ratioArea);

  std::vector<double> surfaceDistance, metricDistance, ratioDistance;
  computeSurfaceDistance(triangulation->getData(), distanceMatrix,
                         surfaceDistance, metricDistance, ratioDistance);

  std::vector<double> surfaceCurvature, metricCurvature, diffCurvature;
  computeSurfaceCurvature(triangulation->getData(), distanceMatrix,
                          surfaceCurvature, metricCurvature, diffCurvature);

  // --------------------------------------------------------------------------
  // Get output object
  // --------------------------------------------------------------------------
  auto outputSurface = vtkPolyData::GetData(outputVector, 0);
  outputSurface->DeepCopy(inputSurface);
  
  auto noPoints = inputSurface->GetNumberOfPoints();
  auto noCells = inputSurface->GetNumberOfCells();

  vtkNew<vtkDoubleArray> surfaceCurvatureArray{};
  surfaceCurvatureArray->SetName("CurvatureSurface");
  surfaceCurvatureArray->SetNumberOfTuples(noPoints);
  vtkNew<vtkDoubleArray> metricCurvatureArray{};
  metricCurvatureArray->SetName("CurvatureMetric");
  metricCurvatureArray->SetNumberOfTuples(noPoints);
  vtkNew<vtkDoubleArray> ratioCurvatureArray{};
  ratioCurvatureArray->SetName("CurvatureDiff");
  ratioCurvatureArray->SetNumberOfTuples(noPoints);

  for(unsigned int i = 0; i < noPoints; ++i) {
    surfaceCurvatureArray->SetTuple1(i, surfaceCurvature[i]);
    metricCurvatureArray->SetTuple1(i, metricCurvature[i]);
    ratioCurvatureArray->SetTuple1(i, diffCurvature[i]);
  }

  outputSurface->GetPointData()->AddArray(surfaceCurvatureArray);
  if(distanceMatrix.size() != 0) {
    outputSurface->GetPointData()->AddArray(metricCurvatureArray);
    outputSurface->GetPointData()->AddArray(ratioCurvatureArray);
  }

  vtkNew<vtkDoubleArray> surfaceAreaArray{};
  surfaceAreaArray->SetName("AreaSurface");
  surfaceAreaArray->SetNumberOfTuples(noCells);
  vtkNew<vtkDoubleArray> metricAreaArray{};
  metricAreaArray->SetName("AreaMetric");
  metricAreaArray->SetNumberOfTuples(noCells);
  vtkNew<vtkDoubleArray> ratioAreaArray{};
  ratioAreaArray->SetName("AreaRatio");
  ratioAreaArray->SetNumberOfTuples(noCells);

  vtkNew<vtkDoubleArray> surfaceDistanceArray{};
  surfaceDistanceArray->SetName("DistanceSurface");
  surfaceDistanceArray->SetNumberOfTuples(noCells);
  vtkNew<vtkDoubleArray> metricDistanceArray{};
  metricDistanceArray->SetName("DistanceMetric");
  metricDistanceArray->SetNumberOfTuples(noCells);
  vtkNew<vtkDoubleArray> ratioDistanceArray{};
  ratioDistanceArray->SetName("DistanceRatio");
  ratioDistanceArray->SetNumberOfTuples(noCells);

  for(unsigned int i = 0; i < noCells; ++i) {
    surfaceAreaArray->SetTuple1(i, surfaceArea[i]);
    metricAreaArray->SetTuple1(i, metricArea[i]);
    ratioAreaArray->SetTuple1(i, ratioArea[i]);
    surfaceDistanceArray->SetTuple1(i, surfaceDistance[i]);
    metricDistanceArray->SetTuple1(i, metricDistance[i]);
    ratioDistanceArray->SetTuple1(i, ratioDistance[i]);
  }

  outputSurface->GetCellData()->AddArray(surfaceAreaArray);
  outputSurface->GetCellData()->AddArray(surfaceDistanceArray);
  if(distanceMatrix.size() != 0) {
    outputSurface->GetCellData()->AddArray(metricAreaArray);
    outputSurface->GetCellData()->AddArray(ratioAreaArray);
    outputSurface->GetCellData()->AddArray(metricDistanceArray);
    outputSurface->GetCellData()->AddArray(ratioDistanceArray);
  }

  // return success
  return 1;
}
