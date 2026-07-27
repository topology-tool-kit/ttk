#include <ttkMetricDistortion.h>

#include <vtkInformation.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
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

ttkMetricDistortion::~ttkMetricDistortion() = default;

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
template <class tableDataType>
int ttkMetricDistortion::run(vtkInformationVector **inputVector) {
  vtkPointSet *inputSurface = vtkPointSet::GetData(inputVector[0]);
  auto triangulation = ttkAlgorithm::GetTriangulation(inputSurface);
  if(!triangulation) {
    printErr("Unable to load triangulation.");
    return -4;
  }
  this->preconditionTriangulation(triangulation);

  // Load Table
  vtkTable *distanceMatrixVTK = vtkTable::GetData(inputVector[1]);
  std::vector<tableDataType *> distanceMatrix;
  if(distanceMatrixVTK) {
    auto noRows = distanceMatrixVTK->GetNumberOfRows();
    distanceMatrix = std::vector<tableDataType *>(noRows);
    for(unsigned int i = 0; i < noRows; ++i) {
      auto row = vtkDataArray::SafeDownCast(distanceMatrixVTK->GetColumn(i));
      if(!row) {
        printErr("Unable to load column " + std::to_string(i)
                 + " in the distance matrix.");
        return -5;
      }
      distanceMatrix[i] = ttkUtils::GetPointer<tableDataType>(row);
    }
  }

  // --------------------------------------------------------------------------
  // Call base
  // --------------------------------------------------------------------------
  computeSurfaceArea(triangulation->getData(), distanceMatrix, surfaceArea_,
                     metricArea_, ratioArea_);

  computeSurfaceDistance(triangulation->getData(), distanceMatrix,
                         surfaceDistance_, metricDistance_, ratioDistance_,
                         surfacePointDistance_, metricPointDistance_,
                         ratioPointDistance_);

  computeSurfaceCurvature(triangulation->getData(), distanceMatrix,
                          surfaceCurvature_, metricCurvature_, diffCurvature_);

  return 1;
}

int ttkMetricDistortion::RequestData(vtkInformation *ttkNotUsed(request),
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {

  // --------------------------------------------------------------------------
  // Get input object from input vector
  // --------------------------------------------------------------------------
  // Load Surface
  vtkPointSet *inputSurface = vtkPointSet::GetData(inputVector[0]);
  if(!inputSurface) {
    printErr("Unable to load input surface.");
    return 0;
  }
  auto noPoints = inputSurface->GetNumberOfPoints();
  if(noPoints == 0) {
    printErr("Input surface should not have 0 points.");
    return -1;
  }
  auto noCells = inputSurface->GetNumberOfCells();

  vtkTable *distanceMatrixVTK = vtkTable::GetData(inputVector[1]);
  bool const validDistanceMatrix
    = (distanceMatrixVTK
       and distanceMatrixVTK->GetNumberOfColumns() == noPoints);
  if(distanceMatrixVTK and not validDistanceMatrix) {
    printErr("Distance matrix should have the same number of rows/columns than "
             "the surface's number of points.");
    return -2;
  }
  if(distanceMatrixVTK
     and distanceMatrixVTK->GetNumberOfColumns()
           != distanceMatrixVTK->GetNumberOfRows()) {
    printErr("Distance matrix should be square.");
    return -3;
  }
  int const tableDataType
    = (validDistanceMatrix ? distanceMatrixVTK->GetColumn(0)->GetDataType()
                           : VTK_DOUBLE);

  surfaceArea_.clear();
  metricArea_.clear();
  ratioArea_.clear();
  surfaceDistance_.clear();
  metricDistance_.clear();
  ratioDistance_.clear();
  surfacePointDistance_.clear();
  metricPointDistance_.clear();
  ratioPointDistance_.clear();
  surfaceCurvature_.clear();
  metricCurvature_.clear();
  diffCurvature_.clear();
  int res = 1;
  switch(tableDataType) {
    vtkTemplateMacro(res = this->run<VTK_TT>(inputVector));
  }
  if(res != 1)
    return res;

  // --------------------------------------------------------------------------
  // Get output object
  // --------------------------------------------------------------------------
  // --- Point Data
  auto outputSurface = vtkPointSet::GetData(outputVector, 0);
  outputSurface->DeepCopy(inputSurface);

  vtkNew<vtkDoubleArray> surfaceCurvature_Array{};
  surfaceCurvature_Array->SetName("SurfaceCurvature");
  surfaceCurvature_Array->SetNumberOfTuples(noPoints);
  vtkNew<vtkDoubleArray> metricCurvature_Array{};
  metricCurvature_Array->SetName("MetricCurvature");
  metricCurvature_Array->SetNumberOfTuples(noPoints);
  vtkNew<vtkDoubleArray> ratioCurvatureArray{};
  ratioCurvatureArray->SetName("CurvatureDiff");
  ratioCurvatureArray->SetNumberOfTuples(noPoints);

  for(unsigned int i = 0; i < noPoints; ++i) {
    surfaceCurvature_Array->SetTuple1(i, surfaceCurvature_[i]);
    metricCurvature_Array->SetTuple1(i, metricCurvature_[i]);
    ratioCurvatureArray->SetTuple1(i, diffCurvature_[i]);
  }

  outputSurface->GetPointData()->AddArray(surfaceCurvature_Array);
  if(validDistanceMatrix) {
    outputSurface->GetPointData()->AddArray(metricCurvature_Array);
    outputSurface->GetPointData()->AddArray(ratioCurvatureArray);
  }

  for(unsigned int i = 0; i < 3; ++i) {
    std::string const type{(i == 0 ? "Min" : (i == 1) ? "Max" : "Avg")};
    vtkNew<vtkDoubleArray> surfaceIDistanceArray{};
    surfaceIDistanceArray->SetName((type + "SurfaceEdgeLength").c_str());
    surfaceIDistanceArray->SetNumberOfTuples(noPoints);
    vtkNew<vtkDoubleArray> metricIDistanceArray{};
    metricIDistanceArray->SetName((type + "MetricEdgeLength").c_str());
    metricIDistanceArray->SetNumberOfTuples(noPoints);
    vtkNew<vtkDoubleArray> ratioIDistanceArray{};
    ratioIDistanceArray->SetName((type + "EdgeLengthRatio").c_str());
    ratioIDistanceArray->SetNumberOfTuples(noPoints);
    for(unsigned int j = 0; j < noPoints; ++j) {
      surfaceIDistanceArray->SetTuple1(j, surfacePointDistance_[j][i]);
      metricIDistanceArray->SetTuple1(j, metricPointDistance_[j][i]);
      ratioIDistanceArray->SetTuple1(j, ratioPointDistance_[j][i]);
    }
    outputSurface->GetPointData()->AddArray(surfaceIDistanceArray);
    if(validDistanceMatrix) {
      outputSurface->GetPointData()->AddArray(metricIDistanceArray);
      outputSurface->GetPointData()->AddArray(ratioIDistanceArray);
    }
  }

  // --- Cell Data
  vtkNew<vtkDoubleArray> surfaceArea_Array{};
  surfaceArea_Array->SetName("SurfaceArea");
  surfaceArea_Array->SetNumberOfTuples(noCells);
  vtkNew<vtkDoubleArray> metricArea_Array{};
  metricArea_Array->SetName("MetricArea");
  metricArea_Array->SetNumberOfTuples(noCells);
  vtkNew<vtkDoubleArray> ratioArea_Array{};
  ratioArea_Array->SetName("AreaRatio");
  ratioArea_Array->SetNumberOfTuples(noCells);

  vtkNew<vtkDoubleArray> surfaceDistance_Array{};
  surfaceDistance_Array->SetName("SurfaceEdgeLength");
  surfaceDistance_Array->SetNumberOfTuples(noCells);
  vtkNew<vtkDoubleArray> metricDistance_Array{};
  metricDistance_Array->SetName("MetricEdgeLength");
  metricDistance_Array->SetNumberOfTuples(noCells);
  vtkNew<vtkDoubleArray> ratioDistance_Array{};
  ratioDistance_Array->SetName("EdgeLengthRatio");
  ratioDistance_Array->SetNumberOfTuples(noCells);

  bool distanceAllNan = true;
  for(unsigned int i = 0; i < noCells; ++i) {
    surfaceArea_Array->SetTuple1(i, surfaceArea_[i]);
    metricArea_Array->SetTuple1(i, metricArea_[i]);
    ratioArea_Array->SetTuple1(i, ratioArea_[i]);
    surfaceDistance_Array->SetTuple1(i, surfaceDistance_[i]);
    metricDistance_Array->SetTuple1(i, metricDistance_[i]);
    ratioDistance_Array->SetTuple1(i, ratioDistance_[i]);
    if(surfaceDistance_[i] == surfaceDistance_[i]) // detect NaN
      distanceAllNan = false;
  }

  outputSurface->GetCellData()->AddArray(surfaceArea_Array);
  if(validDistanceMatrix) {
    outputSurface->GetCellData()->AddArray(metricArea_Array);
    outputSurface->GetCellData()->AddArray(ratioArea_Array);
  }
  if(not distanceAllNan) {
    outputSurface->GetCellData()->AddArray(surfaceDistance_Array);
    if(validDistanceMatrix) {
      outputSurface->GetCellData()->AddArray(metricDistance_Array);
      outputSurface->GetCellData()->AddArray(ratioDistance_Array);
    }
  }

  // return success
  return 1;
}
