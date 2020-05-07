#include "ttkBottleneckDistance.h"

vtkStandardNewMacro(ttkBottleneckDistance);

ttkBottleneckDistance::ttkBottleneckDistance() {
  // settings
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(3);

  // inputs
  this->BenchmarkSize = -1;
  this->UseOutputMatching = false;
  this->UsePersistenceMetric = false;
  this->UseGeometricSpacing = false;
  this->WassersteinMetric = "2";
  this->Alpha = 1.0;
  this->Tolerance = 1.0;
  this->PX = 0;
  this->PY = 0;
  this->PZ = 0;
  this->PE = 1;
  this->PS = 1;
  this->Spacing = 0.0;
  this->Is3D = false;
  this->PVAlgorithm = -1;

  // outputs
  this->result = -1.;
  this->CTPersistenceDiagram1_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
  this->CTPersistenceDiagram2_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
  this->CTPersistenceDiagram3_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
}
ttkBottleneckDistance::~ttkBottleneckDistance(){};

int ttkBottleneckDistance::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0 || port == 1)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  else
    return 0;

  return 1;
}

int ttkBottleneckDistance::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;

  return 1;
}

int ttkBottleneckDistance::doBenchmark() {
  using dataType = double;

  std::vector<diagramTuple> CTDiagram1;
  std::vector<diagramTuple> CTDiagram2;

  int benchmarkSize = BenchmarkSize;
  int status = 0;
  status = generatePersistenceDiagram<double>(CTDiagram1, benchmarkSize);
  if(status < 0)
    return status;
  status = generatePersistenceDiagram<double>(CTDiagram2, 4 * benchmarkSize);
  if(status < 0)
    return status;

  this->setPersistencePercentThreshold(Tolerance);
  this->setPX(PX);
  this->setPY(PY);
  this->setPZ(PZ);
  this->setPE(PE);
  this->setPS(PS);
  this->setCTDiagram1(&CTDiagram1);
  this->setCTDiagram2(&CTDiagram2);

  std::string wassersteinMetric = WassersteinMetric;
  this->setWasserstein(wassersteinMetric);
  std::string algorithm = DistanceAlgorithm;
  this->setAlgorithm(algorithm);
  int pvAlgorithm = PVAlgorithm;
  this->setPVAlgorithm(pvAlgorithm);

  // Empty matchings.
  auto matchings = new std::vector<diagramTuple>();
  this->setOutputMatchings(matchings);

  // Exec.
  bool usePersistenceMetric = UsePersistenceMetric;
  // double alpha = Alpha;
  status = this->execute<dataType>(usePersistenceMetric);

  if(status != 0) {
    return status;
  }

  return 0;
}

int ttkBottleneckDistance::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  using dataType = double;

  int benchmarkSize = BenchmarkSize;
  bool benchmark = benchmarkSize > 0;
  if(benchmark) {
    return doBenchmark();
  }

  // Prepare IO
  auto CTPersistenceDiagram1_ = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto CTPersistenceDiagram2_ = vtkUnstructuredGrid::GetData(inputVector[1]);
  auto outputCT1 = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputCT2 = vtkUnstructuredGrid::GetData(outputVector, 1);
  auto outputCT3 = vtkUnstructuredGrid::GetData(outputVector, 2);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!CTPersistenceDiagram1_ || !CTPersistenceDiagram2_) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }

  if(!outputCT1 || !outputCT2 || !outputCT3) {
    this->printErr("Output pointer is NULL.");
    return -1;
  }

  if(CTPersistenceDiagram1_->GetNumberOfPoints() == 0
     || CTPersistenceDiagram2_->GetNumberOfPoints() == 0) {
    this->printErr("Input has no points.");
    return -1;
  }
#endif

  this->setPersistencePercentThreshold(Tolerance);
  this->setPX(PX);
  this->setPY(PY);
  this->setPZ(PZ);
  this->setPE(PE);
  this->setPS(PS);

  if(!CTPersistenceDiagram1_ || !CTPersistenceDiagram2_ || !outputCT3)
    return -1;

  int dataType1 = CTPersistenceDiagram1_->GetCellData()
                    ->GetArray("Persistence")
                    ->GetDataType();
  int dataType2 = CTPersistenceDiagram2_->GetCellData()
                    ->GetArray("Persistence")
                    ->GetDataType();
  if(dataType1 != dataType2)
    return -1;

  vtkDoubleArray *birthScalars1 = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram1_->GetPointData()->GetArray("Birth"));
  vtkDoubleArray *deathScalars1 = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram1_->GetPointData()->GetArray("Death"));
  vtkDoubleArray *birthScalars2 = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram1_->GetPointData()->GetArray("Birth"));
  vtkDoubleArray *deathScalars2 = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram1_->GetPointData()->GetArray("Death"));
  bool is2D1 = !deathScalars1 && !birthScalars1;
  bool is2D2 = !deathScalars2 && !birthScalars2;
  if(is2D1 != is2D2)
    return -2;
  bool is2D = is2D1;

  // Call package
  int status = 0;

  //  switch (dataType1) {
  //    vtkTemplateMacro(({
  // TODO template my methods
  std::vector<diagramTuple> CTDiagram1;
  std::vector<diagramTuple> CTDiagram2;

  status = getPersistenceDiagram<dataType>(
    CTDiagram1, CTPersistenceDiagram1_, Spacing, 0);
  if(status < 0) {
    return -2;
  }

  status = getPersistenceDiagram<dataType>(
    CTDiagram2, CTPersistenceDiagram2_, Spacing, 1);
  if(status < 0) {
    return -2;
  }

  this->setCTDiagram1(&CTDiagram1);
  this->setCTDiagram2(&CTDiagram2);

  std::string wassersteinMetric = WassersteinMetric;
  this->setWasserstein(wassersteinMetric);
  std::string algorithm = DistanceAlgorithm;
  this->setAlgorithm(algorithm);
  int pvAlgorithm = PVAlgorithm;
  this->setPVAlgorithm(pvAlgorithm);

  // Empty matchings.
  std::vector<matchingTuple> matchings;
  this->setOutputMatchings(&matchings);

  // Exec.
  bool usePersistenceMetric = UsePersistenceMetric;
  // double alpha = Alpha;
  status = this->execute<dataType>(usePersistenceMetric);
  if(status != 0)
    return status;

  // Apply results to outputs 0 and 1.
  status = augmentPersistenceDiagrams<dataType>(
    CTDiagram1, CTDiagram2, matchings, CTPersistenceDiagram1_,
    CTPersistenceDiagram2_);

  bool useOutputMatching = UseOutputMatching;
  bool useGeometricSpacing = UseGeometricSpacing;

  // Apply results to output 2.
  if(useOutputMatching) {
    status = getMatchingMesh<dataType>(
      CTDiagram1, CTDiagram2, matchings, useGeometricSpacing, Spacing, is2D);
  }

  if(status != 0)
    return status;
  //    }));
  //  }

  // Set output.
  outputCT1->ShallowCopy(CTPersistenceDiagram1_);
  outputCT2->DeepCopy(CTPersistenceDiagram2_);
  if(UseGeometricSpacing)
    translateSecondDiagram<dataType>(outputCT2, Spacing);

  if(UseOutputMatching)
    outputCT3->ShallowCopy(CTPersistenceDiagram3_);

  return status;
}
