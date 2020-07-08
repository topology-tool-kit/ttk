#include "ttkBottleneckDistance.h"
#include <ttkMacros.h>
#include <ttkUtils.h>

#include "vtkObjectFactory.h"

vtkStandardNewMacro(ttkBottleneckDistance);

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
  // this->setThreadNumber(thread);

  // Empty matchings.
  auto matchings = new std::vector<diagramTuple>();
  this->setOutputMatchings(matchings);

  // Exec.
  bool usePersistenceMetric = UsePersistenceMetric;
  status = this->execute<dataType>(usePersistenceMetric);

  return status;
}

int ttkBottleneckDistance::RequestData(vtkInformation * /*request*/,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {
  using dataType = double;

  int benchmarkSize = BenchmarkSize;
  bool benchmark = benchmarkSize > 0;
  if(benchmark) {
    return doBenchmark();
  }

  // Prepare IO
  vtkDataSet *input1 = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *input2 = vtkDataSet::GetData(inputVector[1]);

  // vtkDataSet *output1 =
  // outputVector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  // vtkDataSet *output2 =
  // outputVector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));
  // vtkDataSet *output3 =
  // outputVector->GetInformationObject(2)->Get(vtkDataObject::DATA_OBJECT()));

  // #ifndef TTK_ENABLE_KAMIKAZE
  // if(!input1 || !input2) {
  //   cerr << "[ttkBottleneckDistance] Error: input pointer is NULL." << endl;
  //   return -1;
  // }

  // if(!output1 || !output2 || !output3) {
  //   cerr << "[ttkBottleneckDistance] Error: output pointer is NULL." << endl;
  //   return -1;
  // }

  // if(input1->GetNumberOfPoints() == 0 || input2->GetNumberOfPoints() == 0) {
  //   cerr << "[ttkBottleneckDistance] Error: input has no points." << endl;
  //   return -1;
  // }
  // #endif

  auto outputCT1 = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  auto outputCT2 = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));
  auto outputCT3 = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(2)->Get(vtkDataObject::DATA_OBJECT()));

  // Wrap
  // this->setWrapper(this);
  this->setPersistencePercentThreshold(Tolerance);
  this->setPX(PX);
  this->setPY(PY);
  this->setPZ(PZ);
  this->setPE(PE);
  this->setPS(PS);

  CTPersistenceDiagram1_ = vtkUnstructuredGrid::SafeDownCast(input1);
  CTPersistenceDiagram2_ = vtkUnstructuredGrid::SafeDownCast(input2);

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

  // if(status != 0)
  //   return status;

  // Set output.
  outputCT1->ShallowCopy(CTPersistenceDiagram1_);
  outputCT2->DeepCopy(CTPersistenceDiagram2_);
  if(UseGeometricSpacing)
    translateSecondDiagram<dataType>(outputCT2, Spacing);

  if(UseOutputMatching)
    outputCT3->ShallowCopy(CTPersistenceDiagram3_);

  return status;
}
