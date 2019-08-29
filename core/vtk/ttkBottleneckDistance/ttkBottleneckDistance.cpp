#include "ttkBottleneckDistance.h"

vtkStandardNewMacro(ttkBottleneckDistance)

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

  bottleneckDistance_.setPersistencePercentThreshold(Tolerance);
  bottleneckDistance_.setPX(PX);
  bottleneckDistance_.setPY(PY);
  bottleneckDistance_.setPZ(PZ);
  bottleneckDistance_.setPE(PE);
  bottleneckDistance_.setPS(PS);
  bottleneckDistance_.setCTDiagram1(&CTDiagram1);
  bottleneckDistance_.setCTDiagram2(&CTDiagram2);

  std::string wassersteinMetric = WassersteinMetric;
  bottleneckDistance_.setWasserstein(wassersteinMetric);
  std::string algorithm = DistanceAlgorithm;
  bottleneckDistance_.setAlgorithm(algorithm);
  int pvAlgorithm = PVAlgorithm;
  bottleneckDistance_.setPVAlgorithm(pvAlgorithm);
  bottleneckDistance_.setThreadNumber(ThreadNumber);

  // Empty matchings.
  auto matchings = new std::vector<diagramTuple>();
  bottleneckDistance_.setOutputMatchings(matchings);

  // Exec.
  bool usePersistenceMetric = UsePersistenceMetric;
  // double alpha = Alpha;
  status = bottleneckDistance_.execute<dataType>(usePersistenceMetric);

  if(status != 0) {
    return status;
  }

  return 0;
}

int ttkBottleneckDistance::doIt(std::vector<vtkDataSet *> &inputs,
                                std::vector<vtkDataSet *> &outputs) {
  using dataType = double;

  int benchmarkSize = BenchmarkSize;
  bool benchmark = benchmarkSize > 0;
  if(benchmark) {
    return doBenchmark();
  }

  // Prepare IO
  vtkDataSet *input1 = inputs[0];
  vtkDataSet *input2 = inputs[1];
  vtkDataSet *output1 = outputs[0];
  vtkDataSet *output2 = outputs[1];
  vtkDataSet *output3 = outputs[2];

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input1 || !input2) {
    cerr << "[ttkBottleneckDistance] Error: input pointer is NULL." << endl;
    return -1;
  }

  if(!output1 || !output2 || !output3) {
    cerr << "[ttkBottleneckDistance] Error: output pointer is NULL." << endl;
    return -1;
  }

  if(input1->GetNumberOfPoints() == 0 || input2->GetNumberOfPoints() == 0) {
    cerr << "[ttkBottleneckDistance] Error: input has no points." << endl;
    return -1;
  }
#endif

  vtkUnstructuredGrid *outputCT1 = vtkUnstructuredGrid::SafeDownCast(output1);
  vtkUnstructuredGrid *outputCT2 = vtkUnstructuredGrid::SafeDownCast(output2);
  vtkUnstructuredGrid *outputCT3 = vtkUnstructuredGrid::SafeDownCast(output3);

  // Wrap
  bottleneckDistance_.setWrapper(this);
  bottleneckDistance_.setPersistencePercentThreshold(Tolerance);
  bottleneckDistance_.setPX(PX);
  bottleneckDistance_.setPY(PY);
  bottleneckDistance_.setPZ(PZ);
  bottleneckDistance_.setPE(PE);
  bottleneckDistance_.setPS(PS);

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

  bottleneckDistance_.setCTDiagram1(&CTDiagram1);
  bottleneckDistance_.setCTDiagram2(&CTDiagram2);

  std::string wassersteinMetric = WassersteinMetric;
  bottleneckDistance_.setWasserstein(wassersteinMetric);
  std::string algorithm = DistanceAlgorithm;
  bottleneckDistance_.setAlgorithm(algorithm);
  int pvAlgorithm = PVAlgorithm;
  bottleneckDistance_.setPVAlgorithm(pvAlgorithm);

  // Empty matchings.
  std::vector<matchingTuple> matchings;
  bottleneckDistance_.setOutputMatchings(&matchings);

  // Exec.
  bool usePersistenceMetric = UsePersistenceMetric;
  // double alpha = Alpha;
  status = bottleneckDistance_.execute<dataType>(usePersistenceMetric);
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
