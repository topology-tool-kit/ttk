#include                  "ttkBottleneckDistance.h"

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkBottleneckDistance)

#ifndef macroDiagramTuple
#define macroDiagramTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::NodeType, ttk::ftm::idVertex, \
  ttk::ftm::NodeType, double, ttk::ftm::idVertex, \
  double, float, float, float, double, float, float, float>
#endif
#ifndef macroMatchingTuple
#define macroMatchingTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::idVertex, double>
#endif
#ifndef benchDiagramTuple
#define benchDiagramTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::NodeType, ttk::ftm::idVertex, \
  ttk::ftm::NodeType, double, ttk::ftm::idVertex, \
  double, float, float, float, double, float, float, float>
#endif

int ttkBottleneckDistance::doBenchmark()
{
  std::vector<benchDiagramTuple>* CTDiagram1 = new std::vector<benchDiagramTuple>();
  std::vector<benchDiagramTuple>* CTDiagram2 = new std::vector<benchDiagramTuple>();

  int benchmarkSize = BenchmarkSize;
  int status = 0;
  status = generatePersistenceDiagram<double>(CTDiagram1, benchmarkSize);
  if (status < 0) return status;
  status = generatePersistenceDiagram<double>(CTDiagram2, 4 * benchmarkSize);
  if (status < 0) return status;

  bottleneckDistance_.setPersistencePercentThreshold(Tolerance);
  bottleneckDistance_.setPX(PX);
  bottleneckDistance_.setPY(PY);
  bottleneckDistance_.setPZ(PZ);
  bottleneckDistance_.setPE(PE);
  bottleneckDistance_.setPS(PS);
  bottleneckDistance_.setCTDiagram1(CTDiagram1);
  bottleneckDistance_.setCTDiagram2(CTDiagram2);

  std::string wassersteinMetric = WassersteinMetric;
  bottleneckDistance_.setWasserstein(wassersteinMetric);
  std::string algorithm = DistanceAlgorithm;
  bottleneckDistance_.setAlgorithm(algorithm);
  int pvAlgorithm = PVAlgorithm;
  bottleneckDistance_.setPVAlgorithm(pvAlgorithm);
  bottleneckDistance_.setThreadNumber(ThreadNumber);

  // Empty matchings.
  std::vector<benchDiagramTuple>* matchings = new std::vector<benchDiagramTuple>();
  bottleneckDistance_.setOutputMatchings(matchings);

  // Exec.
  bool usePersistenceMetric = UsePersistenceMetric;
  double alpha = Alpha;
  status = bottleneckDistance_.execute<double>(usePersistenceMetric, alpha, Is3D);

  if (status != 0) { return status; }

  return 0;
}

int ttkBottleneckDistance::doIt(
    std::vector<vtkDataSet *> &inputs,
    std::vector<vtkDataSet *> &outputs)
{

  int benchmarkSize = BenchmarkSize;
  bool benchmark = benchmarkSize > 0;
  if (benchmark) {
    return doBenchmark();
  }

  // Prepare IO
  vtkDataSet *input1 = inputs[0];
  vtkDataSet *input2 = inputs[1];
  vtkDataSet *output1 = outputs[0];
  vtkDataSet *output2 = outputs[1];
  vtkDataSet *output3 = outputs[2];

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

  if (!CTPersistenceDiagram1_ || !CTPersistenceDiagram2_ ||
      !outputCT3) return -1;

  int dataType1 = CTPersistenceDiagram1_->GetCellData()->GetArray("Persistence")->GetDataType();
  int dataType2 = CTPersistenceDiagram2_->GetCellData()->GetArray("Persistence")->GetDataType();
  if (dataType1 != dataType2) return -1;

  // Call package
  int status = 0;

//  switch (dataType1) {
//    vtkTemplateMacro(({
  std::vector<macroDiagramTuple>* CTDiagram1 = new std::vector<macroDiagramTuple>();

  std::vector<macroDiagramTuple>* CTDiagram2 = new std::vector<macroDiagramTuple>();

  status = getPersistenceDiagram<double>(
      CTDiagram1, CTPersistenceDiagram1_, Spacing, 0);

  status = getPersistenceDiagram<double>(
      CTDiagram2, CTPersistenceDiagram2_, Spacing, 1);

  bottleneckDistance_.setCTDiagram1(CTDiagram1);
  bottleneckDistance_.setCTDiagram2(CTDiagram2);

  std::string wassersteinMetric = WassersteinMetric;
  bottleneckDistance_.setWasserstein(wassersteinMetric);
  std::string algorithm = DistanceAlgorithm;
  bottleneckDistance_.setAlgorithm(algorithm);
  int pvAlgorithm = PVAlgorithm;
  bottleneckDistance_.setPVAlgorithm(pvAlgorithm);

  // Empty matchings.
  std::vector<macroMatchingTuple>* matchings = new std::vector<macroMatchingTuple>();
  bottleneckDistance_.setOutputMatchings(matchings);

  // Exec.
  bool usePersistenceMetric = UsePersistenceMetric;
  double alpha = Alpha;
  status = bottleneckDistance_.execute<double>(usePersistenceMetric, alpha, Is3D);

  // Apply results to outputs 0 and 1.
  status = augmentPersistenceDiagrams<double>(
      CTDiagram1,
      CTDiagram2,
      matchings,
      CTPersistenceDiagram1_,
      CTPersistenceDiagram2_);

  bool useOutputMatching = UseOutputMatching;
  bool useGeometricSpacing = UseGeometricSpacing;

  // Apply results to output 2.
  if (useOutputMatching) {
    status = getMatchingMesh<double>(
        CTDiagram1, CTDiagram2, matchings,
        useGeometricSpacing, Spacing);
  }

  if (status != 0) { return status; }
//    }));
//  }

  // Set output.
  outputCT1->ShallowCopy(CTPersistenceDiagram1_);
  outputCT2->ShallowCopy(CTPersistenceDiagram2_);
  if (UseOutputMatching)
    outputCT3->ShallowCopy(CTPersistenceDiagram3_);

  return status;
}
