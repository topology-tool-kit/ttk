#ifndef _TTK_TRACKINGFROMF_H
#define _TTK_TRACKINGFROMF_H

#include <tuple>

#include <Wrapper.h>

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

#include <ttkTrackingFromPersistenceDiagrams.h>

#include <TrackingFromFields.h>
#include <TrackingFromPersistenceDiagrams.h>

#include <algorithm>
#include <string>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTrackingFromFields
#else
class ttkTrackingFromFields
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkTrackingFromFields *New();

  vtkTypeMacro(ttkTrackingFromFields, vtkDataSetAlgorithm);

  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  vtkSetMacro(Sampling, int);
  vtkGetMacro(Sampling, int);

  vtkSetMacro(StartTimestep, int);
  vtkGetMacro(StartTimestep, int);

  vtkSetMacro(EndTimestep, int);
  vtkGetMacro(EndTimestep, int);

  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);

  vtkSetMacro(PX, double);
  vtkGetMacro(PX, double);

  vtkSetMacro(PY, double);
  vtkGetMacro(PY, double);

  vtkSetMacro(PZ, double);
  vtkGetMacro(PZ, double);

  vtkSetMacro(PE, double);
  vtkGetMacro(PE, double);

  vtkSetMacro(PS, double);
  vtkGetMacro(PS, double);

  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);

  vtkSetMacro(WassersteinMetric, std::string);
  vtkGetMacro(WassersteinMetric, std::string);

  vtkSetMacro(DistanceAlgorithm, std::string);
  vtkGetMacro(DistanceAlgorithm, std::string);

  vtkSetMacro(PVAlgorithm, int);
  vtkGetMacro(PVAlgorithm, int);

  vtkSetMacro(UseGeometricSpacing, int);
  vtkGetMacro(UseGeometricSpacing, int);

  vtkSetMacro(Spacing, double);
  vtkGetMacro(Spacing, double);
  vtkSetMacro(DoPostProc, int);
  vtkGetMacro(DoPostProc, int);

  vtkSetMacro(PostProcThresh, double);
  vtkGetMacro(PostProcThresh, double);

protected:
  ttkTrackingFromFields() {
    outputMesh_ = nullptr;
    UseAllCores = false;

    DistanceAlgorithm = "ttk";
    PVAlgorithm = -1;
    Alpha = 1.0;
    WassersteinMetric = "2";
    UseGeometricSpacing = false;
    Is3D = true;
    Spacing = 1.0;

    DoPostProc = false;
    PostProcThresh = 0;

    Sampling = 1;
    StartTimestep = 0;
    EndTimestep = -1;

    Tolerance = 1;
    PX = 1;
    PY = 1;
    PZ = 0;
    PE = 0;
    PS = 0;

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }

  ~ttkTrackingFromFields() {
    if(outputMesh_)
      outputMesh_->Delete();
  }

  TTK_SETUP();

  virtual int FillOutputPortInformation(int port,
                                        vtkInformation *info) override;

private:
  // Sampling config.
  int StartTimestep;
  int EndTimestep;
  int Sampling;

  // Filtering config.
  double Tolerance;
  double PX;
  double PY;
  double PZ;
  double PE;
  double PS;

  // Bottleneck config.
  bool UseGeometricSpacing;
  bool Is3D;
  bool DoPostProc;
  double PostProcThresh;
  double Spacing;
  double Alpha;
  std::string DistanceAlgorithm;
  int PVAlgorithm;
  std::string WassersteinMetric;

  vtkUnstructuredGrid *outputMesh_;

  ttk::Triangulation *internalTriangulation_;
  ttkTriangulation triangulation_;
  ttk::TrackingFromFields trackingF_;
  ttk::TrackingFromPersistenceDiagrams tracking_;

  template <typename dataType>
  int trackWithPersistenceMatching(
    vtkDataSet *input,
    vtkUnstructuredGrid *output,
    std::vector<vtkDataArray *> inputScalarFields);
};

// (*) Persistence-driven approach
template <typename dataType>
int ttkTrackingFromFields::trackWithPersistenceMatching(
  vtkDataSet *input,
  vtkUnstructuredGrid *output,
  std::vector<vtkDataArray *> inputScalarFields) {
  unsigned long fieldNumber = inputScalarFields.size();

  // 0. get data
  trackingF_.setThreadNumber(ThreadNumber);
  trackingF_.setTriangulation(internalTriangulation_);
  std::vector<void *> inputFields(fieldNumber);
  for(int i = 0; i < (int)fieldNumber; ++i)
    inputFields[i] = inputScalarFields[i]->GetVoidPointer(0);
  trackingF_.setInputScalars(inputFields);

  // 0'. get offsets
  auto numberOfVertices = (int)input->GetNumberOfPoints();
  vtkIdTypeArray *offsets_ = vtkIdTypeArray::New();
  offsets_->SetNumberOfComponents(1);
  offsets_->SetNumberOfTuples(numberOfVertices);
  offsets_->SetName("OffsetScalarField");
  for(int i = 0; i < numberOfVertices; ++i)
    offsets_->SetTuple1(i, i);
  trackingF_.setInputOffsets(offsets_->GetVoidPointer(0));

  // 1. get persistence diagrams.
  std::vector<std::vector<diagramTuple>> persistenceDiagrams(
    fieldNumber, std::vector<diagramTuple>());

  trackingF_.performDiagramComputation<dataType>(
    (int)fieldNumber, persistenceDiagrams, this);

  // 2. call feature tracking with threshold.
  std::vector<std::vector<matchingTuple>> outputMatchings(
    fieldNumber - 1, std::vector<matchingTuple>());

  double spacing = Spacing;
  std::string algorithm = DistanceAlgorithm;
  double alpha = Alpha;
  double tolerance = Tolerance;
  bool is3D = true; // Is3D;
  std::string wasserstein = WassersteinMetric;

  tracking_.setThreadNumber(ThreadNumber);
  tracking_.performMatchings<dataType>(
    (int)fieldNumber, persistenceDiagrams, outputMatchings,
    algorithm, // Not from paraview, from enclosing tracking plugin
    wasserstein, tolerance, is3D,
    alpha, // Blending
    PX, PY, PZ, PS, PE, // Coefficients
    this // Wrapper for accessing threadNumber
  );

  outputMesh_ = vtkUnstructuredGrid::New();
  vtkUnstructuredGrid *outputMesh = vtkUnstructuredGrid::SafeDownCast(output);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkDoubleArray> persistenceScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> valueScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray> matchingIdScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> lengthScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> timeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> componentIds
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> pointTypeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  persistenceScalars->SetName("Cost");
  valueScalars->SetName("Scalar");
  matchingIdScalars->SetName("MatchingIdentifier");
  lengthScalars->SetName("ComponentLength");
  timeScalars->SetName("TimeStep");
  componentIds->SetName("ConnectedComponentId");
  pointTypeScalars->SetName("CriticalType");

  // (+ vertex id)
  std::vector<trackingTuple> trackingsBase;
  tracking_.setThreadNumber(ThreadNumber);
  tracking_.performTracking<dataType>(
    persistenceDiagrams, outputMatchings, trackingsBase);

  std::vector<std::set<int>> trackingTupleToMerged(
    trackingsBase.size(), std::set<int>());

  if(DoPostProc)
    tracking_.performPostProcess<dataType>(persistenceDiagrams, trackingsBase,
                                           trackingTupleToMerged,
                                           PostProcThresh);

  bool useGeometricSpacing = UseGeometricSpacing;

  // Build mesh.
  ttkTrackingFromPersistenceDiagrams::buildMesh(
    trackingsBase, outputMatchings, persistenceDiagrams, useGeometricSpacing,
    spacing, DoPostProc, trackingTupleToMerged, points, persistenceDiagram,
    persistenceScalars, valueScalars, matchingIdScalars, lengthScalars,
    timeScalars, componentIds, pointTypeScalars);

  outputMesh_->ShallowCopy(persistenceDiagram);
  outputMesh->ShallowCopy(outputMesh_);

  return 0;
}

#endif // _TTK_TRACKINGFROMF_H
