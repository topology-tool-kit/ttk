#pragma once

#include <tuple>

#include <ttkAlgorithm.h>

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

// VTK Module
#include <ttkTrackingFromFieldsModule.h>

#include <TrackingFromFields.h>
#include <TrackingFromPersistenceDiagrams.h>

#include <algorithm>
#include <string>

class TTKTRACKINGFROMFIELDS_EXPORT ttkTrackingFromFields
  : public ttkAlgorithm,
    protected ttk::TrackingFromFields {

public:
  static ttkTrackingFromFields *New();

  vtkTypeMacro(ttkTrackingFromFields, ttkAlgorithm);

  // void SetDebugLevel(int debugLevel) {
  //   setDebugLevel(debugLevel);
  //   Modified();
  // }

  // void SetThreadNumber(int threadNumber) {
  //   ThreadNumber = threadNumber;
  //   SetThreads();
  // }

  // void SetUseAllCores(bool onOff) {
  //   UseAllCores = onOff;
  //   SetThreads();
  // }

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

  vtkSetMacro(UseGeometricSpacing, bool);
  vtkGetMacro(UseGeometricSpacing, bool);

  vtkSetMacro(Spacing, double);
  vtkGetMacro(Spacing, double);
  vtkSetMacro(DoPostProc, bool);
  vtkGetMacro(DoPostProc, bool);

  vtkSetMacro(PostProcThresh, double);
  vtkGetMacro(PostProcThresh, double);

protected:
  ttkTrackingFromFields() {

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }

  ~ttkTrackingFromFields() override {
    if(outputMesh_)
      outputMesh_->Delete();
  }

  virtual int FillInputPortInformation(int port, vtkInformation *info) override;

  virtual int FillOutputPortInformation(int port,
                                        vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // Sampling config.
  int StartTimestep{0};
  int EndTimestep{-1};
  int Sampling{1};

  // Filtering config.
  double Tolerance{1};
  double PX{1};
  double PY{1};
  double PZ{0};
  double PE{0};
  double PS{0};

  // Bottleneck config.
  bool UseGeometricSpacing{false};
  bool Is3D{true};
  bool DoPostProc{false};
  double PostProcThresh{0.0};
  double Spacing{1.0};
  double Alpha{1.0};
  std::string DistanceAlgorithm{"ttk"};
  int PVAlgorithm{-1};
  std::string WassersteinMetric{"2"};

  vtkUnstructuredGrid *outputMesh_{nullptr};

  // ttk::Triangulation *internalTriangulation_;
  // ttkTriangulation triangulation_;
  ttk::TrackingFromFields trackingF_;
  ttk::TrackingFromPersistenceDiagrams tracking_;

  template <class dataType, class triangulationType = AbstractTriangulation>
  int trackWithPersistenceMatching(
    vtkDataSet *input,
    vtkUnstructuredGrid *output,
    std::vector<vtkDataArray *> inputScalarFields,
    triangulationType *triangulation);
};

// (*) Persistence-driven approach
template <class dataType>
int ttkTrackingFromFields::trackWithPersistenceMatching(
  vtkDataSet *input,
  vtkUnstructuredGrid *output,
  std::vector<vtkDataArray *> inputScalarFields) {
  unsigned long fieldNumber = inputScalarFields.size();

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  using trackingTuple = ttk::trackingTuple;

  // 0. get data
  // this->setThreadNumber(ThreadNumber);
  // this->setTriangulation(internalTriangulation_);
  std::vector<void *> inputFields(fieldNumber);
  for(int i = 0; i < (int)fieldNumber; ++i)
    inputFields[i] = inputScalarFields[i]->GetVoidPointer(0);
  this->setInputScalars(inputFields);

  // 0'. get offsets
  auto numberOfVertices = (int)input->GetNumberOfPoints();
  vtkIdTypeArray *offsets_ = vtkIdTypeArray::New();
  offsets_->SetNumberOfComponents(1);
  offsets_->SetNumberOfTuples(numberOfVertices);
  offsets_->SetName("OffsetScalarField");
  for(int i = 0; i < numberOfVertices; ++i)
    offsets_->SetTuple1(i, i);
  this->setInputOffsets(offsets_->GetVoidPointer(0));

  // 1. get persistence diagrams.
  std::vector<std::vector<diagramTuple>> persistenceDiagrams(
    fieldNumber, std::vector<diagramTuple>());

  int status = 0;
  // ttkVtkTemplateMacro(
  //   triangulation->getType(), inputScalarFields[0]->GetDataType(),
  //   (status = this->performDiagramComputation<dataType, TTK_TT>(
  //      (int)fieldNumber, persistenceDiagrams,
  //      (TTK_TT *)triangulation->getData())));
  ttkTemplateMacro(triangulation->getType(),
                   (status = this->performDiagramComputation<dataType, TTK_TT>(
                      (int)fieldNumber, persistenceDiagrams,
                      (TTK_TT *)triangulation->getData())));
  if(status < 0)
    return 0;

  // 2. call feature tracking with threshold.
  std::vector<std::vector<matchingTuple>> outputMatchings(
    fieldNumber - 1, std::vector<matchingTuple>());

  double spacing = Spacing;
  std::string algorithm = DistanceAlgorithm;
  double alpha = Alpha;
  double tolerance = Tolerance;
  bool is3D = true; // Is3D;
  std::string wasserstein = WassersteinMetric;

  // tracking_.setThreadNumber(ThreadNumber);
  tracking_.performMatchings<dataType>(
    (int)fieldNumber, persistenceDiagrams, outputMatchings,
    algorithm, // Not from paraview, from enclosing tracking plugin
    wasserstein, tolerance, is3D,
    alpha, // Blending
    PX, PY, PZ, PS, PE // Coefficients
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
  // tracking_.setThreadNumber(ThreadNumber);
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
  // ttkTrackingFromPersistenceDiagrams::buildMesh(
  //   trackingsBase, outputMatchings, persistenceDiagrams, useGeometricSpacing,
  //   spacing, DoPostProc, trackingTupleToMerged, points, persistenceDiagram,
  //   persistenceScalars, valueScalars, matchingIdScalars, lengthScalars,
  //   timeScalars, componentIds, pointTypeScalars);

  outputMesh_->ShallowCopy(persistenceDiagram);
  outputMesh->ShallowCopy(outputMesh_);

  return 0;
}

