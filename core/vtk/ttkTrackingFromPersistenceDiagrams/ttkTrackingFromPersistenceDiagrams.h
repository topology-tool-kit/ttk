#pragma once

#include <TrackingFromPersistenceDiagrams.h>
#include <ttkAlgorithm.h>

// VTK Module
#include <ttkTrackingFromPersistenceDiagramsModule.h>
#include <vtkUnstructuredGrid.h>

class vtkUnstructuredGrid;

class TTKTRACKINGFROMPERSISTENCEDIAGRAMS_EXPORT
  ttkTrackingFromPersistenceDiagrams
  : public ttkAlgorithm,
    protected ttk::TrackingFromPersistenceDiagrams {

public:
  static ttkTrackingFromPersistenceDiagrams *New();

  vtkTypeMacro(ttkTrackingFromPersistenceDiagrams, ttkAlgorithm);

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

  vtkSetMacro(WassersteinMetric, const std::string &);
  vtkGetMacro(WassersteinMetric, std::string);

  vtkSetMacro(DistanceAlgorithm, const std::string &);
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

  static int buildMesh(
    const std::vector<ttk::trackingTuple> &trackings,
    const std::vector<std::vector<ttk::MatchingType>> &outputMatchings,
    const std::vector<ttk::DiagramType> &inputPersistenceDiagrams,
    const bool useGeometricSpacing,
    const double spacing,
    const bool doPostProc,
    const std::vector<std::set<int>> &trackingTupleToMerged,
    vtkPoints *points,
    vtkUnstructuredGrid *persistenceDiagram,
    vtkDoubleArray *persistenceScalars,
    vtkDoubleArray *valueScalars,
    vtkIntArray *matchingIdScalars,
    vtkIntArray *lengthScalars,
    vtkIntArray *timeScalars,
    vtkIntArray *componentIds,
    vtkIntArray *pointTypeScalars,
    const ttk::Debug &dbg);

protected:
  ttkTrackingFromPersistenceDiagrams();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // Input bottleneck config.
  bool UseGeometricSpacing{false};
  bool DoPostProc{false};
  double PostProcThresh{0.0};
  double Spacing{1.0};
  double Tolerance{1.0};
  double PX{1};
  double PY{1};
  double PZ{1};
  double PE{1};
  double PS{1};
  std::string DistanceAlgorithm{"ttk"};
  int PVAlgorithm{-1};
  std::string WassersteinMetric{"1"};
};
