#pragma once

#include <TrackingFromPersistenceDiagrams.h>
#include <ttkAlgorithm.h>

// VTK Module
#include <ttkTrackingFromPersistenceDiagramsModule.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
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

  template <class triangulationType>
  static int
    buildMesh(const triangulationType *triangulation,
              const std::vector<ttk::trackingTuple> &trackings,
              const std::vector<std::vector<double>> &allTrackingsCosts,
              const std::vector<double> &allTrackingsPersistence,
              const bool useGeometricSpacing,
              const double spacing,
              vtkPoints *points,
              vtkUnstructuredGrid *outputMesh,
              vtkIntArray *pointsCriticalType,
              vtkIntArray *timeScalars,
              vtkIntArray *lengthScalars,
              vtkIntArray *globalVertexIds,
              vtkIntArray *connectedComponentIds,
              vtkDoubleArray *costs,
              vtkDoubleArray *averagePersistence,
              unsigned int *sizes) {

    int pointCpt = 0;
    int edgeCpt = 0;
    for(unsigned int i = 0; i < trackings.size(); i++) {
      ttk::CriticalType currentType = ttk::CriticalType::Local_minimum;
      if(i < sizes[0])
        currentType = ttk::CriticalType::Local_maximum;
      else if(i < sizes[1] && i >= sizes[0])
        currentType = ttk::CriticalType::Saddle1;
      else if(i < sizes[2] && i >= sizes[1])
        currentType = ttk::CriticalType::Saddle2;
      int startTime = std::get<0>(trackings[i]);
      std::vector<ttk::SimplexId> chain = std::get<2>(trackings[i]);

      float x = 0;
      float y = 0;
      float z = 0;
      triangulation->getVertexPoint(chain[0], x, y, z);
      if(useGeometricSpacing)
        z += startTime * spacing;
      points->InsertNextPoint(x, y, z);
      globalVertexIds->InsertTuple1(pointCpt, (int)chain[0]);
      pointsCriticalType->InsertTuple1(pointCpt, (int)currentType);
      timeScalars->InsertTuple1(pointCpt, startTime);
      vtkIdType edge[2];
      for(unsigned int j = 1; j < chain.size(); j++) {
        triangulation->getVertexPoint(chain[j], x, y, z);
        if(useGeometricSpacing)
          z += (j + startTime) * spacing;
        edge[0] = pointCpt;
        pointCpt++;
        edge[1] = pointCpt;
        points->InsertNextPoint(x, y, z);
        globalVertexIds->InsertTuple1(pointCpt, (int)chain[j]);
        outputMesh->InsertNextCell(VTK_LINE, 2, edge);
        pointsCriticalType->InsertTuple1(pointCpt, (int)currentType);
        timeScalars->InsertTuple1(pointCpt, startTime + j);
        lengthScalars->InsertTuple1(edgeCpt, chain.size() - 1);
        connectedComponentIds->InsertTuple1(edgeCpt, i);
        costs->InsertTuple1(edgeCpt, allTrackingsCosts[i][j - 1]);
        averagePersistence->InsertTuple1(edgeCpt, allTrackingsPersistence[i]);
        edgeCpt++;
      }
      pointCpt++;
    }

    outputMesh->SetPoints(points);
    outputMesh->GetCellData()->AddArray(lengthScalars);
    outputMesh->GetCellData()->AddArray(connectedComponentIds);
    outputMesh->GetCellData()->AddArray(averagePersistence);
    outputMesh->GetCellData()->AddArray(costs);
    outputMesh->GetPointData()->AddArray(pointsCriticalType);
    outputMesh->GetPointData()->AddArray(timeScalars);
    outputMesh->GetPointData()->AddArray(globalVertexIds);

    return 0;
  }

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
