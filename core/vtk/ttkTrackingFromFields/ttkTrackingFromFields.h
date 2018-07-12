#ifndef _TTK_TRACKING_H
#define _TTK_TRACKING_H

#include                  <tuple>

#include                  <Wrapper.h>

#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkInformationVector.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkPoints.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkTable.h>
#include                  <vtkUnstructuredGrid.h>
#include                  <vtkCellType.h>
#include                  <vtkCellData.h>

#include                  <ttkBottleneckDistance.h>
#include                  <ttkPersistenceDiagram.h>
#include                  <ttkFTMTree.h>
#include                  <ttkTrackingFromPersistenceDiagrams.h>

#include                  <algorithm>
#include                  <string>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTrackingFromFields
#else
class ttkTrackingFromFields
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper
{

  using idVertex = int;

  public:

    static ttkTrackingFromFields* New();

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

    vtkSetMacro(Is3D, int);
    vtkGetMacro(Is3D, int);

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
      WassersteinMetric = "1";
      UseGeometricSpacing = false;
      Is3D = false;
      Spacing = 1.0;

      DoPostProc = false;
      PostProcThresh = 0;

      Sampling = 1;
      StartTimestep = 0;
      EndTimestep = -1;

      Tolerance = 1;
      PX = 1;
      PY = 1;
      PZ = 1;
      PE = 1;
      PS = 1;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }

    ~ttkTrackingFromFields() {
      if (outputMesh_)
        outputMesh_->Delete();
    }

    TTK_SETUP();

    virtual int FillOutputPortInformation(int port, vtkInformation* info);

  private:

    // Sampling config.
    int                   StartTimestep;
    int                   EndTimestep;
    int                   Sampling;

    // Filtering config.
    double                Tolerance;
    double                PX;
    double                PY;
    double                PZ;
    double                PE;
    double                PS;

    // Bottleneck config.
    bool                  UseGeometricSpacing;
    bool                  Is3D;
    bool                  DoPostProc;
    double                PostProcThresh;
    double                Spacing;
    double                Alpha;
    std::string           DistanceAlgorithm;
    int                   PVAlgorithm;
    std::string           WassersteinMetric;

    vtkUnstructuredGrid   *outputMesh_;

    ttk::Triangulation    *internalTriangulation_;
    ttkTriangulation      triangulation_;

    int trackWithPersistenceMatching(
        vtkDataSet *input,
        vtkUnstructuredGrid *output,
        std::vector<vtkDataArray*> inputScalarFields);

};

#endif // _TTK_TRACKING_H
