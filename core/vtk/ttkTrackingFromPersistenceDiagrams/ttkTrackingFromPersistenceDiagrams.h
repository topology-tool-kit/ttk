#ifndef _TTK_TRACKINGFROMP_H
#define _TTK_TRACKINGFROMP_H

#include                  <tuple>

#include                  <TrackingFromPersistenceDiagrams.h>
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
#include                  <vtkSmartPointer.h>
#include                  <vtkTable.h>
#include                  <vtkUnstructuredGrid.h>
#include                  <vtkCellType.h>
#include                  <vtkCellData.h>
#include                  <vtkIndent.h>

#include                  <ttkBottleneckDistance.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTrackingFromPersistenceDiagrams
#else
class ttkTrackingFromPersistenceDiagrams
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper
{

  public:

    using dataType = double;

    static ttkTrackingFromPersistenceDiagrams* New();

    vtkTypeMacro(ttkTrackingFromPersistenceDiagrams, vtkDataSetAlgorithm);

    vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = ttk::OsCall::getNumberOfCores();
      }
      Modified();
    }

    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }

    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }

    vtkSetMacro(Tolerance, double);
    vtkGetMacro(Tolerance, double);

    vtkSetMacro(Alpha, double);
    vtkGetMacro(Alpha, double);

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

    static int buildMesh(
      std::vector<trackingTuple>& trackings,
      std::vector<std::vector<matchingTuple>*>*& outputMatchings,
      std::vector<std::vector<diagramTuple>*>*& inputPersistenceDiagrams,
      bool useGeometricSpacing,
      double spacing,
      bool doPostProc,
      std::vector<std::set<int>>*& trackingTupleToMerged,
      vtkSmartPointer<vtkPoints>& points,
      vtkSmartPointer<vtkUnstructuredGrid>& persistenceDiagram,
      vtkSmartPointer<vtkDoubleArray>& persistenceScalars,
      vtkSmartPointer<vtkDoubleArray>& valueScalars,
      vtkSmartPointer<vtkIntArray>& matchingIdScalars,
      vtkSmartPointer<vtkIntArray>& timeScalars,
      vtkSmartPointer<vtkIntArray>& componentIds,
      vtkSmartPointer<vtkIntArray>& pointTypeScalars);

  protected:

    ttkTrackingFromPersistenceDiagrams();

    ~ttkTrackingFromPersistenceDiagrams();

    int FillInputPortInformation(int port, vtkInformation *info);
    int FillOutputPortInformation(int port, vtkInformation *info);

    int RequestData(vtkInformation *request,
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);

  private:

    // Input bottleneck config.
    bool                  UseGeometricSpacing;
    bool                  Is3D;
    bool                  DoPostProc;
    double                PostProcThresh;
    double                Spacing;
    double                Alpha;
    double                Tolerance;
    double                PX;
    double                PY;
    double                PZ;
    double                PE;
    double                PS;
    std::string           DistanceAlgorithm;
    int                   PVAlgorithm;
    std::string           WassersteinMetric;

    bool                  UseAllCores;
    int                   ThreadNumber;
    vtkUnstructuredGrid   *outputMesh_;

    int doIt(vtkDataSet **input,
             vtkUnstructuredGrid *outputMean,
             int numInputs);

    bool needsToAbort();

    int updateProgress(const float &progress);

    ttk::TrackingFromPersistenceDiagrams tracking_;
};

#endif // _TTK_TRACKINGFROMP_H
