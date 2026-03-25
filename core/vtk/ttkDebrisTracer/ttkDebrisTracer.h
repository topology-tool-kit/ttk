/// \ingroup vtk
/// \class ttkDebrisTracer
/// \author Théophane Loloum <theophane.loloum@gmail.com> 
/// \date March 2026
///
/// \brief TTK VTK-filter that takes an input tracking mesh and an input time-varying data 
/// set (represented by a list of scalar fields) and which computes linearizes, chains, and
/// surface statistics from tracked debris trajectories..
///
/// \param Input vtkUnstructuredGrid tracking mesh, i.e. ouput of TTKTrackingFromFields 
///	filter
///
/// \param Input time-dependent scalar field, either 2D or 3D, regular
/// grid or triangulation (vtkDataSet); time steps are obtained by
/// GetPointData()->GetArray(i) in increasing time order.
/// 
/// \param Output vtkTable : statistics on tracked object
/// \param Output vtkUnstructuredGrid : linearizes and fuses trajectories
/// \param Output time-dependent scalar field, pointDatai associated with a TrajId 
/// if belonging to a surface or -1
///
/// \sa ttk::DebrisTracer
/// \sa ttkAlgorithm

#pragma once

#include <ttkDebrisTracerModule.h>
#include <ttkAlgorithm.h>
#include <DebrisTracer.h>

class TTKDEBRISTRACER_EXPORT ttkDebrisTracer
  : public ttkAlgorithm  ,
    protected ttk::DebrisTracer 
{
private:
  int frameSurface{0}; 
  double errSurf{0.0};
  double filtreY{1.0};
  double maxRadius{40.0};
  double cosCol{0.96};
  int maxFrameDist{20};
  double spatialScale{1.0};
  double interFrame{1.0};
  bool convertDur{false};
  bool onlyFrameSurface{false};
  double minVx{0.0};
  double maxVx{0.0};
  bool extendTraj{false};
  bool enableFilteringMinVx{false};
  bool enableFilteringCosY{false};
  bool enableFilteringTimeOrigin{false};
  bool enableFilteringDuration{false};
  int duraMin{0};
  int xOrigin{0};
  int minTimeOrigin{0};
  int maxYTimeOrigin{0};
  int minYTimeOrigin{0};
  double persisThresh{0.0};
  int surfaceMethod{0};
  int maxSurfSize{10000};

public:
  /**
   * TODO 6: Automatically generate getters and setters of filter
   *         parameters via vtkMacros.
   */


  static ttkDebrisTracer *New();
  vtkTypeMacro(ttkDebrisTracer, ttkAlgorithm);

  vtkSetMacro(frameSurface, int);
  vtkGetMacro(frameSurface, int);

  vtkSetMacro(errSurf, double);
  vtkGetMacro(errSurf, double);

  vtkSetMacro(filtreY, double);
  vtkGetMacro(filtreY, double);

  vtkSetMacro(cosCol, double);
  vtkGetMacro(cosCol, double);

  vtkSetMacro(maxRadius, double);
  vtkGetMacro(maxRadius, double);
    
  vtkSetMacro(maxFrameDist, int);
  vtkGetMacro(maxFrameDist, int);

  vtkSetMacro(spatialScale, double);
  vtkGetMacro(spatialScale, double);

  vtkSetMacro(interFrame, double);
  vtkGetMacro(interFrame, double);

  vtkSetMacro(convertDur, bool);
  vtkGetMacro(convertDur, bool);

  vtkSetMacro(onlyFrameSurface, bool);
  vtkGetMacro(onlyFrameSurface, bool);
  
  vtkSetMacro(minVx, double);
  vtkGetMacro(minVx, double);

  vtkSetMacro(maxVx, double);
  vtkGetMacro(maxVx, double);

  vtkSetMacro(extendTraj, bool);
  vtkGetMacro(extendTraj, bool);

  vtkSetMacro(enableFilteringMinVx, bool);
  vtkGetMacro(enableFilteringMinVx, bool);

  vtkSetMacro(enableFilteringCosY, bool);
  vtkGetMacro(enableFilteringCosY, bool);

  vtkSetMacro(enableFilteringTimeOrigin, bool);
  vtkGetMacro(enableFilteringTimeOrigin, bool);

  vtkSetMacro(enableFilteringDuration, bool);
  vtkGetMacro(enableFilteringDuration, bool);

  vtkSetMacro(xOrigin, int);
  vtkGetMacro(xOrigin, int);

  vtkSetMacro(minTimeOrigin, int);
  vtkGetMacro(minTimeOrigin, int);
  
  vtkSetMacro(minYTimeOrigin, int);
  vtkGetMacro(minYTimeOrigin, int);

  vtkSetMacro(maxYTimeOrigin, int);
  vtkGetMacro(maxYTimeOrigin, int);

  vtkSetMacro(persisThresh, double);
  vtkGetMacro(persisThresh, double);

  vtkSetMacro(duraMin, int);
  vtkGetMacro(duraMin, int);

  vtkSetMacro(maxSurfSize, int);
  vtkGetMacro(maxSurfSize, int);

  vtkSetMacro(surfaceMethod, int);
  vtkGetMacro(surfaceMethod, int);
protected:

  ttkDebrisTracer();
  ~ttkDebrisTracer() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;


  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int computeAllGradientMagnitudes(
    vtkDataSet *inputDataSet,
    const std::vector<vtkDataArray *> &inputScalarFields,
    std::vector<std::vector<double>> &gradientNorms
  );


  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
