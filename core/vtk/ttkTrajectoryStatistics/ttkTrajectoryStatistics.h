/// TODO 4: Provide your information and **update** the documentation (in
/// particular regarding the order convention if input arrays need to be
/// specified with the standard VTK call SetInputArrayToProcess()).
///
/// \ingroup vtk
/// \class ttkTrajectoryStatistics
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ttk::TrajectoryStatistics module.
///
/// This VTK filter uses the ttk::TrajectoryStatistics module to compute an averaging of
/// the data values of an input point data array defined on the input
/// vtkDataSet.
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the corresponding standalone program for a usage example:
///   - standalone/TrajectoryStatistics/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::TrajectoryStatistics
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkTrajectoryStatisticsModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

/* Note on including VTK modules
 *
 * Each VTK module that you include a header from needs to be specified in this
 * module's vtk.module file, either in the DEPENDS or PRIVATE_DEPENDS (if the
 * header is included in the cpp file only) sections.
 *
 * In order to find the corresponding module, check its location within the VTK
 * source code. The VTK module name is composed of the path to the header. You
 * can also find the module name within the vtk.module file located in the same
 * directory as the header file.
 *
 * For example, vtkSphereSource.h is located in directory VTK/Filters/Sources/,
 * so its corresponding VTK module is called VTK::FiltersSources. In this case,
 * the vtk.module file would need to be extended to
 *
 * NAME
 *   ttkTrajectoryStatistics
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <TrajectoryStatistics.h>

class TTKTRAJECTORYSTATISTICS_EXPORT ttkTrajectoryStatistics
  : public ttkAlgorithm  ,
    protected ttk::TrajectoryStatistics 
{
private:
  int frameSurface{0}; 
  double errSurf{0.0};
  double filtreY{1.0};
  double maxRadus{40.0};
  double cosCol{0.96};
  int maxFrameDist{20};
  double spatialScale{1.0};
  double interFrame{1.0};
  bool convertDur{false};
  bool onlyFrameSurface{false};
  double minVx{0.0};
  double maxVx{0.0};
  double extendTraj{false};
  int enableFilteringMinVx{0};
  int enableFilteringCosY{0};
  int enableFilteringTimeOrigin{0};
  int enableFilteringDuration{0};
  int duraMin{0};
  int xOrigin{0};
  int minTimeOrigin{0};
  int maxYTimeOrigin{0};
  int minYTimeOrigin{0};
  int maxX{-1};
  int maxY{-1};
  int minY{-1};
  int minX{-1};
  double persisThresh{0.0};
  int surfaceMethod{0};
  int maxSurfSize{10000};

public:
  /**
   * TODO 6: Automatically generate getters and setters of filter
   *         parameters via vtkMacros.
   */


  static ttkTrajectoryStatistics *New();
  vtkTypeMacro(ttkTrajectoryStatistics, ttkAlgorithm);

  vtkSetMacro(frameSurface, int);
  vtkGetMacro(frameSurface, int);

  vtkSetMacro(errSurf, double);
  vtkGetMacro(errSurf, double);

  vtkSetMacro(filtreY, double);
  vtkGetMacro(filtreY, double);

  vtkSetMacro(cosCol, double);
  vtkGetMacro(cosCol, double);

  vtkSetMacro(maxRadus, double);
  vtkGetMacro(maxRadus, double);
    
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

  vtkSetMacro(enableFilteringMinVx, int);
  vtkGetMacro(enableFilteringMinVx, int);

  vtkSetMacro(enableFilteringCosY, int);
  vtkGetMacro(enableFilteringCosY, int);

  vtkSetMacro(enableFilteringTimeOrigin, int);
  vtkGetMacro(enableFilteringTimeOrigin, int);

  vtkSetMacro(enableFilteringDuration, int);
  vtkGetMacro(enableFilteringDuration, int);

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

  vtkSetMacro(maxX, int);
  vtkGetMacro(maxX, int);
  vtkSetMacro(maxY, int);
  vtkGetMacro(maxY, int);
  vtkSetMacro(minY, int);
  vtkGetMacro(minY, int);
  vtkSetMacro(minX, int);
  vtkGetMacro(minX, int);



  vtkSetMacro(surfaceMethod, int);
  vtkGetMacro(surfaceMethod, int);
protected:

  ttkTrajectoryStatistics();
  ~ttkTrajectoryStatistics() override = default;

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
