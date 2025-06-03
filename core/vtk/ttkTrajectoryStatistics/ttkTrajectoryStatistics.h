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
  bool addThresh{true};
  bool gradThresh{false};

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

  vtkSetMacro(addThresh, bool);
  vtkGetMacro(addThresh, bool);

  vtkSetMacro(gradThresh, bool);
  vtkGetMacro(gradThresh, bool);


protected:

  ttkTrajectoryStatistics();
  ~ttkTrajectoryStatistics() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;


  int FillOutputPortInformation(int port, vtkInformation *info) override;


  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
