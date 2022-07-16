/// \ingroup vtk
/// \class ttkMetricDistortion
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2022.
///
/// \brief TTK VTK-filter that wraps the ttk::MetricDistortion module.
///
/// This VTK filter uses the ttk::MetricDistortion module to compute distance,
/// area and curvature information about a surface and an optionnal distance
/// matrix (giving the distance between the points of the surface in a metric
/// space).
///
/// \param Input vtkPolyData.
/// \param Input vtkTable (optionnal)
/// \param Output vtkDataSet.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the corresponding standalone program for a usage example:
///   - standalone/MetricDistortion/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MetricDistortion
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkMetricDistortionModule.h>

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
 *   ttkMetricDistortion
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <MetricDistortion.h>

class TTKMETRICDISTORTION_EXPORT ttkMetricDistortion
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::MetricDistortion // and we inherit from the base class
{
private:
  // Filled by the algorithm
  // Area
  std::vector<double> surfaceArea_, metricArea_, ratioArea_;
  // Distance
  std::vector<double> surfaceDistance_, metricDistance_, ratioDistance_;
  std::vector<std::array<double, 3>> surfacePointDistance_,
    metricPointDistance_, ratioPointDistance_;
  // Curvature
  std::vector<double> surfaceCurvature_, metricCurvature_, diffCurvature_;

  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMetricDistortion *New();
  vtkTypeMacro(ttkMetricDistortion, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkMetricDistortion();
  ~ttkMetricDistortion() override;

  /**
   * Specify the input data type of each input port
   * (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * Specify the data object type of each output port
   * (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * Pass VTK data to the base code and convert base code output to VTK
   * (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <class tableDataType>
  int run(vtkInformationVector **inputVector);
};
