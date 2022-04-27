/// \ingroup vtk
/// \class ttkTrackingFromFields
/// \author Maxime Soler <soler.maxime@total.com>
/// \date August 2018.
///
/// \brief TTK VTK-filter that takes an input time-varying data set (represented
/// by a list of scalar fields) and which computes a tracking mesh.
///
/// \param Input Input time-dependent scalar field, either 2D or 3D, regular
/// grid or triangulation (vtkDataSet); time steps are obtained by
/// GetPointData()->GetArray(i) in increasing time order.
/// \param Output Output persistence diagram (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Lifted Wasserstein Matcher for Fast and Robust Topology Tracking" \n
/// Maxime Soler, Melanie Plainchault, Bruno Conche, Julien Tierny \n
/// Proc. of IEEE Symposium on Large Data Analysis and Visualization, 2018
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/timeTracking/">Time
///   tracking example</a>
///

#pragma once

#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <TrackingFromFields.h>
#include <ttkAlgorithm.h>
#include <ttkTrackingFromFieldsModule.h>

#include <algorithm>
#include <string>

class TTKTRACKINGFROMFIELDS_EXPORT ttkTrackingFromFields
  : public ttkAlgorithm,
    protected ttk::TrackingFromFields {

public:
  static ttkTrackingFromFields *New();

  vtkTypeMacro(ttkTrackingFromFields, ttkAlgorithm);

  /// @brief Temporal sampling (take every N timestep).
  /// @{
  vtkSetMacro(Sampling, int);
  vtkGetMacro(Sampling, int);
  /// @}

  /// @brief First timestep.
  /// @{
  vtkSetMacro(StartTimestep, int);
  vtkGetMacro(StartTimestep, int);
  /// @}

  /// @brief Last timestep (-1 to use the last timestep available).
  /// @{
  vtkSetMacro(EndTimestep, int);
  vtkGetMacro(EndTimestep, int);
  /// @}

  /// @brief Discard pairs below this threshold (percentage of the function
  /// span).
  /// @{
  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);
  /// @}

  /// @brief Importance weight for the X component of the extremum.
  /// @{
  vtkSetMacro(PX, double);
  vtkGetMacro(PX, double);
  /// @}

  /// @brief Importance weight for the Y component of the extremum.
  /// @{
  vtkSetMacro(PY, double);
  vtkGetMacro(PY, double);
  /// @}

  /// @brief Importance weight for the Z component of the extremum.
  /// @{
  vtkSetMacro(PZ, double);
  vtkGetMacro(PZ, double);
  /// @}

  /// @brief Importance weight for extrema.
  /// @{
  vtkSetMacro(PE, double);
  vtkGetMacro(PE, double);
  /// @}

  /// @brief Importance weight for saddles.
  /// @{
  vtkSetMacro(PS, double);
  vtkGetMacro(PS, double);
  /// @}

  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);

  /// @brief Value of the parameter p for the Wp (p-th Wasserstein) distance
  /// computation (type "inf" for the Bottleneck distance).
  /// @{
  vtkSetMacro(WassersteinMetric, const std::string &);
  vtkGetMacro(WassersteinMetric, std::string);
  /// @}

  vtkSetMacro(DistanceAlgorithm, const std::string &);
  vtkGetMacro(DistanceAlgorithm, std::string);

  /// @brief Method for computing matchings.
  ///
  /// 0: sparse Munkres (Wasserstein), Gabow-Tarjan (Bottleneck)
  /// @{
  vtkSetMacro(PVAlgorithm, int);
  vtkGetMacro(PVAlgorithm, int);
  /// @}

  /// @brief For the translation of the second set of critical points even the
  /// persistence diagrams are embedded in the original domain. This is useful
  /// to visualize the matching between the diagrams of two 2D scalar fields.
  /// @{
  vtkSetMacro(UseGeometricSpacing, bool);
  vtkGetMacro(UseGeometricSpacing, bool);
  /// @}

  /// @brief Translation on the Z axis between the output representations of the
  /// persistence diagrams.
  /// @{
  vtkSetMacro(Spacing, double);
  vtkGetMacro(Spacing, double);
  /// @}

  /// @brief Do post-processing.
  /// @{
  vtkSetMacro(DoPostProc, bool);
  vtkGetMacro(DoPostProc, bool);
  /// @}

  /// @brief Threshold for merging/splitting trajectories in connected
  /// components array.
  /// @{
  vtkSetMacro(PostProcThresh, double);
  vtkGetMacro(PostProcThresh, double);
  /// @}

protected:
  ttkTrackingFromFields();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
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
  bool DoPostProc{false};
  double PostProcThresh{0.0};
  double Spacing{1.0};
  double Alpha{1.0};
  std::string DistanceAlgorithm{"ttk"};
  int PVAlgorithm{-1};
  std::string WassersteinMetric{"2"};

  template <class dataType, class triangulationType>
  int trackWithPersistenceMatching(vtkUnstructuredGrid *output,
                                   unsigned long fieldNumber,
                                   const triangulationType *triangulation);
};
