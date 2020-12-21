/// \ingroup vtk
/// \class ttkPointMerger
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2018.
///
/// \brief TTK VTK-filter for point merging.
///
/// This filter merges the points of a mesh whose distance is lower than a user
/// defined threshold.
///
/// \param Input Input data set (vtkDataSet)
/// \param Output Output data set (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
#pragma once

// VTK Module
#include <ttkPointMergerModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class TTKPOINTMERGER_EXPORT ttkPointMerger : public ttkAlgorithm {
public:
  static ttkPointMerger *New();
  vtkTypeMacro(ttkPointMerger, ttkAlgorithm);

  vtkSetMacro(BoundaryOnly, bool);
  vtkGetMacro(BoundaryOnly, bool);

  vtkSetMacro(DistanceThreshold, double);
  vtkGetMacro(DistanceThreshold, double);

protected:
  ttkPointMerger();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool BoundaryOnly{true};
  double DistanceThreshold{0.001};
};
