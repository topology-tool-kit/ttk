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

#include <ttkAlgorithm.h>
#include <ttkPointMergerModule.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKPOINTMERGER_EXPORT ttkPointMerger : public ttkAlgorithm {
private:
  bool BoundaryOnly{true};
  double DistanceThreshold{0.001};

public:
  vtkSetMacro(BoundaryOnly, bool);
  vtkGetMacro(BoundaryOnly, bool);

  vtkSetMacro(DistanceThreshold, double);
  vtkGetMacro(DistanceThreshold, double);

  static ttkPointMerger *New();
  vtkTypeMacro(ttkPointMerger, ttkAlgorithm);

protected:
  ttkPointMerger();
  ~ttkPointMerger();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
