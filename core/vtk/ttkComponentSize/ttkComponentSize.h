/// \ingroup vtk
/// \class ttkComponentSize
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date September 2015.
///
/// \brief TTK VTK-filter that computes the connected components of a data-set,
/// and that computes their size (number of vertices, number of cells, etc.).
///
/// This filter computes the connected component of a point-set data-set and
/// computes their size (number of vertices, number of cells, etc). The size
/// information is attached on the output, either as point or cell data.
/// The identifier of each connected component is also attached to the geometry
/// as point or cell data.
///
/// This filter is useful when used in conjunction with some thresholding, to
/// only display the largest connected components of a data-set.
///
/// \param Input Input data-set (vtkPointSet)
/// \param Output Output data-set (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///

#pragma once

// VTK Module
#include <ttkComponentSizeModule.h>

// TTK Include
#include <ttkAlgorithm.h>

class TTKCOMPONENTSIZE_EXPORT ttkComponentSize : public ttkAlgorithm {

public:
  static ttkComponentSize *New();
  vtkTypeMacro(ttkComponentSize, ttkAlgorithm);

protected:
  ttkComponentSize();
  ~ttkComponentSize();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
