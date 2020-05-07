/// \ingroup vtk
/// \class ttkFiber
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date May 2016.
///
/// \brief TTK VTK-filter for fiber computation on bivariate volumetric data.
///
/// Given a point in the range, this filter computes its fiber (i.e. pre-image)
/// on bivariate volumetric data. The bivariate input data must be provided as
/// two independent scalar fields attached as point data to the input geometry.
///
/// \param Input Input bivariate volumetric data-set, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output Fiber (vtkPolyData)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Fast and Exact Fiber Surface Extraction for Tetrahedral Meshes" \n
/// Pavol Klacansky, Julien Tierny, Hamish Carr, Zhao Geng \n
/// IEEE Transactions on Visualization and Computer Graphics, 2016.
///
/// \sa ttkFiberSurface
///
#pragma once

#include <ttkAlgorithm.h>
#include <ttkFiberModule.h>

class TTKFIBER_EXPORT ttkFiber : public ttkAlgorithm {

private:
  double UValue{0};
  double VValue{0};

public:
  vtkGetMacro(UValue, double);
  vtkSetMacro(UValue, double);

  vtkGetMacro(VValue, double);
  vtkSetMacro(VValue, double);

  vtkTypeMacro(ttkFiber, ttkAlgorithm);
  static ttkFiber *New();

protected:
  ttkFiber();
  ~ttkFiber();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};