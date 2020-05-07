/// \ingroup vtk
/// \class ttkProjectionFromField
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK VTK-filter which projects a data-set to 2D given two point-data
/// scalar fields to be used as 2D coordinates.
///
/// \param Input Input data-set, with at least two point data scalar fields or
/// texture coordinates (vtkPointSet)
/// \param Output Output projected data-set (vtkPointSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkTextureMapFromField
///
#pragma once

#include <ttkAlgorithm.h>
#include <ttkProjectionFromFieldModule.h>

class TTKPROJECTIONFROMFIELD_EXPORT ttkProjectionFromField
  : public ttkAlgorithm {
private:
  bool UseTextureCoordinates;

public:
  vtkSetMacro(UseTextureCoordinates, bool);
  vtkGetMacro(UseTextureCoordinates, bool);

  static ttkProjectionFromField *New();
  vtkTypeMacro(ttkProjectionFromField, ttkAlgorithm);

protected:
  ttkProjectionFromField();
  ~ttkProjectionFromField();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};