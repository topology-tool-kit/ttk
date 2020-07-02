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
/// The input data array for the first component (u) needs to be specified via
/// the standard VTK call vtkAlgorithm::SetInputArrayToProcess() with the
/// following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The input data array for the second component (v) needs to be specified via
/// the standard VTK call vtkAlgorithm::SetInputArrayToProcess() with the
/// following parameters:
/// \param idx 1 (FIXED: the second array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkTextureMapFromField
///
#pragma once

// VTK includes
#include <vtkSmartPointer.h>

// VTK Module
#include <ttkProjectionFromFieldModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class TTKPROJECTIONFROMFIELD_EXPORT ttkProjectionFromField
  : public ttkAlgorithm {

public:
  static ttkProjectionFromField *New();

  vtkTypeMacro(ttkProjectionFromField, ttkAlgorithm);

  vtkSetMacro(UseTextureCoordinates, bool);
  vtkGetMacro(UseTextureCoordinates, bool);

protected:
  ttkProjectionFromField();

  ~ttkProjectionFromField() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool UseTextureCoordinates{false};
  vtkSmartPointer<vtkPoints> pointSet_;
};
