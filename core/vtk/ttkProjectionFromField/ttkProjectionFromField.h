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
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/builtInExample2/">
///   Builtin example 2</a> \n

#pragma once

// VTK Module
#include <ttkProjectionFromFieldModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class vtkUnstructuredGrid;

class TTKPROJECTIONFROMFIELD_EXPORT ttkProjectionFromField
  : public ttkAlgorithm {

public:
  static ttkProjectionFromField *New();

  vtkTypeMacro(ttkProjectionFromField, ttkAlgorithm);

  vtkSetMacro(UseTextureCoordinates, bool);
  vtkGetMacro(UseTextureCoordinates, bool);

  vtkSetMacro(Use3DCoordinatesArray, bool);
  vtkGetMacro(Use3DCoordinatesArray, bool);

  vtkSetMacro(ProjectPersistenceDiagram, bool);
  vtkGetMacro(ProjectPersistenceDiagram, bool);

protected:
  ttkProjectionFromField();

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * @brief Switch a given Persistence Diagram representation
   */
  int projectPersistenceDiagram(vtkUnstructuredGrid *const inputDiagram,
                                vtkUnstructuredGrid *const outputDiagram);
  /**
   * @brief Generate the spatial embedding of a given Persistence Diagram
   */
  int projectDiagramInsideDomain(vtkUnstructuredGrid *const inputDiagram,
                                 vtkUnstructuredGrid *const outputDiagram);
  /**
   * @brief Generate the 2D embedding of a given Persistence Diagram
   */
  template <typename VTK_TT>
  int projectDiagramIn2D(vtkUnstructuredGrid *const inputDiagram,
                         vtkUnstructuredGrid *const outputDiagram,
                         const VTK_TT *const births,
                         const VTK_TT *const deaths);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool ProjectPersistenceDiagram{false};
  bool UseTextureCoordinates{false};
  bool Use3DCoordinatesArray{false};
};
