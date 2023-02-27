/// \ingroup vtk
/// \class ttkIntegralLines
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date March 2016
/// \date MPI implementation: December 2022
///
/// \brief TTK VTK-filter for the computation of edge-based integral lines of
/// the gradient of an input scalar field.
///
/// The filter takes on its input a scalar field attached as point data to an
/// input geometry (either 2D or 3D, either regular grids or triangulations)
/// and computes the forward or backward integral lines along the edges of the
/// input mesh, given a list of input sources.
/// The sources are specified with a vtkPointSet on which is attached as point
/// data a scalar field that represent the vertex identifiers of the sources in
/// the input geometry.
///
/// \param Input0 Input scalar field, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet)
/// \param Input1 Input sources (vtkPointSet)
/// \param Output Output integral lines (vtkUnstructuredGrid)
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The optional offset array can be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 1 (FIXED: the second array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the offset array)
/// \note: To use this optional array, `ForceInputOffsetScalarField` needs to be
/// enabled with the setter `setForceInputOffsetScalarField()'.
///
/// The vertex identifier array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 2 (FIXED: the third array the algorithm requires)
/// \param port 1 (FIXED: second port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the vertex identifier array)
/// \note: To use this optional array, `ForceInputVertexScalarField` needs to be
/// enabled with the setter `setForceInputVertexScalarField()'.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::IntegralLines
/// \sa vtkIdentifiers

#pragma once

// VTK Module
#include <ttkIntegralLinesModule.h>

// ttk code includes
#include <IntegralLines.h>
#include <ttkAlgorithm.h>

class vtkUnstructuredGrid;

class TTKINTEGRALLINES_EXPORT ttkIntegralLines : public ttkAlgorithm,
                                                 protected ttk::IntegralLines {

public:
  static ttkIntegralLines *New();

  vtkTypeMacro(ttkIntegralLines, ttkAlgorithm);

  vtkGetMacro(Direction, int);
  vtkSetMacro(Direction, int);

  vtkSetMacro(ForceInputVertexScalarField, bool);
  vtkGetMacro(ForceInputVertexScalarField, bool);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  vtkSetMacro(EnableForking, bool);
  vtkGetMacro(EnableForking, bool);

  /**
   * This method converts the output data in VTK data types.
   * It constructs the new unstructured grid and associates scalar data
   * to its points and cells.
   */
  template <typename triangulationType>
  int getTrajectories(
    vtkDataSet *input,
    const triangulationType *triangulation,
    const std::vector<ttk::ArrayLinkedList<ttk::intgl::IntegralLine,
                                           INTEGRAL_LINE_TABULAR_SIZE>>
      &integralLines,
#ifdef TTK_ENABLE_MPI
    const std::vector<ttk::SimplexId> &globalVertexId,
    const std::vector<ttk::SimplexId> &globalCellId,
#endif
    vtkUnstructuredGrid *output);

protected:
  ttkIntegralLines();
  ~ttkIntegralLines() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int Direction{0};
  bool ForceInputVertexScalarField{false};
  bool ForceInputOffsetScalarField{false};
};
