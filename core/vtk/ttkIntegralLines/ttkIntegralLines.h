/// \ingroup vtk
/// \class ttkIntegralLines
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date March 2016
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
/// This module respects the following convention regarding the order of the
/// input arrays to process (SetInputArrayToProcess()):
/// \param idx 0: input data array to average.
/// \param idx 1: offset scalar field.
/// \param idx 2: vertex identifiers.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::IntegralLines
/// \sa vtkIdentifiers

#ifndef _TTK_DISCRETESTREAMLINE_H
#define _TTK_DISCRETESTREAMLINE_H

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

  int getTrajectories(vtkDataSet *input,
                      ttk::Triangulation *triangulation,
                      std::vector<std::vector<ttk::SimplexId>> &trajectories,
                      vtkUnstructuredGrid *output);

  template <typename VTK_TT, typename TTK_TT>
  int dispatch(int, const TTK_TT*);

protected:
  ttkIntegralLines();
  ~ttkIntegralLines() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
private:
  int Direction {0};
  bool ForceInputVertexScalarField {false};
  bool ForceInputOffsetScalarField {false};
};

#endif // _TTK_DISCRETESTREAMLINE_H
