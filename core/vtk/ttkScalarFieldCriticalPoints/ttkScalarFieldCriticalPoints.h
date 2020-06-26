/// \ingroup vtk
/// \class ttkScalarFieldCriticalPoints
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date June 2015.
///
/// \brief TTK VTK-filter for the computation of critical points in PL
/// scalar fields defined on PL manifolds.
///
/// This filter computes the list of critical points of the input scalar field
/// and classify them according to their type.
///
/// \param Input Input PL scalar field (vtkDataSet)
/// \param Output Output critical points (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// The name of the input scalar field to consider can be specified with the
/// standard VTK call SetInputArrayToProcess(), with the following parameters:
/// \param port 0
/// \param connection 0
/// \param fieldAssociation 0 (point data)
/// \param arrayName (const char* string representing the name of the VTK array)
///
/// This module respects the following convention regarding the order of the
/// input arrays to process (SetInputArrayToProcess()):
/// \param idx 0: input scalar field
/// \param idx 1 (optional): input offset scalar field
///
/// See the corresponding standalone program for a usage example:
///   - standalone/ScalarFieldCriticalPoints/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Critical points and curvature for embedded polyhedral surfaces" \n
/// Thomas Banchoff \n
/// American Mathematical Monthly, 1970.
///
/// \sa ttk::ScalarFieldCriticalPoints
///
#pragma once

// VTK Module
#include <ttkScalarFieldCriticalPointsModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// ttk baseCode includes
#include <ScalarFieldCriticalPoints.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKSCALARFIELDCRITICALPOINTS_EXPORT ttkScalarFieldCriticalPoints
  : public ttkAlgorithm,
    protected ttk::ScalarFieldCriticalPoints {

public:
  static ttkScalarFieldCriticalPoints *New();

  vtkTypeMacro(ttkScalarFieldCriticalPoints, ttkAlgorithm);

  vtkGetMacro(VertexBoundary, bool);
  vtkSetMacro(VertexBoundary, bool);

  vtkGetMacro(VertexIds, bool);
  vtkSetMacro(VertexIds, bool);

  vtkGetMacro(VertexScalars, bool);
  vtkSetMacro(VertexScalars, bool);

  vtkGetMacro(ForceInputOffsetScalarField, bool);
  vtkSetMacro(ForceInputOffsetScalarField, bool);

protected:
  ttkScalarFieldCriticalPoints();

  ~ttkScalarFieldCriticalPoints() override;
  
  int FillInputPortInformation(int port, vtkInformation *info) override;
  
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  
  int RequestData(vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector) override;

private:
  
  bool ForceInputOffsetScalarField{false};
  bool VertexIds{true}, VertexScalars{true}, VertexBoundary{true};
  
  std::vector<std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>>
    vertexLinkEdgeList_;
  std::vector<std::pair<ttk::SimplexId, char>> criticalPoints_;
  std::vector<ttk::SimplexId> sosOffsets_;
};
