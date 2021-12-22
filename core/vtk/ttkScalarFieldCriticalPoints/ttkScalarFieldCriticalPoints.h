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
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
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
///  Progressive Approach used by default
///
/// \b Related \b publication \n
/// "A Progressive Approach to Scalar Field Topology" \n
/// Jules Vidal, Pierre Guillou, Julien Tierny\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// \sa ttk::ScalarFieldCriticalPoints
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/dragon/">Dragon
///   example</a>
///   - <a
///   href="https://topology-tool-kit.github.io/examples/uncertainStartingVortex/">
///   Uncertain Starting Vortex example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/interactionSites/">
///   Interaction sites</a> \n
///
#pragma once

// VTK Module
#include <ttkScalarFieldCriticalPointsModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

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

  ttkSetEnumMacro(BackEnd, BACKEND);
  vtkGetEnumMacro(BackEnd, BACKEND);

  vtkGetMacro(StartingResolutionLevel, int);
  vtkSetMacro(StartingResolutionLevel, int);

  vtkGetMacro(StoppingResolutionLevel, int);
  vtkSetMacro(StoppingResolutionLevel, int);

  vtkGetMacro(IsResumable, bool);
  vtkSetMacro(IsResumable, bool);

  vtkGetMacro(TimeLimit, double);
  vtkSetMacro(TimeLimit, double);

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
};
