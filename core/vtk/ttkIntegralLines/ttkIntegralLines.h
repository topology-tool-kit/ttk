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

#ifdef TTK_ENABLE_MPI
namespace ttk {
  namespace intgl {
    /**
     * @brief Struct used for sorting ghost data during the generation of global
     * ids. A GhostElementsToSort object is created for each segment of integral
     * line that ends or begins with a ghost vertex. minLocalVertexId,
     * seedIdentifier and forkIdentifier are used to identify this segment of
     * integral line on different processes.
     */
    struct GhostElementsToSort {
      ttk::SimplexId ownedGlobalId;
      ttk::SimplexId minLocalVertexId;
      ttk::SimplexId seedIdentifier;
      ttk::SimplexId forkIdentifier{-1};
      ttk::SimplexId globalEdgeId{-1};
      ttk::SimplexId ghostVertexLocalId;
      ttk::SimplexId ghostEdgeLocalId;
    };

    bool operator<(const GhostElementsToSort &left,
                   const GhostElementsToSort &right) {
      if(left.seedIdentifier != right.seedIdentifier) {
        return left.seedIdentifier < right.seedIdentifier;
      }
      if(left.forkIdentifier != right.forkIdentifier) {
        return left.forkIdentifier < right.forkIdentifier;
      }
      return left.minLocalVertexId < right.minLocalVertexId;
    };
  } // namespace intgl
} // namespace ttk
#endif
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

  template <typename triangulationType>
  int getTrajectories(
    vtkDataSet *input,
    triangulationType *triangulation,
    std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
      &trajectories,
    std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
      &forkIdentifiers,
    std::vector<ttk::ArrayLinkedList<std::vector<double>, TABULAR_SIZE>>
      &distancesFromSeed,
    std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
      &seedIdentifiers,
#ifdef TTK_ENABLE_MPI
    std::vector<ttk::SimplexId> &globalVertexId,
    std::vector<ttk::SimplexId> &globalCellId,
#endif
    vtkUnstructuredGrid *output);
#ifdef TTK_ENABLE_MPI
  /**
   * @brief Constructs the global identifiers for vertices and edges.
   *
   * @tparam triangulationType
   * @param globalVertexId vector of global identifiers for vertices
   * @param globalCellId vector of global identifiers for cells
   * @param seedIdentifiers linked list of identifiers of seed, one for each
   * integral line
   * @param forkIdentifiers linked list of identifiers of forks, one for each
   * integral line
   * @param trajectories linked list of identifiers of forks, one for each
   * integral line
   * @param localVertexIdentifiers linked list of identifiers, one for each
   * vertex of an integral lines. It represents the number of the vertex in the
   * integral line.
   * @param triangulation Triangulation
   * @return int
   */
  template <typename triangulationType>
  int getGlobalIdentifiers(
    std::vector<ttk::SimplexId> &globalVertexId,
    std::vector<ttk::SimplexId> &globalCellId,
    std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
      &seedIdentifiers,
    std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
      &forkIdentifiers,
    std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
      &trajectories,
    std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
      &localVertexIdentifiers,
    triangulationType *triangulation);

  /**
   * @brief Sorts ghosts and exchanges their global identifiers between
   * processes.
   *
   * @param unmatchedGhosts Ghosts vertices and edges for which the global
   * identifier is to be determined
   * @param globalVertexId vector of global identifiers for vertices
   * @param globalCellId vector of global identifiers for cells
   * @return int 0 for success
   */
  int exchangeGhosts(
    std::vector<std::vector<ttk::intgl::GhostElementsToSort>> &unmatchedGhosts,
    std::vector<ttk::SimplexId> &globalVertexId,
    std::vector<ttk::SimplexId> &globalCellId);

#endif // TTK_ENABLE_MPI

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
