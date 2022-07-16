/// \ingroup vtk
/// \class ttkJacobiSet
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date May 2015.
///
/// \brief TTK VTK-filter that computes the Jacobi set of a bivariate
/// volumetric data-set.
///
/// Given a bivariate scalar field defined on a PL 3-manifold, this filter
/// produces the list of Jacobi edges (each entry is a pair given by the edge
/// identifier and the Jacobi edge type).
///
/// The input bivariate data must be provided as two independent scalar fields
/// attached as point data to the input geometry.
///
/// \param Input Input bivariate volumetric data (vtkDataSet)
/// \param Output Output Jacobi set (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// The input data arrays needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 for the U Component, 1 for the V Component
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The optional offset arrays can be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 2 for the U Offset Field, 3 for the V Offset Field
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the offset array)
/// \note: To use this optional array, `ForceInputOffsetScalarField` needs to be
/// enabled with the setter `setForceInputOffsetScalarField()'.
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Jacobi sets of multiple Morse functions" \n
/// Herbert Edelsbrunner, John Harer \n
/// Foundations of Computational Mathematics. Cambridge University Press, 2002.
///
/// \sa ttk::JacobiSet
/// \sa vtkReebSpace
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/builtInExample2/">
///   Builtin example 2</a> \n

#pragma once

// VTK Module
#include <ttkJacobiSetModule.h>

// ttk code includes
#include <JacobiSet.h>
#include <ttkAlgorithm.h>

class TTKJACOBISET_EXPORT ttkJacobiSet : public ttkAlgorithm,
                                         protected ttk::JacobiSet {
public:
  static ttkJacobiSet *New();
  vtkTypeMacro(ttkJacobiSet, ttkAlgorithm);

  vtkGetMacro(ForceInputOffsetScalarField, bool);
  vtkSetMacro(ForceInputOffsetScalarField, bool);

  vtkSetMacro(EdgeIds, bool);
  vtkGetMacro(EdgeIds, bool);

  vtkSetMacro(VertexScalars, bool);
  vtkGetMacro(VertexScalars, bool);

protected:
  ttkJacobiSet();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <class dataTypeU, class dataTypeV>
  int dispatch(const dataTypeU *const uField,
               const dataTypeV *const vField,
               ttk::Triangulation *const triangulation);

private:
  bool ForceInputOffsetScalarField{false};
  bool EdgeIds{false}, VertexScalars{false};
  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> edgeList_{};
  // for each edge, one skeleton of its triangle fan
  std::vector<std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>>
    edgeFanLinkEdgeList_{};
  // for each edge, the one skeleton of its triangle fan
  std::vector<std::vector<ttk::SimplexId>> edgeFans_{};
  std::vector<std::pair<ttk::SimplexId, char>> jacobiSet_{};
  std::vector<char> isPareto_{};
};
