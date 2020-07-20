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

  vtkSetMacro(Ucomponent, std::string);
  vtkGetMacro(Ucomponent, std::string);

  vtkSetMacro(Vcomponent, std::string);
  vtkGetMacro(Vcomponent, std::string);

  vtkSetMacro(UcomponentId, int);
  vtkGetMacro(UcomponentId, int);

  vtkSetMacro(VcomponentId, int);
  vtkGetMacro(VcomponentId, int);

  vtkSetMacro(UoffsetId, int);
  vtkGetMacro(UoffsetId, int);

  vtkSetMacro(VoffsetId, int);
  vtkGetMacro(VoffsetId, int);

  vtkGetMacro(ForceInputOffsetScalarField, bool);
  vtkSetMacro(ForceInputOffsetScalarField, bool);

  vtkGetMacro(OffsetFieldU, std::string);
  vtkSetMacro(OffsetFieldU, std::string);

  vtkGetMacro(OffsetFieldV, std::string);
  vtkSetMacro(OffsetFieldV, std::string);

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

private:
  bool ForceInputOffsetScalarField{false};
  bool EdgeIds{false}, VertexScalars{false};
  int UcomponentId{0}, VcomponentId{1}, UoffsetId{-1}, VoffsetId{-1};
  std::string Ucomponent{}, Vcomponent{}, OffsetFieldU{}, OffsetFieldV{};
  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> edgeList_{};
  // for each edge, one skeleton of its triangle fan
  std::vector<std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>>
    edgeFanLinkEdgeList_{};
  // for each edge, the one skeleton of its triangle fan
  std::vector<std::vector<ttk::SimplexId>> edgeFans_{};
  std::vector<std::pair<ttk::SimplexId, char>> jacobiSet_{};
  std::vector<ttk::SimplexId> sosOffsetsU_{}, sosOffsetsV_{};

  template <class dataTypeU, class dataTypeV>
  int baseCall(vtkDataSet *input, vtkDataArray *uField, vtkDataArray *vField);
};
