/// \ingroup vtk
/// \class ttkReebSpace
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2015.
///
/// \brief TTK VTK-filter that efficiently computes the Reeb space of
/// bivariate volumetric data.
///
/// The Reeb space is a useful topological abstraction of bivariate scalar
/// fields for data segmentation purposes. Intuitively, it allows the automatic
/// separation of volumetric regions that project to the same areas in the
/// range. This class also implements a simplification heuristic for progressive
/// coarsening. Used in conjunction with continuous scatterplots, this class
/// enables continuous scattterplot peeling for instance.
///
/// The input data must be provided as two independent point data scalar fields
/// defined on the geometry.
/// \warning Only tetrahedral meshes are supported.
///
/// \param Input Input bivariate data (vtkUnstructuredGrid)
/// \param Output0 Output 0-sheets (vtkUnstructuredGrid)
/// \param Output1 Output 1-sheets (vtkUnstructuredGrid)
/// \param Output2 Output 2-sheets (vtkUnstructuredGrid)
/// \param Output3 Output 3-sheets (vtkUnstructuredGrid)
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
/// "Jacobi Fiber Surfaces for Bivariate Reeb Space Computation" \n
/// Julien Tierny, Hamish Carr \n
/// Proc. of IEEE VIS 2016.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2016.
///
/// \sa vtkContinuousScatterplot
/// \sa vtkJacobiSet
/// \sa vtkFiberSurface
/// \sa vtkRangePolygon
/// \sa ttk::ReebSpace
///

#pragma once

// VTK Module
#include <ttkReebSpaceModule.h>

// ttk code includes
#include <ReebSpace.h>
#include <ttkAlgorithm.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKREEBSPACE_EXPORT ttkReebSpace : public ttkAlgorithm,
                                         protected ttk::ReebSpace {

public:
  static ttkReebSpace *New();
  vtkTypeMacro(ttkReebSpace, ttkAlgorithm);

  vtkGetMacro(ForceInputOffsetScalarField, bool);
  vtkSetMacro(ForceInputOffsetScalarField, bool);

  vtkGetMacro(UseOctreeAcceleration, bool);
  vtkSetMacro(UseOctreeAcceleration, bool);

  // 0-sheet options
  vtkGetMacro(ZeroSheetId, bool);
  vtkSetMacro(ZeroSheetId, bool);

  vtkGetMacro(ZeroSheetType, bool);
  vtkSetMacro(ZeroSheetType, bool);

  vtkGetMacro(ZeroSheetValue, bool);
  vtkSetMacro(ZeroSheetValue, bool);

  vtkGetMacro(ZeroSheetVertexId, bool);
  vtkSetMacro(ZeroSheetVertexId, bool);

  // 1-sheet options
  vtkGetMacro(OneSheetId, bool);
  vtkSetMacro(OneSheetId, bool);

  vtkGetMacro(OneSheetType, bool);
  vtkSetMacro(OneSheetType, bool);

  vtkGetMacro(OneSheetValue, bool);
  vtkSetMacro(OneSheetValue, bool);

  vtkGetMacro(OneSheetVertexId, bool);
  vtkSetMacro(OneSheetVertexId, bool);

  vtkGetMacro(OneSheetEdgeId, bool);
  vtkSetMacro(OneSheetEdgeId, bool);

  // 2-sheet options
  vtkGetMacro(TwoSheets, bool);
  vtkSetMacro(TwoSheets, bool);

  vtkGetMacro(TwoSheetCaseId, bool);
  vtkSetMacro(TwoSheetCaseId, bool);

  vtkGetMacro(TwoSheetValue, bool);
  vtkSetMacro(TwoSheetValue, bool);

  vtkGetMacro(TwoSheetParameterization, bool);
  vtkSetMacro(TwoSheetParameterization, bool);

  vtkGetMacro(TwoSheetId, bool);
  vtkSetMacro(TwoSheetId, bool);

  vtkGetMacro(TwoSheetEdgeId, bool);
  vtkSetMacro(TwoSheetEdgeId, bool);

  vtkGetMacro(TwoSheetTetId, bool);
  vtkSetMacro(TwoSheetTetId, bool);

  vtkGetMacro(TwoSheetEdgeType, bool);
  vtkSetMacro(TwoSheetEdgeType, bool);

  vtkGetMacro(ThreeSheetTetNumber, bool);
  vtkSetMacro(ThreeSheetTetNumber, bool);

  vtkGetMacro(ThreeSheetVertexNumber, bool);
  vtkSetMacro(ThreeSheetVertexNumber, bool);

  vtkGetMacro(ThreeSheetExpansion, bool);
  vtkSetMacro(ThreeSheetExpansion, bool);

  vtkGetMacro(ThreeSheetDomainVolume, bool);
  vtkSetMacro(ThreeSheetDomainVolume, bool);

  vtkGetMacro(ThreeSheetRangeArea, bool);
  vtkSetMacro(ThreeSheetRangeArea, bool);

  vtkGetMacro(ThreeSheetHyperVolume, bool);
  vtkSetMacro(ThreeSheetHyperVolume, bool);

  vtkGetMacro(SimplificationThreshold, double);
  vtkSetMacro(SimplificationThreshold, double);

  vtkGetMacro(SimplificationCriterion, int);
  vtkSetMacro(SimplificationCriterion, int);

protected:
  ttkReebSpace();

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
  bool ZeroSheetValue{true}, ZeroSheetVertexId{true}, ZeroSheetType{true},
    ZeroSheetId{true};
  bool OneSheetValue{true}, OneSheetVertexId{true}, OneSheetType{true},
    OneSheetId{true}, OneSheetEdgeId{true};
  bool TwoSheets{true}, TwoSheetCaseId{true}, TwoSheetValue{true},
    TwoSheetParameterization{true}, TwoSheetId{true}, TwoSheetEdgeId{true},
    TwoSheetTetId{true}, TwoSheetEdgeType{true};
  bool ThreeSheetVertexNumber{true}, ThreeSheetTetNumber{true},
    ThreeSheetExpansion{true}, ThreeSheetDomainVolume{true},
    ThreeSheetRangeArea{true}, ThreeSheetHyperVolume{true};

  bool ForceInputOffsetScalarField{false};
  bool UseOctreeAcceleration{true};
  int SimplificationCriterion{1};
  double SimplificationThreshold{0.0};
};
