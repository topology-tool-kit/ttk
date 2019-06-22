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
#ifndef _TTK_REEBSPACE_H
#define _TTK_REEBSPACE_H

// VTK includes -- to adapt
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// ttk code includes
#include <ReebSpace.h>
#include <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkReebSpace
#else
class ttkReebSpace
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkReebSpace *New();

  vtkTypeMacro(ttkReebSpace, vtkDataSetAlgorithm);

  // default ttk setters
  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  vtkSetMacro(Ucomponent, std::string);
  vtkGetMacro(Ucomponent, std::string);

  vtkSetMacro(Vcomponent, std::string);
  vtkGetMacro(Vcomponent, std::string);

  vtkSetMacro(UcomponentId, int);
  vtkGetMacro(UcomponentId, int);

  vtkSetMacro(VcomponentId, int);
  vtkGetMacro(VcomponentId, int);

  vtkGetMacro(ForceInputOffsetScalarField, bool);
  vtkSetMacro(ForceInputOffsetScalarField, bool);

  vtkGetMacro(UseOctreeAcceleration, bool);
  vtkSetMacro(UseOctreeAcceleration, bool);

  vtkGetMacro(OffsetFieldU, std::string);
  vtkSetMacro(OffsetFieldU, std::string);

  vtkGetMacro(OffsetFieldV, std::string);
  vtkSetMacro(OffsetFieldV, std::string);

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

  ~ttkReebSpace();

  virtual int FillOutputPortInformation(int port,
                                        vtkInformation *info) override {

    if(port == 0) {
      // 0-sheets, corners of jacobi set segments
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    } else if(port == 1) {
      // 1-sheets, jacobi sets
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    } else if(port == 2) {
      // 2-sheets, fiber surfaces of jacobi sets
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    }

    return 1;
  }

  TTK_SETUP();

private:
  int UcomponentId, VcomponentId;

  bool ZeroSheetValue, ZeroSheetVertexId, ZeroSheetType, ZeroSheetId;
  bool OneSheetValue, OneSheetVertexId, OneSheetType, OneSheetId,
    OneSheetEdgeId;
  bool TwoSheets, TwoSheetCaseId, TwoSheetValue, TwoSheetParameterization,
    TwoSheetId, TwoSheetEdgeId, TwoSheetTetId, TwoSheetEdgeType;
  bool ThreeSheetVertexNumber, ThreeSheetTetNumber, ThreeSheetExpansion,
    ThreeSheetDomainVolume, ThreeSheetRangeArea, ThreeSheetHyperVolume;

  bool ForceInputOffsetScalarField;
  bool UseOctreeAcceleration;
  int SimplificationCriterion;
  double SimplificationThreshold;
  std::string Ucomponent, Vcomponent, OffsetFieldU, OffsetFieldV;

  vtkDataArray *uComponent_, *vComponent_, *offsetFieldU_, *offsetFieldV_;
  std::vector<ttk::SimplexId> sosOffsetsU_, sosOffsetsV_;

  // core data-structure
  ttk::ReebSpace reebSpace_;

  // template base call
  template <class dataTypeU, class dataTypeV>
  int baseCall(vtkDataSet *input,
               vtkDataArray *uField,
               vtkDataArray *offsetFieldU,
               vtkDataArray *vField,
               vtkDataArray *offsetFieldV);

  template <class dataTypeU, class dataTypeV>
  int preProcess(vtkDataSet *input, vtkDataArray *uField, vtkDataArray *vField);
};

#endif // _TTK_REEBSPACE_H
