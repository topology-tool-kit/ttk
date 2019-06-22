/// \ingroup vtk
/// \class ttkTopologicalSimplification
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date February 2016
///
/// \brief TTK VTK-filter for the topological simplification of scalar
/// data.
///
/// Given an input scalar field and a list of critical points to remove, this
/// filter minimally edits the scalar field such that the listed critical points
/// disappear. This procedure is useful to speedup subsequent topological data
/// analysis when outlier critical points can be easily identified. It is
/// also useful for data simplification.
///
/// The list of critical points to remove must be associated with a point data
/// scalar field that represent the vertex global identifiers in the input
/// geometry.
///
/// Note that this filter will also produce an output vertex offset scalar field
/// that can be used for further topological data analysis tasks to disambiguate
/// vertices on flat plateaus. For instance, this output vertex offset field
/// can specified to the ttkFTMTree, vtkIntegralLines, or
/// vtkScalarFieldCriticalPoints filters.
///
/// Also, this filter can be given a specific input vertex offset.
///
/// \param Input0 Input scalar field, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet)
/// \param Input1 List of critical point constraints (vtkPointSet)
/// \param Output Output simplified scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Generalized Topological Simplification of Scalar Fields on Surfaces" \n
/// Julien Tierny, Valerio Pascucci \n
/// Proc. of IEEE VIS 2012.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2012.
///
/// \sa ttkScalarFieldCriticalPoints
/// \sa ttkIntegralLines
/// \sa ttkFTMTree
/// \sa ttkIdentifiers
/// \sa ttk::TopologicalSimplification
#ifndef _TTK_TOPOLOGICALSIMPLIFICATION_H
#define _TTK_TOPOLOGICALSIMPLIFICATION_H

// VTK includes -- to adapt
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>

// ttk code includes
#include <TopologicalSimplification.h>
#include <ttkWrapper.h>

#include <ttkTriangulation.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTopologicalSimplification
#else
class ttkTopologicalSimplification
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkTopologicalSimplification *New();
  vtkTypeMacro(ttkTopologicalSimplification, vtkDataSetAlgorithm);

  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(ForceInputOffsetScalarField, int);
  vtkGetMacro(ForceInputOffsetScalarField, int);

  vtkSetMacro(ConsiderIdentifierAsBlackList, int);
  vtkGetMacro(ConsiderIdentifierAsBlackList, int);

  vtkSetMacro(AddPerturbation, int);
  vtkGetMacro(AddPerturbation, int);

  vtkSetMacro(InputOffsetScalarFieldName, std::string);
  vtkGetMacro(InputOffsetScalarFieldName, std::string);

  vtkSetMacro(OutputOffsetScalarFieldName, std::string);
  vtkGetMacro(OutputOffsetScalarFieldName, std::string);

  vtkSetMacro(ForceInputVertexScalarField, int);
  vtkGetMacro(ForceInputVertexScalarField, int);

  vtkSetMacro(InputVertexScalarFieldName, std::string);
  vtkGetMacro(InputVertexScalarFieldName, std::string);

  int getTriangulation(vtkDataSet *input);
  int getScalars(vtkDataSet *input);
  int getIdentifiers(vtkPointSet *input);
  int getOffsets(vtkDataSet *input);

protected:
  ttkTopologicalSimplification();

  ~ttkTopologicalSimplification();

  TTK_SETUP();

  int FillInputPortInformation(int port, vtkInformation *info) override;

private:
  int ScalarFieldId;
  int OffsetFieldId;
  std::string ScalarField;
  std::string InputOffsetScalarFieldName;
  std::string OutputOffsetScalarFieldName;
  bool ForceInputVertexScalarField;
  std::string InputVertexScalarFieldName;
  bool ForceInputOffsetScalarField;
  bool ConsiderIdentifierAsBlackList;
  bool AddPerturbation;
  bool hasUpdatedMesh_;

  ttk::TopologicalSimplification topologicalSimplification_;
  ttk::Triangulation *triangulation_;
  vtkDataArray *identifiers_;
  vtkDataArray *inputScalars_;
  vtkDataArray *offsets_;
  vtkDataArray *inputOffsets_;
};

#endif // _TTK_TOPOLOGICALSIMPLIFICATION_H
