/// \ingroup vtk
/// \class ttkQuadrangulationSubdivision
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2019
///
/// \brief TTK VTK-filter for surface quadrangulation.
///
/// The current filter transforms a triangulated surface into a
/// quadrangulated one.
///
/// \param Input0 Input triangular surface (2D) geometry (vtkDataSet)
/// \param Output Quadrangular mesh (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttkScalarFieldCriticalPoints
/// \sa ttkIntegralLines
/// \sa ttkFTMTree
/// \sa ttkIdentifiers
/// \sa ttk::QuadrangulationSubdivision

#pragma once

// VTK includes -- to adapt
#include <vtkCellData.h>
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
#include <QuadrangulationSubdivision.h>
#include <ttkWrapper.h>

#include <ttkTriangulation.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkQuadrangulationSubdivision
#else
class ttkQuadrangulationSubdivision
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkQuadrangulationSubdivision *New();
  vtkTypeMacro(ttkQuadrangulationSubdivision, vtkDataSetAlgorithm);

  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  vtkSetMacro(InputIdentifiersFieldName, std::string);
  vtkGetMacro(InputIdentifiersFieldName, std::string);

  vtkSetMacro(ForceInputIdentifiersField, bool);
  vtkGetMacro(ForceInputIdentifiersField, bool);

  vtkSetMacro(ForceInputOffsetIdentifiersField, bool);
  vtkGetMacro(ForceInputOffsetIdentifiersField, bool);

  vtkSetMacro(SubdivisionLevel, unsigned int);
  vtkGetMacro(SubdivisionLevel, unsigned int);

  vtkSetMacro(RelaxationIterations, unsigned int);
  vtkGetMacro(RelaxationIterations, unsigned int);

  vtkSetMacro(LockInputExtrema, bool);
  vtkGetMacro(LockInputExtrema, bool);

  vtkSetMacro(LockAllInputVertices, bool);
  vtkGetMacro(LockAllInputVertices, bool);

  vtkSetMacro(ReverseProjection, bool);
  vtkGetMacro(ReverseProjection, bool);

  vtkSetMacro(HausdorffLevel, float);
  vtkGetMacro(HausdorffLevel, float);

  vtkSetMacro(ShowResError, bool);
  vtkGetMacro(ShowResError, bool);

  vtkSetMacro(QuadStatistics, bool);
  vtkGetMacro(QuadStatistics, bool);

  // default copy constructor
  ttkQuadrangulationSubdivision(const ttkQuadrangulationSubdivision &) = delete;
  // default move constructor
  ttkQuadrangulationSubdivision(ttkQuadrangulationSubdivision &&) = delete;
  // default copy assignment operator
  ttkQuadrangulationSubdivision &
    operator=(const ttkQuadrangulationSubdivision &)
    = delete;
  // default move assignment operator
  ttkQuadrangulationSubdivision &operator=(ttkQuadrangulationSubdivision &&)
    = delete;

protected:
  ttkQuadrangulationSubdivision();

  ~ttkQuadrangulationSubdivision() override = default;

  TTK_SETUP();

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int getTriangulation(vtkUnstructuredGrid *input);

  int getQuadVertices(vtkUnstructuredGrid *input);

private:
  // user-defined input identifier (SimplexId) scalar field name
  std::string InputIdentifiersFieldName{ttk::VertexScalarFieldName};
  // let the user choose a different identifier scalar field
  bool ForceInputIdentifiersField{false};
  // let the user choose an offset identifier scalar field
  bool ForceInputOffsetIdentifiersField{false};
  // number of subdivisions of the Morse-Smale Complex cells
  unsigned int SubdivisionLevel{1};
  // number of relaxation iterations
  unsigned int RelaxationIterations{10};
  // lock input extrema
  bool LockInputExtrema{false};
  // lock all input vertices
  bool LockAllInputVertices{false};
  // projection method
  bool ReverseProjection{false};
  // Hausdorff warning level
  float HausdorffLevel{200.F};
  // show result despite error
  bool ShowResError{false};
  // display quadrangle statistics
  bool QuadStatistics{false};

  // base worker object
  ttk::QuadrangulationSubdivision baseWorker_{};
};
