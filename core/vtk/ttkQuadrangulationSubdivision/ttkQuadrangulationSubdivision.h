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

  // fill triangulation_ array
  int getTriangulation(vtkDataSet *input);

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
  unsigned int SubdivisionLevel{3};
  // number of relaxation iterations
  unsigned int RelaxationIterations{100};
  // lock input extrema
  bool LockInputExtrema{true};
  // lock all input vertices
  bool LockAllInputVertices{true};
  // projection method
  bool ReverseProjection{false};

  // base worker object
  ttk::QuadrangulationSubdivision baseWorker_{};
  // output vector of interleaved quadrangles
  std::vector<vtkIdType> outQuadrangles_{};
  // output vector of quadrangle vertices
  std::vector<float> outVertices_{};
  // output vector of quadrangle vertex valences
  std::vector<ttk::SimplexId> outVertexValences_{};
  // output vector of vertex (debug) info
  // 0 - input (critical) point
  // 1 - edge middle
  // 2 - quadrangle barycenter
  std::vector<ttk::SimplexId> outVertexType_{};
  // subdivision level that generated the points
  std::vector<ttk::SimplexId> outSubdvisionLevel_{};
  // number of triangles checked during last projection
  std::vector<ttk::SimplexId> trianglesChecked_{};
  // success of last projection
  // 0 - not projected (critical point)
  // 1 - projection alongside quadrangle normal
  // 2 - projection alongside triangle normal
  // 3 - failed projection
  std::vector<ttk::SimplexId> projSucceeded_{};
};
