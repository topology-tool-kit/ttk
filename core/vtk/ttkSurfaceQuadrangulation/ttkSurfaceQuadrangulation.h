/// \ingroup vtk
/// \class ttkSurfaceQuadrangulation
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
/// \sa ttk::SurfaceQuadrangulation

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
#include <SurfaceQuadrangulation.h>
#include <ttkWrapper.h>

#include <ttkTriangulation.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkSurfaceQuadrangulation
#else
class ttkSurfaceQuadrangulation
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkSurfaceQuadrangulation *New();
  vtkTypeMacro(ttkSurfaceQuadrangulation, vtkDataSetAlgorithm);

  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  vtkGetMacro(DualQuadrangulation, bool);
  vtkSetMacro(DualQuadrangulation, bool);

  // get data from Morse-Smale Complex
  int getCriticalPoints(vtkUnstructuredGrid *input);
  int getSeparatrices(vtkUnstructuredGrid *input);
  int getSegmentation(vtkUnstructuredGrid *input);

  // default copy constructor
  ttkSurfaceQuadrangulation(const ttkSurfaceQuadrangulation &) = delete;
  // default move constructor
  ttkSurfaceQuadrangulation(ttkSurfaceQuadrangulation &&) = delete;
  // default copy assignment operator
  ttkSurfaceQuadrangulation &operator=(const ttkSurfaceQuadrangulation &)
    = delete;
  // default move assignment operator
  ttkSurfaceQuadrangulation &operator=(ttkSurfaceQuadrangulation &&) = delete;

protected:
  ttkSurfaceQuadrangulation();

  ~ttkSurfaceQuadrangulation() override = default;

  TTK_SETUP();

  int FillInputPortInformation(int port, vtkInformation *info) override;

private:
  // output vector of interleaved quadrangles
  std::vector<vtkIdType> outQuadrangles_{};
  // output vector of quadrangle vertices (critical points + ...)
  std::vector<float> outQuadPoints_{};
  // output vector of quadrange vertices ids
  std::vector<ttk::SimplexId> outPointsIds_{};
  // triangulation
  ttk::Triangulation *triangulation_{};
  // if dual quadrangulation
  bool DualQuadrangulation{false};
  // worker object
  ttk::SurfaceQuadrangulation surfaceQuadrangulation_{
    outQuadrangles_, outQuadPoints_, outPointsIds_};
};
