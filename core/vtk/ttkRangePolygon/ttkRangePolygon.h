/// \ingroup vtk
/// \class ttkRangePolygon
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date May 2016.
///
/// \brief TTK VTK-filter which produces a valid range polygon for fiber
/// surface extraction.
///
/// Given an input 2D selection, this filter produces a polygon to be used as
/// an input to vtkFiberSurface. Typically, users generate a 2D selection from
/// the continuous scatterplot and this filter translates this selection into
/// a valid range polygon.
///
/// The user can either select:
///   -# Cells on the continuous scatterplot ("Select Cells with Polygon" in
/// ParaView)
/// \warning If used from ParaView, users may need to zoom in sufficiently so
/// that the "Extract Selection" indeed captures all of the user brushed cells.
/// \warning In this case, the generated polygon may count A LOT of edges,
/// which will seriously increase run-times. The next alternative is the
/// default recommendation.
///   -# Points on the continuous scatterplot ("Interactive Select Points on"
/// in ParaView)
/// \warning This feature will only work properly with the TTK-branded ParaView
/// (ParaView sources need to be patched with TTK fixes, see the documentation)
///
/// -
///
/// \param Input Input 2D selection, typically "Extract Selection" in ParaView
/// (vtkUnstructuredGrid)
/// \param Output Range polygon to be used with vtkFiberSurface
/// (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Fast and Exact Fiber Surface Extraction for Tetrahedral Meshes" \n
/// Pavol Klacansky, Julien Tierny, Hamish Carr, Zhao Geng \n
/// IEEE Transactions on Visualization and Computer Graphics, 2016.
///
/// \sa vtkContinuousScatterplot
/// \sa vtkFiberSurface
/// \sa ttk::FiberSurface
/// \sa vtkReebSpace
///
#ifndef _TTK_RANGEPOLYGON_H
#define _TTK_RANGEPOLYGON_H

// VTK includes -- to adapt
#include <vtkCleanPolyData.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// ttk code includes
#include <ScalarFieldSmoother.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkRangePolygon
#else
class ttkRangePolygon
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkRangePolygon *New();

  vtkTypeMacro(ttkRangePolygon, vtkDataSetAlgorithm);

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

  vtkGetMacro(ClosedLoop, bool);
  vtkSetMacro(ClosedLoop, bool);

  vtkGetMacro(NumberOfIterations, int);
  vtkSetMacro(NumberOfIterations, int);

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }

protected:
  ttkRangePolygon();

  ~ttkRangePolygon();

  TTK_SETUP();

private:
  bool ClosedLoop;
  int NumberOfIterations;

  int processPoints(vtkUnstructuredGrid *input, vtkUnstructuredGrid *output);

  int processTriangles(vtkUnstructuredGrid *input, vtkUnstructuredGrid *output);
};

#endif // _TTK_RANGEPOLYGON_H
