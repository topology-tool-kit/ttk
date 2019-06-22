/// \ingroup vtk
/// \class ttkFiberSurface
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2015.
///
/// \brief TTK VTK-filter that computes fiber surfaces.
///
/// Fiber surfaces are defined as the pre-images of curves drawn in the range
/// of bivariate volumetric functions, typically on top of the continuous
/// scatterplot. Fiber surfaces generalize the segmentation features of
/// isosurfaces to bivariate data.
/// This filter implements an exact, parallel and fast algorithm for fiber
/// surface computation on (explicit or implicit) tetrahedral meshes.
///
/// The input bivariate data must be provided as two independent scalar fields
/// attached as point data to the input geometry. The input range polygon must
/// be provided as a vtkUnstructuredGrid with the actual 2D locations of the
/// vertices also provided as two independent scalar fields attached as point
/// data to the geometry. See vtkRangePolygon to create such an input polygon
/// from sparse user inputs.
///
/// \param Input0 Input bivariate volumetric data, either regular grid or
/// triangulation (vtkDataSet)
/// \param Input1 Input range polygon (vtkUnstructuredGrid)
/// \param Output Output fiber surface (vtkPolyData)
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
/// \sa ttk::FiberSurface
/// \sa vtkRangePolygon

#ifndef _TTK_FIBERSURFACE_H
#define _TTK_FIBERSURFACE_H

// VTK includes -- to adapt
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// ttk code includes
#include <FiberSurface.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkFiberSurface
#else
class ttkFiberSurface
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkFiberSurface *New();

  vtkTypeMacro(ttkFiberSurface, vtkDataSetAlgorithm);

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

  vtkGetMacro(DataUcomponent, std::string);
  vtkSetMacro(DataUcomponent, std::string);

  vtkGetMacro(DataVcomponent, std::string);
  vtkSetMacro(DataVcomponent, std::string);

  vtkGetMacro(PolygonUcomponent, std::string);
  vtkSetMacro(PolygonUcomponent, std::string);

  vtkGetMacro(PolygonVcomponent, std::string);
  vtkSetMacro(PolygonVcomponent, std::string);

  vtkGetMacro(RangeCoordinates, bool);
  vtkSetMacro(RangeCoordinates, bool);

  vtkGetMacro(EdgeParameterization, bool);
  vtkSetMacro(EdgeParameterization, bool);

  vtkGetMacro(EdgeIds, bool);
  vtkSetMacro(EdgeIds, bool);

  vtkGetMacro(TetIds, bool);
  vtkSetMacro(TetIds, bool);

  vtkGetMacro(CaseIds, bool);
  vtkSetMacro(CaseIds, bool);

  vtkGetMacro(PointMerge, bool);
  vtkSetMacro(PointMerge, bool);

  vtkGetMacro(RangeOctree, bool);
  vtkSetMacro(RangeOctree, bool);

  vtkGetMacro(PointMergeDistanceThreshold, double);
  vtkSetMacro(PointMergeDistanceThreshold, double);

  int FillInputPortInformation(int port, vtkInformation *info) override {

    switch(port) {
      case 0:
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
        return 1;
      case 1:
        info->Set(
          vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
        return 1;
    }

    return 0;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }

protected:
  ttkFiberSurface();

  ~ttkFiberSurface();

  TTK_SETUP();

private:
  bool RangeCoordinates, EdgeParameterization, EdgeIds, TetIds, CaseIds,
    RangeOctree, PointMerge;

  double PointMergeDistanceThreshold;

  std::string DataUcomponent, DataVcomponent, PolygonUcomponent,
    PolygonVcomponent;

  // NOTE: we assume here that this guy is small and that making a copy from
  // VTK is not an issue.
  std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>
    inputPolygon_;

  std::vector<ttk::FiberSurface::Vertex> outputVertexList_;
  std::vector<std::vector<ttk::FiberSurface::Vertex>> threadedVertexList_;
  std::vector<std::vector<ttk::FiberSurface::Triangle>> threadedTriangleList_;

  ttk::FiberSurface fiberSurface_;
};

#endif // _TTK_RANGEPOLYGON_H
