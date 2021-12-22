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
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/builtInExample2/">
///   Builtin example 2</a> \n

#pragma once

// VTK Module
#include <ttkFiberSurfaceModule.h>

// ttk code includes
#include <FiberSurface.h>
#include <ttkAlgorithm.h>

class TTKFIBERSURFACE_EXPORT ttkFiberSurface : public ttkAlgorithm,
                                               protected ttk::FiberSurface {

public:
  static ttkFiberSurface *New();
  vtkTypeMacro(ttkFiberSurface, ttkAlgorithm);

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

protected:
  ttkFiberSurface();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <typename VTK_T1, typename VTK_T2>
  int dispatch(ttk::Triangulation *const triangulation);

private:
  bool RangeCoordinates{true}, EdgeParameterization{true}, EdgeIds{true},
    TetIds{true}, CaseIds{true}, RangeOctree{true}, PointMerge{false};

  double PointMergeDistanceThreshold{0.000001};

  // NOTE: we assume here that this guy is small and that making a copy from
  // VTK is not an issue.
  std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>
    inputPolygon_{};

  std::vector<ttk::FiberSurface::Vertex> outputVertexList_{};
  std::vector<std::vector<ttk::FiberSurface::Vertex>> threadedVertexList_{};
  std::vector<std::vector<ttk::FiberSurface::Triangle>> threadedTriangleList_{};
};
