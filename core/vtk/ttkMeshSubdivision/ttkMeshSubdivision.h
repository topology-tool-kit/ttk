/// \ingroup vtk
/// \class ttkMeshSubdivision
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date December 2015.
///
/// \brief TTK VTK-filter for the computation of mesh subdivisions (without
/// vertex displacement).
///
/// \param Input Input mesh (vtkUnstructuredGrid)
/// \param Output Output mesh (vtkUnstructuredGrid)
///
/// This filter subdivides an input mesh with a strategy inspired by Discrete
/// Morse theory. It does not modify the position of the original vertices.
///
/// The filter inserts a new vertex at the barycenter of each d-cell C and
/// connects it with an edge to the new vertices inserted in the (d-1) and
/// (d+1)-cells that are faces or co-faces of C.
///
/// In practice, for surface meshes, this subdivision scheme corresponds to a
/// Catmull-Clark subdivision (without vertex displacement). It will turn a
/// triangle-mesh into a quad mesh, a quad mesh into a (finer) quad-mesh, a
/// tetrahedral mesh into a hexahedral mesh, a hexahedral mesh into a (finer)
/// hexadrehal mesh, etc. Generally, it will turn any 2D mesh into a quad mesh
/// and any 3D mesh into a hexadrehal mesh.
///
/// This filter assumes that all the cells of the input mesh are of the same
/// type. Also, the filter creates duplicate points, to be merged after the
/// fact with "Clean to Grid" under ParaView or vtkMergePoints for instance.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
#pragma once

#include <ttkAlgorithm.h>
#include <ttkMeshSubdivisionModule.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKMESHSUBDIVISION_EXPORT ttkMeshSubdivision : public ttkAlgorithm {
private:
  int IterationNumber;

public:
  vtkSetMacro(IterationNumber, int);
  vtkGetMacro(IterationNumber, int);

  vtkTypeMacro(ttkMeshSubdivision, ttkAlgorithm);
  static ttkMeshSubdivision *New();

protected:
  ttkMeshSubdivision();
  ~ttkMeshSubdivision();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};