/// \ingroup vtk
/// \class ttkStableManifoldPersistence
/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
/// \date June 2021
///
/// \brief TTK VTK-filter for attaching to an input stable manifold (given by
/// the Morse-Smale complex module) its persistence (given by the persistence
/// diagram).
///
/// To change
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
/// \sa ttkMorseSmaleComplex
/// \sa ttkPersistenceDiagram

#pragma once

// VTK Module
#include <ttkStableManifoldPersistenceModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class vtkDataArray;
class vtkPolyData;
class vtkUnstructuredGrid;

class TTKSTABLEMANIFOLDPERSISTENCE_EXPORT ttkStableManifoldPersistence
  : public ttkAlgorithm {

public:
  static ttkStableManifoldPersistence *New();
  vtkTypeMacro(ttkStableManifoldPersistence, ttkAlgorithm);

protected:
  ttkStableManifoldPersistence();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int AttachPersistence(const std::vector<double> &vertex2persistence,
                        vtkDataSet *output) const;

  int BuildSimplex2PersistenceMap(
    vtkPolyData *criticalPoints,
    vtkUnstructuredGrid *persistenceDiagram,
    std::vector<double> &vertex2persistence) const;
};
