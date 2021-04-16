/// \ingroup vtk
/// \class ttkMorseSmaleComplex
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK VTK-filter that wraps the morseSmaleComplex processing package.
///
/// TTK module for the computation of Morse-Smale complexes.
/// Morse-Smale complexes are useful topological abstractions of scalar
/// fields for data segmentation, feature extraction, etc.
///
/// \b Related \b publication \n
/// "Parallel Computation of 3D Morse-Smale Complexes" \n
/// Nithin Shivashankar, Vijay Natarajan \n
/// Proc. of EuroVis 2012. \n
/// Computer Graphics Forum, 2012.
///
/// \param Input Input scalar field, defined as a point data scalar field
/// attached to a geometry, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output0 Output critical points (vtkPolyData)
/// \param Output1 Output 1-separatrices (vtkPolyData)
/// \param Output2 Output 2-separatrices (vtkPolyData)
/// \param Output3 Output data segmentation (vtkDataSet)
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The optional offset array can be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 1 (FIXED: the second array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the offset array)
/// \note: To use this optional array, `ForceInputOffsetScalarField` needs to be
/// enabled with the setter `setForceInputOffsetScalarField()'.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MorseSmaleComplex
///

#pragma once

// VTK Module
#include <ttkMorseSmaleComplexModule.h>

// ttk code includes
#include <MorseSmaleComplex.h>
#include <ttkAlgorithm.h>

class vtkPolyData;

class TTKMORSESMALECOMPLEX_EXPORT ttkMorseSmaleComplex
  : public ttkAlgorithm,
    protected ttk::MorseSmaleComplex {

public:
  static ttkMorseSmaleComplex *New();

  vtkTypeMacro(ttkMorseSmaleComplex, ttkAlgorithm);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  vtkSetMacro(IterationThreshold, int);
  vtkGetMacro(IterationThreshold, int);

  vtkSetMacro(ComputeCriticalPoints, bool);
  vtkGetMacro(ComputeCriticalPoints, bool);

  vtkSetMacro(ComputeAscendingSeparatrices1, bool);
  vtkGetMacro(ComputeAscendingSeparatrices1, bool);

  vtkSetMacro(ComputeDescendingSeparatrices1, bool);
  vtkGetMacro(ComputeDescendingSeparatrices1, bool);

  vtkSetMacro(ComputeSaddleConnectors, bool);
  vtkGetMacro(ComputeSaddleConnectors, bool);

  vtkSetMacro(ComputeAscendingSeparatrices2, bool);
  vtkGetMacro(ComputeAscendingSeparatrices2, bool);

  vtkSetMacro(ComputeDescendingSeparatrices2, bool);
  vtkGetMacro(ComputeDescendingSeparatrices2, bool);

  vtkSetMacro(ComputeAscendingSegmentation, bool);
  vtkGetMacro(ComputeAscendingSegmentation, bool);

  vtkSetMacro(ComputeDescendingSegmentation, bool);
  vtkGetMacro(ComputeDescendingSegmentation, bool);

  vtkSetMacro(ComputeFinalSegmentation, bool);
  vtkGetMacro(ComputeFinalSegmentation, bool);

  vtkSetMacro(ReturnSaddleConnectors, int);
  vtkGetMacro(ReturnSaddleConnectors, int);

  vtkSetMacro(SaddleConnectorsPersistenceThreshold, double);
  vtkGetMacro(SaddleConnectorsPersistenceThreshold, double);

protected:
  template <typename scalarType, typename triangulationType>
  int dispatch(vtkDataArray *const inputScalars,
               vtkDataArray *const inputOffsets,
               vtkPolyData *const outputCriticalPoints,
               vtkPolyData *const outputSeparatrices1,
               vtkPolyData *const outputSeparatrices2,
               const triangulationType &triangulation);

  ttkMorseSmaleComplex();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool ForceInputOffsetScalarField{};
  int IterationThreshold{-1};
  bool ComputeCriticalPoints{true};
  bool ComputeAscendingSeparatrices1{true};
  bool ComputeDescendingSeparatrices1{true};
  bool ComputeSaddleConnectors{true};
  bool ComputeAscendingSeparatrices2{false};
  bool ComputeDescendingSeparatrices2{false};
  bool ComputeAscendingSegmentation{true};
  bool ComputeDescendingSegmentation{true};
  bool ComputeFinalSegmentation{true};
  int ReturnSaddleConnectors{false};
  double SaddleConnectorsPersistenceThreshold{0.0};

  // critical points
  std::vector<std::array<float, 3>> criticalPoints_points{};
  std::vector<char> criticalPoints_points_cellDimensions{};
  std::vector<ttk::SimplexId> criticalPoints_points_cellIds{};
  std::vector<char> criticalPoints_points_isOnBoundary{};
  std::vector<ttk::SimplexId> criticalPoints_points_PLVertexIdentifiers{};
  std::vector<ttk::SimplexId> criticalPoints_points_manifoldSize{};

  // 1-separatrices data
  std::vector<float> separatrices1_points{};
  std::vector<char> separatrices1_points_smoothingMask{};
  std::vector<char> separatrices1_points_cellDimensions{};
  std::vector<ttk::SimplexId> separatrices1_points_cellIds{};
  std::vector<ttk::SimplexId> separatrices1_cells_connectivity{};
  std::vector<ttk::SimplexId> separatrices1_cells_sourceIds{};
  std::vector<ttk::SimplexId> separatrices1_cells_destinationIds{};
  std::vector<ttk::SimplexId> separatrices1_cells_separatrixIds{};
  std::vector<char> separatrices1_cells_separatrixTypes{};
  std::vector<char> separatrices1_cells_isOnBoundary{};

  // 2-separatrices data
  std::vector<float> separatrices2_points{};
  std::vector<ttk::SimplexId> separatrices2_cells_offsets{};
  std::vector<ttk::SimplexId> separatrices2_cells_connectivity{};
  std::vector<ttk::SimplexId> separatrices2_cells_sourceIds{};
  std::vector<ttk::SimplexId> separatrices2_cells_separatrixIds{};
  std::vector<char> separatrices2_cells_separatrixTypes{};
  std::vector<char> separatrices2_cells_isOnBoundary{};
};
