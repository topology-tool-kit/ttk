/// \ingroup vtk
/// \class ttkTopologicalSkeleton
/// \author Tanner Finken <finkent@arizona.edu>
/// \author Joshua A. Levine <josh@cs.arizona.edu>
/// \date Summer 2024.
///
/// \brief TTK VTK-filter that wraps the topologicalSkeleton processing package.
///
/// TTK module for the computation of Topological Skeleton.
/// The Topological Skeleton is a useful topological abstractions of vector
/// fields for data segmentation, feature extraction, etc.
///
/// \param Input Input vector field, defined as a point data vector field
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
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \b Related \b publication \n
/// "Localized Evaluation for Constructing Discrete Vector Fields" \n
/// Tanner Finken, Julien Tierny, Joshua A. Levine \n
/// IEEE VIS 2024.
///
/// \sa ttk::TopologicalSkeleton
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/discreteVectorFieldTopology/">Discrete
///   Vector Field Topology example</a> \n
///

#pragma once

// VTK Module
#include <ttkTopologicalSkeletonModule.h>

// ttk code includes
#include <TopologicalSkeleton.h>
#include <ttkAlgorithm.h>

class vtkPolyData;

class TTKTOPOLOGICALSKELETON_EXPORT ttkTopologicalSkeleton
  : public ttkAlgorithm,
    protected ttk::TopologicalSkeleton {

public:
  static ttkTopologicalSkeleton *New();

  vtkTypeMacro(ttkTopologicalSkeleton, ttkAlgorithm);

  vtkSetMacro(ComputeCriticalPoints, bool);
  vtkGetMacro(ComputeCriticalPoints, bool);

  vtkSetMacro(ComputeAscendingSeparatrices1, bool);
  vtkGetMacro(ComputeAscendingSeparatrices1, bool);

  vtkSetMacro(ComputeDescendingSeparatrices1, bool);
  vtkGetMacro(ComputeDescendingSeparatrices1, bool);

  vtkSetMacro(ComputeSaddleConnectors, bool);
  vtkGetMacro(ComputeSaddleConnectors, bool);

  vtkSetMacro(ComputeAttractingCycles1, bool);
  vtkGetMacro(ComputeAttractingCycles1, bool);

  vtkSetMacro(ComputeRepellingCycles1, bool);
  vtkGetMacro(ComputeRepellingCycles1, bool);

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

  vtkSetMacro(IterationThreshold, int);
  vtkGetMacro(IterationThreshold, int);

  vtkSetMacro(RunSimplification, int);
  vtkGetMacro(RunSimplification, int);

  vtkSetMacro(ReverseFullOrbit, bool);
  vtkGetMacro(ReverseFullOrbit, bool);

  vtkSetMacro(SimplificationThreshold, double);
  vtkGetMacro(SimplificationThreshold, double);

protected:
  template <typename scalarType, typename triangulationType>
  int dispatch(vtkDataArray *const inputVectors,
               vtkPolyData *const outputCriticalPoints,
               vtkPolyData *const outputSeparatrices1,
               vtkPolyData *const outputSeparatrices2,
               const triangulationType &triangulation);

  ttkTopologicalSkeleton();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int IterationThreshold{-1};
  OutputCriticalPoints criticalPoints_{};
  Output1Separatrices separatrices1_{};
  Output2Separatrices separatrices2_{};
  OutputManifold segmentations_{};
};
