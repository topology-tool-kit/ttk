/// \ingroup vtk
/// \class ttkTopologicalSimplification
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date February 2016
///
/// \brief TTK VTK-filter for the topological simplification of scalar
/// data.
///
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
/// can specified to the ttkMergeTree, vtkIntegralLines, or
/// vtkScalarFieldCriticalPoints filters.
///
/// Also, this filter can be given a specific input vertex offset.
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
/// \b Related \b publications \n
/// "Generalized Topological Simplification of Scalar Fields on Surfaces" \n
/// Julien Tierny, Valerio Pascucci \n
/// IEEE Transactions on Visualization and Computer Graphics.\n
/// Proc. of IEEE VIS 2012.
///
/// "Localized Topological Simplification of Scalar Data" \n
/// Jonas Lukasczyk, Christoph Garth, Ross Maciejewski, Julien Tierny \n
/// IEEE Transactions on Visualization and Computer Graphics.\n
/// Proc. of IEEE VIS 2020.
///
/// "A Practical Solver for Scalar Data Topological Simplification"\n
/// Mohamed Kissi, Mathieu Pont, Joshua A. Levine, Julien Tierny\n
/// IEEE Transactions on Visualization and Computer Graphics.\n
/// Proc. of IEEE VIS 2024.
///
/// \sa ttkTopologicalSimplificationByPersistence
/// \sa ttkScalarFieldCriticalPoints
/// \sa ttkIntegralLines
/// \sa ttkMergeTree
/// \sa ttkMorseSmaleComplex
/// \sa ttkIdentifiers
/// \sa ttk::TopologicalSimplification
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearning/">1-Manifold
///   Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearningCircles/">1-Manifold
///   Learning Circles example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/2manifoldLearning/">
///   2-Manifold Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/BuiltInExample1/">BuiltInExample1
///   example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/contourTreeAlignment/">Contour
///   Tree Alignment example</a> \n
///   - <a href="https://topology-tool-kit.github.io/examples/ctBones/">CT Bones
///   example</a> \n
///   - <a href="https://topology-tool-kit.github.io/examples/dragon/">Dragon
///   example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/harmonicSkeleton/">
///   Harmonic Skeleton example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/imageProcessing/">Image
///   Processing example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/interactionSites/">
///   Interaction sites</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/karhunenLoveDigits64Dimensions/">Karhunen-Love
///   Digits 64-Dimensions example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morsePersistence/">Morse
///   Persistence example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleQuadrangulation/">Morse-Smale
///   Quadrangulation example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 0 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 1 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 2 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 3 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 4 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/tectonicPuzzle/">Tectonic
///   Puzzle example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/topologicalOptimization_darkSky/">Topological
///   Optimization DarkSky example</a>\n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/topologicalOptimization_pegasus/">Topological
///   Optimization for Pegasus Genus Repair example</a>\n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/topologicalOptimization_torus/">Topological
///   Optimization for Torus Repair example</a>\n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/tribute/">Tribute
///   example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/uncertainStartingVortex/">
///   Uncertain Starting Vortex example</a> \n
///

#pragma once

// VTK Module
#include <ttkTopologicalSimplificationModule.h>

// ttk code includes
#include <TopologicalSimplification.h>
#include <ttkAlgorithm.h>
#include <ttkPersistenceDiagramUtils.h>
#include <ttkUtils.h>

class vtkDataArray;

class TTKTOPOLOGICALSIMPLIFICATION_EXPORT ttkTopologicalSimplification
  : public ttkAlgorithm,
    protected ttk::TopologicalSimplification {

public:
  static ttkTopologicalSimplification *New();
  vtkTypeMacro(ttkTopologicalSimplification, ttkAlgorithm);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  vtkSetMacro(ConsiderIdentifierAsBlackList, bool);
  vtkGetMacro(ConsiderIdentifierAsBlackList, bool);

  vtkSetMacro(AddPerturbation, bool);
  vtkGetMacro(AddPerturbation, bool);

  vtkSetMacro(ForceInputVertexScalarField, bool);
  vtkGetMacro(ForceInputVertexScalarField, bool);

  vtkSetMacro(Method, int);
  vtkGetMacro(Method, int);

  vtkSetMacro(PersistenceThreshold, double);
  vtkGetMacro(PersistenceThreshold, double);

  vtkSetMacro(UseFastPersistenceUpdate, bool);
  vtkGetMacro(UseFastPersistenceUpdate, bool);

  vtkSetMacro(FastAssignmentUpdate, bool);
  vtkGetMacro(FastAssignmentUpdate, bool);

  vtkSetMacro(EpochNumber, int);
  vtkGetMacro(EpochNumber, int);

  vtkSetMacro(PDCMethod, int);
  vtkGetMacro(PDCMethod, int);

  vtkSetMacro(MethodOptimization, int);
  vtkGetMacro(MethodOptimization, int);

  vtkSetMacro(FinePairManagement, int);
  vtkGetMacro(FinePairManagement, int);

  vtkSetMacro(ChooseLearningRate, bool);
  vtkGetMacro(ChooseLearningRate, bool);

  vtkSetMacro(LearningRate, double);
  vtkGetMacro(LearningRate, double);

  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);

  vtkSetMacro(CoefStopCondition, double);
  vtkGetMacro(CoefStopCondition, double);

  vtkSetMacro(OptimizationWithoutMatching, bool);
  vtkGetMacro(OptimizationWithoutMatching, bool);

  vtkSetMacro(ThresholdMethod, int);
  vtkGetMacro(ThresholdMethod, int);

  vtkSetMacro(Threshold, double);
  vtkGetMacro(Threshold, double);

  vtkSetMacro(LowerThreshold, int);
  vtkGetMacro(LowerThreshold, int);

  vtkSetMacro(UpperThreshold, int);
  vtkGetMacro(UpperThreshold, int);

  vtkSetMacro(PairTypeToDelete, int);
  vtkGetMacro(PairTypeToDelete, int);

  vtkSetMacro(ConstraintAveraging, bool);
  vtkGetMacro(ConstraintAveraging, bool);

  vtkSetMacro(PrintFrequency, int);
  vtkGetMacro(PrintFrequency, int);

protected:
  ttkTopologicalSimplification();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool ForceInputVertexScalarField{false};
  bool ForceInputOffsetScalarField{false};
  bool ConsiderIdentifierAsBlackList{false};
  bool AddPerturbation{false};
  int Method{0};
  double PersistenceThreshold{0};
};
