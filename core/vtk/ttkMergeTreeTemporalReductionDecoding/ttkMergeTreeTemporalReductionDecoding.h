/// \ingroup vtk
/// \class ttkMergeTreeTemporalReductionDecoding
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// \brief TTK VTK-filter that wraps the ttk::MergeTreeTemporalReductionDecoding
/// module.
///
/// This VTK filter uses the ttk::MergeTreeTemporalReductionDecoding module to
/// compute the reconstruction of a reduced sequence of merge trees.
///
/// \param Input vtkMultiBlockDataSet Key frames.
/// \param Input vtkTable Interpolation coefficients.
/// \param Output vtkMultiBlockDataSet Input trees
///
/// See the corresponding standalone program for a usage example:
///   - standalone/MergeTreeTemporalReductionDecoding/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MergeTreeTemporalReductionDecoding
/// \sa ttkAlgorithm
///
/// \b Related \b publication \n
/// "Wasserstein Distances, Geodesics and Barycenters of Merge Trees" \n
/// Mathieu Pont, Jules Vidal, Julie Delon, Julien Tierny.\n
/// Proc. of IEEE VIS 2021.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeTemporalReduction/">Merge
///   Tree Temporal Reduction</a> \n

#pragma once

// VTK Module
#include <ttkMergeTreeTemporalReductionDecodingModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

/* Note on including VTK modules
 *
 * Each VTK module that you include a header from needs to be specified in this
 * module's vtk.module file, either in the DEPENDS or PRIVATE_DEPENDS (if the
 * header is included in the cpp file only) sections.
 *
 * In order to find the corresponding module, check its location within the VTK
 * source code. The VTK module name is composed of the path to the header. You
 * can also find the module name within the vtk.module file located in the same
 * directory as the header file.
 *
 * For example, vtkSphereSource.h is located in directory VTK/Filters/Sources/,
 * so its corresponding VTK module is called VTK::FiltersSources. In this case,
 * the vtk.module file would need to be extended to
 *
 * NAME
 *   ttkMergeTreeTemporalReductionDecoding
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <MergeTreeTemporalReductionDecoding.h>

// VTK Includes
#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>

class TTKMERGETREETEMPORALREDUCTIONDECODING_EXPORT
  ttkMergeTreeTemporalReductionDecoding
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::MergeTreeTemporalReductionDecoding // and we inherit from the
                                                      // base class
{
private:
  using idNode = ttk::ftm::idNode;

  // Output options
  bool OutputTrees = true;
  bool PlanarLayout = false;
  bool BranchDecompositionPlanarLayout = false;
  double BranchSpacing = 1.;
  bool RescaleTreesIndividually = false;
  double DimensionSpacing = 1.;
  int DimensionToShift = 0;
  double ImportantPairs = 50.;
  int MaximumImportantPairs = 0;
  int MinimumImportantPairs = 0;
  double ImportantPairsSpacing = 1.;
  double NonImportantPairsSpacing = 1.;
  double NonImportantPairsProximity = 0.05;

  // ----------------------
  // Data for visualization
  // ----------------------
  // Trees
  std::vector<vtkUnstructuredGrid *> treesNodes;
  std::vector<vtkUnstructuredGrid *> treesArcs;
  std::vector<vtkDataSet *> treesSegmentation;
  // Output
  std::vector<std::vector<int>> treesNodeCorrMesh;
  std::vector<ttk::ftm::MergeTree<double>> intermediateSTrees;
  std::vector<std::vector<std::tuple<idNode, idNode, double>>> allMatching;

  void setDataVisualization(int numInputs) {
    // Trees
    treesNodes = std::vector<vtkUnstructuredGrid *>(numInputs);
    treesArcs = std::vector<vtkUnstructuredGrid *>(numInputs);
    treesSegmentation = std::vector<vtkDataSet *>(numInputs);
  }

  void resetDataVisualization() {
    setDataVisualization(0);
    treesNodeCorrMesh = std::vector<std::vector<int>>();
    intermediateSTrees = std::vector<ttk::ftm::MergeTree<double>>();
    allMatching
      = std::vector<std::vector<std::tuple<idNode, idNode, double>>>();
  }

  bool isDataVisualizationFilled() {
    return treesNodeCorrMesh.size() != 0 and intermediateSTrees.size() != 0
           and allMatching.size() != 0;
  }

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */
  // Input Options
  void SetAssignmentSolver(int assignmentSolver) {
    assignmentSolverID_ = assignmentSolver;
    Modified();
    resetDataVisualization();
  }
  int GetAssignmentSolver() {
    return assignmentSolverID_;
  }

  // Output Options
  vtkSetMacro(OutputTrees, bool);
  vtkGetMacro(OutputTrees, bool);

  vtkSetMacro(PlanarLayout, bool);
  vtkGetMacro(PlanarLayout, bool);

  vtkSetMacro(BranchDecompositionPlanarLayout, bool);
  vtkGetMacro(BranchDecompositionPlanarLayout, bool);

  vtkSetMacro(BranchSpacing, double);
  vtkGetMacro(BranchSpacing, double);

  vtkSetMacro(RescaleTreesIndividually, bool);
  vtkGetMacro(RescaleTreesIndividually, bool);

  vtkSetMacro(DimensionSpacing, double);
  vtkGetMacro(DimensionSpacing, double);

  vtkSetMacro(DimensionToShift, int);
  vtkGetMacro(DimensionToShift, int);

  vtkSetMacro(ImportantPairs, double);
  vtkGetMacro(ImportantPairs, double);

  vtkSetMacro(MaximumImportantPairs, int);
  vtkGetMacro(MaximumImportantPairs, int);

  vtkSetMacro(MinimumImportantPairs, int);
  vtkGetMacro(MinimumImportantPairs, int);

  vtkSetMacro(ImportantPairsSpacing, double);
  vtkGetMacro(ImportantPairsSpacing, double);

  vtkSetMacro(NonImportantPairsSpacing, double);
  vtkGetMacro(NonImportantPairsSpacing, double);

  vtkSetMacro(NonImportantPairsProximity, double);
  vtkGetMacro(NonImportantPairsProximity, double);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreeTemporalReductionDecoding *New();
  vtkTypeMacro(ttkMergeTreeTemporalReductionDecoding, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkMergeTreeTemporalReductionDecoding();
  ~ttkMergeTreeTemporalReductionDecoding() override;

  /**
   * Specify the input data type of each input port
   * (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * Specify the data object type of each output port
   * (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * Pass VTK data to the base code and convert base code output to VTK
   * (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <class dataType>
  int run(vtkInformationVector *outputVector,
          std::vector<vtkMultiBlockDataSet *> &inputTrees,
          std::vector<std::tuple<double, int, int, int, int>> &coefs,
          std::vector<bool> &interpolatedTrees);

  template <class dataType>
  int runCompute(std::vector<vtkMultiBlockDataSet *> &inputTrees,
                 std::vector<std::tuple<double, int, int, int, int>> &coefs);

  template <class dataType>
  int runOutput(vtkInformationVector *outputVector,
                std::vector<vtkMultiBlockDataSet *> &inputTrees,
                std::vector<std::tuple<double, int, int, int, int>> &coefs,
                std::vector<bool> &interpolatedTrees);
};
