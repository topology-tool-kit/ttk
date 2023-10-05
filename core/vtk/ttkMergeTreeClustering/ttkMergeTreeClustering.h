/// \ingroup vtk
/// \class ttkMergeTreeClustering
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// \brief TTK VTK-filter that wraps the ttk::MergeTreeClustering module.
///
/// This VTK filter uses the ttk::MergeTreeClustering module to compute the edit
/// distance between two merge trees.
///
/// \param Input vtkMultiBlockDataSet Input trees
/// \param Input (optional) vtkMultiBlockDataSet Input trees
/// \param Output vtkMultiBlockDataSet Input trees (processed)
/// \param Output vtkMultiBlockDataSet Centroids trees
/// \param Output vtkMultiBlockDataSet Matchings
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Wasserstein Distances, Geodesics and Barycenters of Merge Trees" \n
/// Mathieu Pont, Jules Vidal, Julie Delon, Julien Tierny.\n
/// Proc. of IEEE VIS 2021.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// \sa ttk::MergeTreeClustering
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeClustering/">Merge
///   Tree Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeFeatureTracking/">Merge
///   Tree Feature Tracking example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n

#pragma once

// VTK Module
#include <ttkMergeTreeClusteringModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <MergeTreeBarycenter.h>
#include <MergeTreeClustering.h>
#include <MergeTreeDistance.h>

class TTKMERGETREECLUSTERING_EXPORT ttkMergeTreeClustering
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
//, protected ttk::MergeTreeDistance // and we inherit from the base class
{
private:
  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */
  // Input Options
  bool Epsilon1UseFarthestSaddle = false;
  double EpsilonTree1 = 5.;
  double EpsilonTree2 = 5.;
  double Epsilon2Tree1 = 95.;
  double Epsilon2Tree2 = 95.;
  double Epsilon3Tree1 = 90.;
  double Epsilon3Tree2 = 90.;
  double PersistenceThreshold = 0.;
  bool DeleteMultiPersPairs = false;
  bool UseMinMaxPair = true;
  bool IsPersistenceDiagram = false;

  // Execution Options
  int Backend = 0;
  double Alpha = 0.5;
  int AssignmentSolver = 0;
  bool BranchDecomposition = true;
  bool NormalizedWasserstein = true;
  bool KeepSubtree = false;
  bool oldBD = BranchDecomposition;
  bool oldNW = NormalizedWasserstein;
  bool oldKS = KeepSubtree;
  double JoinSplitMixtureCoefficient = 0.5;
  bool ComputeBarycenter = false;
  unsigned int NumberOfBarycenters = 1;
  double BarycenterSizeLimitPercent = 0.0;
  bool Deterministic = false;
  int pathMetric = 0;
  int branchMetric = 0;
  int baseModule = 0;
  double NonMatchingWeight = 1.0;

  // Output Options
  bool OutputTrees = true;
  bool OutputSegmentation = false;
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
  bool BarycenterPositionAlpha = false;
  std::string ExcludeImportantPairsLower = "";
  std::string ExcludeImportantPairsHigher = "";

  // ----------------------
  // Data for visualization
  // ----------------------
  // Trees
  std::vector<ttk::ftm::MergeTree<double>> intermediateSTrees,
    intermediateSTrees2;
  std::vector<vtkUnstructuredGrid *> treesNodes, treesNodes2;
  std::vector<vtkUnstructuredGrid *> treesArcs, treesArcs2;
  std::vector<vtkDataSet *> treesSegmentation, treesSegmentation2;

  // Matching
  std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>
    outputMatching;
  std::vector<std::vector<
    std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>>>
    outputMatchingBarycenter, outputMatchingBarycenter2;

  // Barycenter
  std::vector<ttk::ftm::MergeTree<double>> barycentersS, barycentersS2;
  std::vector<int> clusteringAssignment;

  // Node correspondence
  std::vector<std::vector<int>> trees1NodeCorrMesh, trees2NodeCorrMesh;

  // Distances
  std::vector<double> finalDistances;

  void setDataVisualization(int numInputs, int numInputs2) {
    // Trees
    intermediateSTrees = std::vector<ttk::ftm::MergeTree<double>>(numInputs);
    intermediateSTrees2 = std::vector<ttk::ftm::MergeTree<double>>(numInputs2);
    treesNodes = std::vector<vtkUnstructuredGrid *>(numInputs);
    treesNodes2 = std::vector<vtkUnstructuredGrid *>(numInputs2);
    treesArcs = std::vector<vtkUnstructuredGrid *>(numInputs);
    treesArcs2 = std::vector<vtkUnstructuredGrid *>(numInputs2);
    treesSegmentation = std::vector<vtkDataSet *>(numInputs);
    treesSegmentation2 = std::vector<vtkDataSet *>(numInputs2);

    // Matching
    outputMatchingBarycenter = std::vector<std::vector<
      std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>>>(
      NumberOfBarycenters,
      std::vector<
        std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>>(
        numInputs));

    outputMatchingBarycenter2 = std::vector<std::vector<
      std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>>>(
      NumberOfBarycenters,
      std::vector<
        std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>>(
        numInputs2));

    // Barycenter
    barycentersS
      = std::vector<ttk::ftm::MergeTree<double>>(NumberOfBarycenters);
    clusteringAssignment = std::vector<int>(numInputs, 0);
  }

  void resetDataVisualization() {
    setDataVisualization(0, 0);
    trees1NodeCorrMesh = std::vector<std::vector<int>>();
    trees2NodeCorrMesh = std::vector<std::vector<int>>();
    finalDistances = std::vector<double>();
  }

  bool isDataVisualizationFilled() {
    return trees1NodeCorrMesh.size() != 0 and finalDistances.size() != 0;
  }

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */
  // Input Options
  void SetEpsilon1UseFarthestSaddle(bool epsilon1UseFarthestSaddle) {
    Epsilon1UseFarthestSaddle = epsilon1UseFarthestSaddle;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(Epsilon1UseFarthestSaddle, bool);

  void SetEpsilonTree1(double epsilonTree1) {
    EpsilonTree1 = epsilonTree1;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(EpsilonTree1, double);

  void SetEpsilon2Tree1(double epsilon2Tree1) {
    Epsilon2Tree1 = epsilon2Tree1;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(Epsilon2Tree1, double);

  void SetEpsilon3Tree1(double epsilon3Tree1) {
    Epsilon3Tree1 = epsilon3Tree1;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(Epsilon3Tree1, double);

  void SetPersistenceThreshold(double persistenceThreshold) {
    PersistenceThreshold = persistenceThreshold;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(PersistenceThreshold, double);

  void SetUseMinMaxPair(bool useMinMaxPair) {
    UseMinMaxPair = useMinMaxPair;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(UseMinMaxPair, bool);

  void SetDeleteMultiPersPairs(bool deleteMultiPersPairs) {
    DeleteMultiPersPairs = deleteMultiPersPairs;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(DeleteMultiPersPairs, bool);

  // Execution Options
  void SetBackend(int newBackend) {
    if(Backend == 2) { // Custom
      oldBD = BranchDecomposition;
      oldNW = NormalizedWasserstein;
      oldKS = KeepSubtree;
    }
    if(newBackend == 2) { // Custom
      BranchDecomposition = oldBD;
      NormalizedWasserstein = oldNW;
      KeepSubtree = oldKS;
    }
    Backend = newBackend;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(Backend, int);

  void SetAlpha(double alpha) {
    Alpha = 1 - alpha;
    Alpha = std::min(1 - 1e-6, Alpha);
    Alpha = std::max(1e-6, Alpha);
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(Alpha, double);

  void SetAssignmentSolver(int assignmentSolver) {
    AssignmentSolver = assignmentSolver;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(AssignmentSolver, int);

  void SetBranchDecomposition(bool branchDecomposition) {
    BranchDecomposition = branchDecomposition;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(BranchDecomposition, bool);

  void SetDeterministic(bool deterministic) {
    Deterministic = deterministic;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(Deterministic, bool);

  void SetNormalizedWasserstein(bool normalizedWasserstein) {
    NormalizedWasserstein = normalizedWasserstein;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(NormalizedWasserstein, bool);

  void SetKeepSubtree(bool keepSubtree) {
    KeepSubtree = keepSubtree;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(KeepSubtree, bool);

  void SetJoinSplitMixtureCoefficient(double joinSplitMixtureCoefficient) {
    JoinSplitMixtureCoefficient = joinSplitMixtureCoefficient;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(JoinSplitMixtureCoefficient, double);

  void SetComputeBarycenter(bool computeBarycenter) {
    ComputeBarycenter = computeBarycenter;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(ComputeBarycenter, bool);

  void SetNumberOfBarycenters(unsigned int numberOfBarycenters) {
    NumberOfBarycenters = numberOfBarycenters;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(NumberOfBarycenters, unsigned int);

  void SetBarycenterSizeLimitPercent(double percent) {
    BarycenterSizeLimitPercent = percent;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(BarycenterSizeLimitPercent, double);

  void SetBranchMetric(int m) {
    branchMetric = m;
    Modified();
  }

  void SetPathMetric(int m) {
    pathMetric = m;
    Modified();
  }

  void SetNonMatchingWeight(double weight) {
    NonMatchingWeight = weight;
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(NonMatchingWeight, double);

  // Output Options
  vtkSetMacro(BarycenterPositionAlpha, bool);
  vtkGetMacro(BarycenterPositionAlpha, bool);

  vtkSetMacro(OutputTrees, bool);
  vtkGetMacro(OutputTrees, bool);

  vtkSetMacro(OutputSegmentation, bool);
  vtkGetMacro(OutputSegmentation, bool);

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

  vtkSetMacro(ExcludeImportantPairsLower, const std::string &);
  vtkGetMacro(ExcludeImportantPairsLower, std::string);

  vtkSetMacro(ExcludeImportantPairsHigher, const std::string &);
  vtkGetMacro(ExcludeImportantPairsHigher, std::string);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreeClustering *New();
  vtkTypeMacro(ttkMergeTreeClustering, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkMergeTreeClustering();
  ~ttkMergeTreeClustering() override;

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
          std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
          std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2);

  template <class dataType>
  int runCompute(
    vtkInformationVector *outputVector,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2);

  template <class dataType>
  int runOutput(
    vtkInformationVector *outputVector,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2);
};
