/// \ingroup vtk
/// \class ttkMergeTreeTemporalReductionEncoding
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// \brief TTK VTK-filter that wraps the ttk::MergeTreeTemporalReductionEncoding
/// module.
///
/// This VTK filter uses the ttk::MergeTreeTemporalReductionEncoding module to
/// compute a temporal reduction of a sequence of merge trees.
///
/// \param Input vtkMultiBlockDataSet Input trees
/// \param Output vtkMultiBlockDataSet Key frames.
/// \param Output vtkTable Interpolation coefficients.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MergeTreeTemporalReductionEncoding
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkMergeTreeTemporalReductionEncodingModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <MergeTreeTemporalReductionEncoding.h>

class TTKMERGETREETEMPORALREDUCTIONENCODING_EXPORT
  ttkMergeTreeTemporalReductionEncoding
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::MergeTreeTemporalReductionEncoding // and we inherit from the
                                                      // base class
{
private:
  // Execution options
  std::string TimeVariableName{""};

  // Input Options
  bool DoResampleToImage = false;

  // Output options
  double DistanceAxisStretch = 1.;
  bool OutputTrees = true;
  bool OutputSegmentation = false;
  bool PlanarLayout = false;
  bool BranchDecompositionPlanarLayout = false;
  double BranchSpacing = 1.;
  bool RescaleTreesIndividually = false;
  double DimensionSpacing = 1.;
  int DimensionToShift = 0;
  double ImportantPairs = 50.;
  double ImportantPairsSpacing = 1.;
  double NonImportantPairsSpacing = 1.;
  double NonImportantPairsProximity = 0.05;

  bool TemporalSubsamplingMDS = false;

  // ----------------------
  // Data for visualization
  // ----------------------
  // Trees
  std::vector<vtkUnstructuredGrid *> treesNodes;
  std::vector<vtkUnstructuredGrid *> treesArcs;
  std::vector<vtkDataSet *> treesSegmentation;
  // Output
  std::vector<std::vector<int>> treesNodeCorrMesh;
  std::vector<double> emptyTreeDistances;
  std::vector<MergeTree<double>> keyFrames;
  std::vector<int> removed;

  void setDataVisualization(int numInputs) {
    // Trees
    treesNodes = std::vector<vtkUnstructuredGrid *>(numInputs);
    treesArcs = std::vector<vtkUnstructuredGrid *>(numInputs);
    treesSegmentation = std::vector<vtkDataSet *>(numInputs);
  }

  void resetDataVisualization() {
    setDataVisualization(0);
    treesNodeCorrMesh = std::vector<std::vector<int>>();
    emptyTreeDistances = std::vector<double>();
    keyFrames = std::vector<MergeTree<double>>();
    removed = std::vector<int>();
  }

  bool isDataVisualizationFilled() {
    return treesNodeCorrMesh.size() != 0 and keyFrames.size() != 0
           and emptyTreeDistances.size() != 0 and removed.size() != 0;
  }

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */
  // Input Options
  void SetEpsilon1UseFarthestSaddle(bool epsilon1UseFarthestSaddle) {
    epsilon1UseFarthestSaddle_ = epsilon1UseFarthestSaddle;
    Modified();
    resetDataVisualization();
  }
  bool GetEpsilon1UseFarthestSaddle() {
    return epsilon1UseFarthestSaddle_;
  }

  void SetEpsilonTree1(double epsilonTree1) {
    epsilonTree1_ = epsilonTree1;
    Modified();
    resetDataVisualization();
  }
  double SetEpsilonTree1() {
    return epsilonTree1_;
  }

  void SetEpsilon2Tree1(double epsilon2Tree1) {
    epsilon2Tree1_ = epsilon2Tree1;
    Modified();
    resetDataVisualization();
  }
  double SetEpsilon2Tree1() {
    return epsilon2Tree1_;
  }

  void SetEpsilon3Tree1(double epsilon3Tree1) {
    epsilon3Tree1_ = epsilon3Tree1;
    Modified();
    resetDataVisualization();
  }
  double SetEpsilon3Tree1() {
    return epsilon3Tree1_;
  }

  void SetPersistenceThreshold(double persistenceThreshold) {
    persistenceThreshold_ = persistenceThreshold;
    Modified();
    resetDataVisualization();
  }
  double SetPersistenceThreshold() {
    return persistenceThreshold_;
  }

  void SetUseMinMaxPair(bool useMinMaxPair) {
    useMinMaxPair_ = useMinMaxPair;
    Modified();
    resetDataVisualization();
  }
  bool SetUseMinMaxPair() {
    return useMinMaxPair_;
  }

  void SetDeleteMultiPersPairs(bool deleteMultiPersPairs) {
    deleteMultiPersPairs_ = deleteMultiPersPairs;
    Modified();
    resetDataVisualization();
  }
  bool SetDeleteMultiPersPairs() {
    return deleteMultiPersPairs_;
  }

  // Execution Options
  void SetAssignmentSolver(int assignmentSolver) {
    assignmentSolverID_ = assignmentSolver;
    Modified();
    resetDataVisualization();
  }
  int GetAssignmentSolver() {
    return assignmentSolverID_;
  }

  void SetRemovalPercentage(double removePerc) {
    removalPercentage_ = removePerc;
    Modified();
    resetDataVisualization();
  }
  double GetRemovalPercentage() {
    return removalPercentage_;
  }

  void SetUseL2Distance(double useL2) {
    useL2Distance_ = useL2;
    Modified();
    resetDataVisualization();
  }
  double GetUseL2Distance() {
    return useL2Distance_;
  }

  void SetUseCustomTimeVariable(bool useCustomTime) {
    useCustomTimeVariable_ = useCustomTime;
    Modified();
    resetDataVisualization();
  }
  bool GetUseCustomTimeVariable() {
    return useCustomTimeVariable_;
  }

  void SetTimeVariableName(
    int idx, int port, int connection, int fieldAssociation, const char *name) {
    SetInputArrayToProcess(idx, port, connection, fieldAssociation, name);
    TimeVariableName = std::string(name);
    Modified();
    resetDataVisualization();
  }
  vtkGetMacro(TimeVariableName, std::string);

  // Output Options
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

  vtkSetMacro(ImportantPairsSpacing, double);
  vtkGetMacro(ImportantPairsSpacing, double);

  vtkSetMacro(NonImportantPairsSpacing, double);
  vtkGetMacro(NonImportantPairsSpacing, double);

  vtkSetMacro(NonImportantPairsProximity, double);
  vtkGetMacro(NonImportantPairsProximity, double);

  vtkSetMacro(DistanceAxisStretch, double);
  vtkGetMacro(DistanceAxisStretch, double);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreeTemporalReductionEncoding *New();
  vtkTypeMacro(ttkMergeTreeTemporalReductionEncoding, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkMergeTreeTemporalReductionEncoding();
  ~ttkMergeTreeTemporalReductionEncoding() override;

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
          std::vector<vtkMultiBlockDataSet *> &inputTrees);

  template <class dataType>
  int runCompute(std::vector<vtkMultiBlockDataSet *> &inputTrees);

  template <class dataType>
  int runOutput(vtkInformationVector *outputVector,
                std::vector<vtkMultiBlockDataSet *> &inputTrees);
};
