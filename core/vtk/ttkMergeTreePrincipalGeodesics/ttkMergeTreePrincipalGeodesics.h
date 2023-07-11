/// \ingroup vtk
/// \class ttkMergeTreePrincipalGeodesics
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2022.
///
/// \brief TTK VTK-filter that wraps the ttk::MergeTreePrincipalGeodesics
/// module.
///
/// This VTK filter uses the ttk::MergeTreePrincipalGeodesics module to compute
/// Principal Geodesic Analysis on the space of merge trees or persistence
/// diagrams, that is, a set of orthogonal geodesic axes defining a basis with
/// the barycenter as origin.
///
/// \param Input vtkMultiBlockDataSet Input trees
/// \param Input (optional) vtkMultiBlockDataSet Input trees
/// \param Output vtkMultiBlockDataSet Barycenter
/// \param Output vtkTable Coefficients
/// \param Output vtkTable Geodesics Vectors
/// \param Output vtkTable Correlation Matrix
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the corresponding standalone program for a usage example:
///   - standalone/MergeTreePrincipalGeodesics/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication: \n
/// "Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)" \n
/// Mathieu Pont, Jules Vidal, Julien Tierny.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2022
///
/// \sa ttk::MergeTreePrincipalGeodesics
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n

#pragma once

// VTK Module
#include <ttkMergeTreePrincipalGeodesicsModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <FTMTree.h>
#include <MergeTreePrincipalGeodesics.h>
#include <ttkMergeTreeVisualization.h>

class TTKMERGETREEPRINCIPALGEODESICS_EXPORT ttkMergeTreePrincipalGeodesics
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::MergeTreePrincipalGeodesics // and we inherit from the base
                                               // class
{
private:
  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */
  // Input options
  double oldEpsilonTree1;
  // Output options

  // ----------------------
  // Data for visualization
  // ----------------------
  // Trees
  std::vector<ttk::ftm::MergeTree<double>> intermediateDTrees;
  std::vector<vtkUnstructuredGrid *> treesNodes, treesNodes2;
  std::vector<vtkUnstructuredGrid *> treesArcs, treesArcs2;
  std::vector<vtkDataSet *> treesSegmentation, treesSegmentation2;
  // Output
  // -> base class

  void setDataVisualization(int ttkNotUsed(numInputs),
                            int ttkNotUsed(numInputs2)) {
  }

  void resetDataVisualization() {
    setDataVisualization(0, 0);
    if(not keepState_) {
      barycenter_ = ttk::ftm::MergeTree<double>();
      allTs_.clear();
    }
  }

  bool isDataVisualizationFilled() {
    return allTs_.size() != 0 and not keepState_;
  }

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */
  // Input Options
  void SetNormalizedWasserstein(bool nW) {
    normalizedWasserstein_ = nW;
    Modified();
    resetDataVisualization();
  }
  bool GetNormalizedWasserstein() {
    return normalizedWasserstein_;
  }

  void SetNumberOfGeodesics(unsigned int numberOfGeodesics) {
    numberOfGeodesics_ = numberOfGeodesics;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetNumberOfGeodesics() {
    return numberOfGeodesics_;
  }

  void SetNumberOfProjectionIntervals(unsigned int intervals) {
    k_ = intervals;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetNumberOfProjectionIntervals() {
    return k_;
  }

  void SetNumberOfProjectionSteps(unsigned int steps) {
    noProjectionStep_ = steps;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetNumberOfProjectionSteps() {
    return noProjectionStep_;
  }

  void SetBarycenterSizeLimitPercent(double percent) {
    barycenterSizeLimitPercent_ = percent;
    Modified();
    resetDataVisualization();
  }
  double GetBarycenterSizeLimitPercent() {
    return barycenterSizeLimitPercent_;
  }

  void SetDeterministic(bool deterministic) {
    deterministic_ = deterministic;
    Modified();
    resetDataVisualization();
  }
  bool GetDeterministic() {
    return deterministic_;
  }

  void SetJoinSplitMixtureCoefficient(double joinSplitMixtureCoefficient) {
    mixtureCoefficient_ = joinSplitMixtureCoefficient;
    Modified();
    resetDataVisualization();
  }
  double GetJoinSplitMixtureCoefficient() {
    return mixtureCoefficient_;
  }

  void SetKeepState(bool keepState) {
    keepState_ = keepState;
    Modified();
    resetDataVisualization();
  }
  bool GetKeepState() {
    return keepState_;
  }

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
    oldEpsilonTree1 = epsilonTree1_;
    Modified();
    resetDataVisualization();
  }
  double GetEpsilonTree1() {
    return epsilonTree1_;
  }

  void SetEpsilon2Tree1(double epsilon2Tree1) {
    epsilon2Tree1_ = epsilon2Tree1;
    Modified();
    resetDataVisualization();
  }
  double GetEpsilon2Tree1() {
    return epsilon2Tree1_;
  }

  void SetEpsilon3Tree1(double epsilon3Tree1) {
    epsilon3Tree1_ = epsilon3Tree1;
    Modified();
    resetDataVisualization();
  }
  double GetEpsilon3Tree1() {
    return epsilon3Tree1_;
  }

  void SetPersistenceThreshold(double persistenceThreshold) {
    persistenceThreshold_ = persistenceThreshold;
    Modified();
    resetDataVisualization();
  }
  double GetPersistenceThreshold() {
    return persistenceThreshold_;
  }

  void SetDeleteMultiPersPairs(bool delMultiPersPairs) {
    deleteMultiPersPairs_ = delMultiPersPairs;
    Modified();
    resetDataVisualization();
  }
  bool GetDeleteMultiPersPairs() {
    return deleteMultiPersPairs_;
  }

  // Output options
  void SetComputeReconstructionError(bool b) {
    doComputeReconstructionError_ = b;
    Modified();
    resetDataVisualization();
  }
  bool GetComputeReconstructionError() {
    return doComputeReconstructionError_;
  }

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreePrincipalGeodesics *New();
  vtkTypeMacro(ttkMergeTreePrincipalGeodesics, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkMergeTreePrincipalGeodesics();

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

  template <class dataType>
  void makeBarycenterOutput(ttk::ftm::MergeTree<dataType> &barycenter,
                            int blockId,
                            vtkMultiBlockDataSet *output_barycenter);
};
