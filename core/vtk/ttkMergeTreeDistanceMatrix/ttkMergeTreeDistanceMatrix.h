/// \ingroup vtk
/// \class ttkMergeTreeDistanceMatrix
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2021.
///
/// \brief TTK VTK-filter that wraps the ttk::MergeTreeDistanceMatrix module.
///
/// This VTK filter uses the ttk::MergeTreeDistanceMatrix module to compute the
/// distance matrix of a group of merge trees.
///
/// \param Input vtkMultiBlockDataset
/// \param Output vtkTable
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MergeTreeDistanceMatrix
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeClustering/">Merge
///   Tree Clustering example</a> \n

#pragma once

// VTK Module
#include <ttkMergeTreeDistanceMatrixModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>

// TTK Base Includes
#include <MergeTreeDistanceMatrix.h>

class TTKMERGETREEDISTANCEMATRIX_EXPORT ttkMergeTreeDistanceMatrix
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::MergeTreeDistanceMatrix // and we inherit from the base class
{
private:
  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */
  // Execution Options
  int Backend = 0;

  bool UseFieldDataParameters = false;

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */
  // Input Options
  void SetEpsilon1UseFarthestSaddle(bool epsilon1UseFarthestSaddle) {
    epsilon1UseFarthestSaddle_ = epsilon1UseFarthestSaddle;
    Modified();
  }
  bool GetEpsilon1UseFarthestSaddle() {
    return epsilon1UseFarthestSaddle_;
  }

  void SetEpsilonTree1(double epsilonTree1) {
    epsilonTree1_ = epsilonTree1;
    Modified();
  }
  double SetEpsilonTree1() {
    return epsilonTree1_;
  }

  void SetEpsilon2Tree1(double epsilon2Tree1) {
    epsilon2Tree1_ = epsilon2Tree1;
    Modified();
  }
  double SetEpsilon2Tree1() {
    return epsilon2Tree1_;
  }

  void SetEpsilon3Tree1(double epsilon3Tree1) {
    epsilon3Tree1_ = epsilon3Tree1;
    Modified();
  }
  double SetEpsilon3Tree1() {
    return epsilon3Tree1_;
  }

  void SetPersistenceThreshold(double persistenceThreshold) {
    persistenceThreshold_ = persistenceThreshold;
    Modified();
  }
  double SetPersistenceThreshold() {
    return persistenceThreshold_;
  }

  void SetDeleteMultiPersPairs(bool doDelete) {
    deleteMultiPersPairs_ = doDelete;
    Modified();
  }
  bool SetDeleteMultiPersPairs() {
    return deleteMultiPersPairs_;
  }

  // Execution Options
  vtkSetMacro(Backend, int);
  vtkGetMacro(Backend, int);

  void SetAssignmentSolver(int assignmentSolver) {
    assignmentSolverID_ = assignmentSolver;
    Modified();
  }
  int GetAssignmentSolver() {
    return assignmentSolverID_;
  }

  void SetBranchDecomposition(bool branchDecomposition) {
    branchDecomposition_ = branchDecomposition;
    Modified();
  }
  int GetBranchDecomposition() {
    return branchDecomposition_;
  }

  void SetNormalizedWasserstein(bool normalizedWasserstein) {
    normalizedWasserstein_ = normalizedWasserstein;
    Modified();
  }
  int GetNormalizedWasserstein() {
    return normalizedWasserstein_;
  }

  void SetKeepSubtree(bool keepSubtree) {
    keepSubtree_ = keepSubtree;
    Modified();
  }
  int GetKeepSubtree() {
    return keepSubtree_;
  }

  void SetDistanceSquared(bool distanceSquared) {
    distanceSquared_ = distanceSquared;
    Modified();
  }
  int GetDistanceSquared() {
    return distanceSquared_;
  }

  vtkSetMacro(UseFieldDataParameters, bool);
  vtkGetMacro(UseFieldDataParameters, bool);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreeDistanceMatrix *New();
  vtkTypeMacro(ttkMergeTreeDistanceMatrix, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   *         (see cpp file)
   */
  ttkMergeTreeDistanceMatrix();
  ~ttkMergeTreeDistanceMatrix() override;

  /**
   * Specify the input data type of each input port
   *         (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * Specify the data object type of each output port
   *         (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * Pass VTK data to the base code and convert base code output to VTK
   *          (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <class dataType>
  int run(vtkInformationVector *outputVector,
          std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
          std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2);
};
