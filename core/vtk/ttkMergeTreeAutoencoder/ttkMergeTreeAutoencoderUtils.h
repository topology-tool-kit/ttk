/// \ingroup vtk
/// \class ttk::ttkMergeTreeAutoencoderUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2023.
///
/// Utils function for manipulating merge and contour tree classes

#pragma once

#include <MergeTreeTorchUtils.h>
#include <ttkMergeTreeVisualization.h>
#include <ttkUtils.h>

#include <vtkMultiBlockDataSet.h>

#ifdef TTK_ENABLE_TORCH
namespace ttk {
  namespace wae {
    void makeOneOutput(
      ttk::ftm::MergeTree<float> &tree,
      vtkUnstructuredGrid *treeNodes,
      std::vector<int> &treeNodeCorr,
      vtkDataSet *treeSegmentation,
      vtkSmartPointer<vtkUnstructuredGrid> &vtkOutputNode,
      vtkSmartPointer<vtkUnstructuredGrid> &vtkOutputArc,
      vtkSmartPointer<vtkDataSet> &vtkOutputSegmentation,
      unsigned int treeID,
      std::vector<std::tuple<std::string, std::vector<int>>> &customIntArrays,
      std::vector<std::tuple<std::string, std::vector<double>>>
        &customDoubleArrays,
      bool outputSegmentation,
      double mixtureCoefficient,
      bool isPersistenceDiagram,
      bool convertToDiagram,
      int debugLevel);

    void makeManyOutput(
      std::vector<ttk::ftm::MergeTree<float> *> &trees,
      std::vector<vtkUnstructuredGrid *> &treesNodesT,
      std::vector<std::vector<int>> &treesNodeCorr,
      std::vector<vtkDataSet *> &treesSegmentationT,
      vtkSmartPointer<vtkMultiBlockDataSet> &output,
      std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
        &customIntArrays,
      std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
        &customDoubleArrays,
      double mixtureCoefficient,
      bool isPersistenceDiagram,
      bool convertToDiagram,
      int debugLevel);

    void makeManyOutput(
      std::vector<ttk::ftm::MergeTree<float> *> &trees,
      vtkSmartPointer<vtkMultiBlockDataSet> &output,
      std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
        &customIntArrays,
      std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
        &customDoubleArrays,
      double mixtureCoefficient,
      bool isPersistenceDiagram,
      bool convertToDiagram,
      int debugLevel);

    void makeManyOutput(std::vector<ttk::ftm::MergeTree<float> *> &trees,
                        vtkSmartPointer<vtkMultiBlockDataSet> &output,
                        double mixtureCoefficient,
                        bool isPersistenceDiagram,
                        bool convertToDiagram,
                        int debugLevel);

    void computeTrackingInformation(
      std::vector<MergeTreeTorchUtils::TorchMergeTree<float>> &origins,
      std::vector<MergeTreeTorchUtils::TorchMergeTree<float>> &originsPrime,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &invOriginsMatchingVectorT,
      bool isPersistenceDiagram,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVector,
      std::vector<std::vector<double>> &originsPersPercent,
      std::vector<std::vector<double>> &originsPersDiff,
      std::vector<double> &originPersPercent,
      std::vector<double> &originPersDiff,
      std::vector<int> &originPersistenceOrder);

    void computeCustomArrays(
      std::vector<std::vector<MergeTreeTorchUtils::TorchMergeTree<float>>>
        &recs,
      std::vector<std::vector<double>> &persCorrelationMatrix,
      std::vector<std::vector<std::vector<ttk::ftm::idNode>>>
        &invDataMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &invReconstMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVector,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVectorT,
      std::vector<std::vector<double>> &originsPersPercent,
      std::vector<std::vector<double>> &originPersDiff,
      std::vector<int> &originPersistenceOrder,
      unsigned int l,
      unsigned int lShift,
      std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
        &customIntArrays,
      std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
        &customDoubleArrays);
  } // namespace wae
} // namespace ttk
#endif
