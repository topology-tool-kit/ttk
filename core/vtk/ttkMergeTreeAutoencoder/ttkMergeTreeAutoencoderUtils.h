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
    /**
     * @brief Proxy function to use ttkMergeTreeVisualization to create the vtk
     * objects of a merge tree.
     *
     * @param[in] tree input tree to process.
     * @param[in] treeNodes original vtk nodes object of the tree.
     * @param[in] treeNodeCorr correspondence between the input tree and the vtk
     * object.
     * @param[in] treeSegmentation original segmentation.
     * @param[out] vtkOutputNode output vtk object: node.
     * @param[out] vtkOutputArc output vtk object: arc.
     * @param[out] vtkOutputSegmentation output vtk object: segmentation.
     * @param[in] treeID identifier of the tree.
     * @param[in] customIntArrays custom arrays to add on nodes.
     * @param[in] customDoubleArrays custom arrays to add on nodes.
     * @param[in] outputSegmentation choose to output segmentation.
     * @param[in] mixtureCoefficient weight of pairs type.
     * @param[in] isPersistenceDiagram if input is a diagram.
     * @param[in] convertToDiagram if input is a diagram but the vtk object is a
     * tree.
     * @param[in] debugLevel debug level.
     */
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

    /**
     * @brief Proxy function to use ttkMergeTreeVisualization to create the vtk
     * objects of a set of merge trees.
     *
     * @param[in] trees input trees to process.
     * @param[in] treesNodesT original vtk nodes object of the trees.
     * @param[in] treesNodeCorr correspondence between the input trees and the
     * vtk objects.
     * @param[in] treesSegmentationT original segmentations.
     * @param[out] output output blocks.
     * @param[in] customIntArrays custom arrays to add on nodes.
     * @param[in] customDoubleArrays custom arrays to add on nodes.
     * @param[in] mixtureCoefficient weight of pairs type.
     * @param[in] isPersistenceDiagram if input is a diagram.
     * @param[in] convertToDiagram if input is a diagram but the vtk object is a
     * tree.
     * @param[in] debugLevel debug level.
     */
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

    /**
     * @brief Proxy function to use ttkMergeTreeVisualization to create the vtk
     * objects of a set of merge trees.
     *
     * @param[in] trees input trees to process.
     * @param[in] treesNodesT original vtk nodes object of the trees.
     * @param[in] treesNodeCorr correspondence between the input trees and the
     * vtk objects.
     * @param[out] output output blocks.
     * @param[in] customIntArrays custom arrays to add on nodes.
     * @param[in] customDoubleArrays custom arrays to add on nodes.
     * @param[in] mixtureCoefficient weight of pairs type.
     * @param[in] isPersistenceDiagram if input is a diagram.
     * @param[in] convertToDiagram if input is a diagram but the vtk object is a
     * tree.
     * @param[in] debugLevel debug level.
     */
    void makeManyOutput(
      std::vector<ttk::ftm::MergeTree<float> *> &trees,
      std::vector<vtkUnstructuredGrid *> &treesNodesT,
      std::vector<std::vector<int>> &treesNodeCorr,
      vtkSmartPointer<vtkMultiBlockDataSet> &output,
      std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
        &customIntArrays,
      std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
        &customDoubleArrays,
      double mixtureCoefficient,
      bool isPersistenceDiagram,
      bool convertToDiagram,
      int debugLevel);

    /**
     * @brief Proxy function to use ttkMergeTreeVisualization to create the vtk
     * objects of a set of merge trees.
     *
     * @param[in] trees input trees to process.
     * @param[out] output output blocks.
     * @param[in] customIntArrays custom arrays to add on nodes.
     * @param[in] customDoubleArrays custom arrays to add on nodes.
     * @param[in] mixtureCoefficient weight of pairs type.
     * @param[in] isPersistenceDiagram if input is a diagram.
     * @param[in] convertToDiagram if input is a diagram but the vtk object is a
     * tree.
     * @param[in] debugLevel debug level.
     */
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

    /**
     * @brief Proxy function to use ttkMergeTreeVisualization to create the vtk
     * objects of a set of merge trees.
     *
     * @param[in] trees input trees to process.
     * @param[out] output output blocks.
     * @param[in] mixtureCoefficient weight of pairs type.
     * @param[in] isPersistenceDiagram if input is a diagram.
     * @param[in] convertToDiagram if input is a diagram but the vtk object is a
     * tree.
     * @param[in] debugLevel debug level.
     */
    void makeManyOutput(std::vector<ttk::ftm::MergeTree<float> *> &trees,
                        vtkSmartPointer<vtkMultiBlockDataSet> &output,
                        double mixtureCoefficient,
                        bool isPersistenceDiagram,
                        bool convertToDiagram,
                        int debugLevel);

    /**
     * @brief Compute information between the origins given their matching, such
     * as the latent importance.
     *
     * @param[in] origins origins of each input basis.
     * @param[in] originsPrime origins of each output basis.
     * @param[in] originsMatchingVectorT matchings from the origin of layer l
     * to the one of layer l+1.
     * @param[in] invOriginsMatchingVectorT matchings from the origin of
     * layer l+1 to the one of layer l.
     * @param[in] isPersistenceDiagram if input is a diagram.
     * @param[out] originsMatchingVector matchings from the origin of layer l
     * to the one of layer 0.
     * @param[out] originsPersPercent percentage of the persistence for each
     * pair of the origin of each layer compared to its matched pair in the
     * origin of first layer.
     * @param[out] originsPersDiff difference of the persistence for each
     * pair of the origin of each layer compared to its matched pair in the
     * origin of first layer.
     * @param[out] originPersPercent percentage of each pair in the origin of
     * the first layer compared to its matched pair of the origin in the latent
     * layer.
     * @param[out] originPersDiff difference of each pair in the origin of
     * the first layer compared to its matched pair of the origin in the latent
     * layer.
     * @param[out] originPersistenceOrder array assigning for each node the
     * index of its persistence pairs (ordered by persistence).
     */
    void computeTrackingInformation(
      std::vector<mtu::TorchMergeTree<float>> &origins,
      std::vector<mtu::TorchMergeTree<float>> &originsPrime,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVectorT,
      std::vector<std::vector<ttk::ftm::idNode>> &invOriginsMatchingVectorT,
      bool isPersistenceDiagram,
      std::vector<std::vector<ttk::ftm::idNode>> &originsMatchingVector,
      std::vector<std::vector<double>> &originsPersPercent,
      std::vector<std::vector<double>> &originsPersDiff,
      std::vector<double> &originPersPercent,
      std::vector<double> &originPersDiff,
      std::vector<int> &originPersistenceOrder);

    /**
     * @brief Compute the custom arrays to add to the vtk outputs.
     *
     * @param[in] recs merge trees at each layer.
     * @param[in] persCorrelationMatrix correlation matrix.
     * @param[in] invDataMatchingVectorT matchings from trees of layer l to the
     * origin of layer l.
     * @param[in] invReconstMatchingVectorT matchings from output trees to input
     * trees.
     * @param[in] originsMatchingVector matchings from the origin of layer l to
     * the one of layer 0.
     * @param[in] originsMatchingVectorT matchings from the origin of layer l
     * to the one of layer l+1.
     * @param[in] originsPersPercent percentage of the persistence for each
     * pair of the origin of each layer compared to its matched pair in the
     * origin of first layer.
     * @param[in] originPersDiff difference of each pair in the origin of
     * the first layer compared to its matched pair of the origin in the latent
     * layer.
     * @param[in] originPersistenceOrder array assigning for each node the
     * index of its persistence pairs (ordered by persistence).
     * @param[in] l layer index.
     * @param[in] lShift layer index shift.
     * @param[out] customIntArrays custom arrays to add on nodes.
     * @param[out] customDoubleArrays custom arrays to add on nodes.
     */
    void computeCustomArrays(
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &recs,
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
