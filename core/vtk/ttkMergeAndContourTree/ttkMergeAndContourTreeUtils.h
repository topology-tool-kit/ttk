/// \ingroup vtk
/// \class ttk::ttkMergeAndContourTreeUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating MergeAndContourTree class

#pragma once

#include <FTMTree.h>
#include <FTMTreeUtils.h>

#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

namespace ttk {
  namespace ftm {

    template <class dataType>
    MergeTree<dataType> makeTree(vtkUnstructuredGrid *treeNodes,
                                 vtkUnstructuredGrid *treeArcs) {
      auto treeNodeIdArray = treeNodes->GetPointData()->GetArray("TreeNodeId");

      // Init Scalars
      auto scalars = std::make_shared<Scalars>();
      vtkSmartPointer<vtkDataArray> const nodesScalar
        = treeNodes->GetPointData()->GetArray("Scalar"); // 1: Scalar
      scalars->size = nodesScalar->GetNumberOfTuples();
      auto scalarsValues = std::make_shared<std::vector<dataType>>(
        nodesScalar->GetNumberOfTuples());
      for(int i = 0; i < nodesScalar->GetNumberOfTuples(); ++i) {
        int const index = (treeNodeIdArray ? treeNodeIdArray->GetTuple1(i) : i);
        (*scalarsValues)[index] = nodesScalar->GetTuple1(i);
      }
      scalars->values = (void *)(scalarsValues->data());

      // Init Tree
      auto params = std::make_shared<Params>();
      params->treeType = Join_Split;
      MergeTree<dataType> mergeTree(scalars, scalarsValues, params);

      // Add Nodes
      vtkSmartPointer<vtkDataArray> const nodesId
        = treeNodes->GetPointData()->GetArray("NodeId"); // 0: NodeId
      vtkIdType const nodesNumTuples = nodesId->GetNumberOfTuples();
      for(vtkIdType i = 0; i < nodesNumTuples; ++i) {
        mergeTree.tree.makeNode(i);
      }

      // Add Arcs
      vtkSmartPointer<vtkDataArray> const arcsUp
        = treeArcs->GetCellData()->GetArray("upNodeId"); // 1: upNodeId
      vtkSmartPointer<vtkDataArray> const arcsDown
        = treeArcs->GetCellData()->GetArray("downNodeId"); // 2: downNodeId
      vtkIdType const arcsNumTuples = arcsUp->GetNumberOfTuples();
      vtkSmartPointer<vtkDataArray> const dummyArcArray
        = treeArcs->GetCellData()->GetArray("isDummyArc");
      std::set<std::tuple<double, double>> added_arcs; // Avoid duplicates
      for(vtkIdType i = 0; i < arcsNumTuples; ++i) {
        if(dummyArcArray != nullptr and dummyArcArray->GetTuple1(i) == 1)
          continue;
        double downId = arcsDown->GetTuple1(i);
        double upId = arcsUp->GetTuple1(i);
        if(treeNodeIdArray) {
          downId = treeNodeIdArray->GetTuple1(downId);
          upId = treeNodeIdArray->GetTuple1(upId);
        }
        auto it = added_arcs.find(std::make_tuple(downId, upId));
        if(it == added_arcs.end()) { // arc not added yet
          mergeTree.tree.makeSuperArc(downId, upId); // (down, Up)
          added_arcs.insert(std::make_tuple(downId, upId));
        }
      }

      // Manage inconsistent arcs
      manageInconsistentArcsMultiParent(&(mergeTree.tree));

      // Remove self link
      removeSelfLink(&(mergeTree.tree));

      return mergeTree;
    }

    // Returns a branch decomposition tree
    template <class dataType>
    MergeTree<dataType>
      makeBDTreeFromPDGrid(vtkUnstructuredGrid *persistenceDiagram,
                           bool useSadMaxPairs = true) {
      auto birthArray
        = persistenceDiagram->GetCellData()->GetArray(PersistenceBirthName);
      auto persArray
        = persistenceDiagram->GetCellData()->GetArray(PersistenceName);
      auto pairTypeArray
        = persistenceDiagram->GetCellData()->GetArray(PersistencePairTypeName);
      auto criticalTypeArray = persistenceDiagram->GetPointData()->GetArray(
        PersistenceCriticalTypeName);

      auto treeNodeIdArray
        = persistenceDiagram->GetPointData()->GetArray("TreeNodeId");
      auto treeNodeId2Array
        = persistenceDiagram->GetPointData()->GetArray("TreeNodeIdOrigin");
      bool const gotNodeArrays = (treeNodeIdArray and treeNodeId2Array);

      auto noPairs = birthArray->GetNumberOfTuples();
      int noNodes = noPairs * 2;
      if(gotNodeArrays) {
        for(vtkIdType i = 0; i < treeNodeIdArray->GetNumberOfTuples(); ++i) {
          int const val = std::max(treeNodeIdArray->GetTuple1(i),
                                   treeNodeId2Array->GetTuple1(i))
                          + 1;
          noNodes = std::max(noNodes, val);
        }
      }
      std::vector<dataType> scalarsVector(noNodes);

      // Init Tree
      MergeTree<dataType> mergeTree
        = ttk::ftm::createEmptyMergeTree<dataType>(scalarsVector.size());

      // Init critical type enum values
      auto locMin = static_cast<int>(CriticalType::Local_minimum);
      auto locMax = static_cast<int>(CriticalType::Local_maximum);

      // Get Min-Max pair index
      int minMaxPairIndex = -1;
      for(vtkIdType i = 0; i < noPairs; ++i) {
        vtkIdType npts;
        vtkIdType const *pts;
        persistenceDiagram->GetCellPoints(i, npts, pts);
        auto ct1 = criticalTypeArray->GetTuple1(pts[0]);
        auto ct2 = criticalTypeArray->GetTuple1(pts[1]);
        if((ct1 == locMin and ct2 == locMax)
           or (ct1 == locMax and ct2 == locMin))
          minMaxPairIndex = i;
      }

      // Init Tree Structure
      for(vtkIdType i = 0; i < noNodes; ++i)
        mergeTree.tree.makeNode(i);
      for(vtkIdType i = 0; i < noPairs; ++i) {
        auto pairType = pairTypeArray->GetTuple1(i);
        vtkIdType npts;
        vtkIdType const *pts;
        persistenceDiagram->GetCellPoints(i, npts, pts);
        auto ct1 = criticalTypeArray->GetTuple1(pts[0]);
        auto ct2 = criticalTypeArray->GetTuple1(pts[1]);
        if((pairType == -1
            or (useSadMaxPairs and ct1 != locMax and ct2 != locMax)
            or (not useSadMaxPairs and ct1 != locMin and ct2 != locMin))
           and i != minMaxPairIndex)
          continue;
        int const index1
          = (gotNodeArrays ? treeNodeId2Array->GetTuple1(pts[0]) : i * 2);
        int const index2
          = (gotNodeArrays ? treeNodeIdArray->GetTuple1(pts[0]) : i * 2 + 1);
        mergeTree.tree.getNode(index1)->setOrigin(index2);
        mergeTree.tree.getNode(index2)->setOrigin(index1);
        scalarsVector[index1] = birthArray->GetTuple1(i);
        scalarsVector[index2]
          = birthArray->GetTuple1(i) + persArray->GetTuple1(i);

        if(i != minMaxPairIndex) {
          auto up = (gotNodeArrays ? treeNodeIdArray->GetTuple1(minMaxPairIndex)
                                   : minMaxPairIndex * 2 + 1);
          // mergeTree.tree.makeSuperArc(index1, up);
          mergeTree.tree.makeSuperArc(index2, up);
        }
      }

      // Init scalars
      ttk::ftm::setTreeScalars<dataType>(mergeTree, scalarsVector);

      return mergeTree;
    }

    inline void
      loadBlocks(std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
                 vtkMultiBlockDataSet *blocks) {
      if(blocks != nullptr) {
        if(blocks->GetBlock(0)->IsA("vtkMultiBlockDataSet"))
          inputTrees.resize(
            vtkMultiBlockDataSet::SafeDownCast(blocks->GetBlock(0))
              ->GetNumberOfBlocks());
        else if(blocks->GetBlock(0)->IsA("vtkUnstructuredGrid"))
          inputTrees.resize(blocks->GetNumberOfBlocks());
        for(size_t i = 0; i < inputTrees.size(); ++i) {
          if(blocks->GetBlock(0)->IsA("vtkMultiBlockDataSet")) {
            vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlock
              = vtkSmartPointer<vtkMultiBlockDataSet>::New();
            vtkBlock->SetNumberOfBlocks(blocks->GetNumberOfBlocks());
            for(unsigned int j = 0; j < blocks->GetNumberOfBlocks(); ++j)
              vtkBlock->SetBlock(
                j, vtkMultiBlockDataSet::SafeDownCast(blocks->GetBlock(j))
                     ->GetBlock(i));
            inputTrees[i] = vtkBlock;
          } else if(blocks->GetBlock(0)->IsA("vtkUnstructuredGrid")) {
            vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlock
              = vtkSmartPointer<vtkMultiBlockDataSet>::New();
            vtkBlock->SetNumberOfBlocks(1);
            vtkBlock->SetBlock(
              0, vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i)));
            inputTrees[i] = vtkBlock;
          }
        }
      }
    }

    template <class dataType>
    bool constructTrees(
      std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
      std::vector<MergeTree<dataType>> &intermediateTrees,
      std::vector<vtkUnstructuredGrid *> &treesNodes,
      std::vector<vtkUnstructuredGrid *> &treesArcs,
      std::vector<vtkDataSet *> &treesSegmentation,
      std::vector<bool> useSadMaxPairs) {
      bool isPersistenceDiagram = false;
      const int numInputs = inputTrees.size();
      intermediateTrees.resize(numInputs);
      treesNodes.resize(numInputs);
      treesArcs.resize(numInputs);
      treesSegmentation.resize(numInputs);
      for(int i = 0; i < numInputs; i++) {
        if(inputTrees[i]->GetNumberOfBlocks() >= 2) {
          treesNodes[i]
            = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(0));
          treesArcs[i]
            = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(1));
          if(inputTrees[i]->GetNumberOfBlocks() > 2)
            treesSegmentation[i]
              = vtkDataSet::SafeDownCast(inputTrees[i]->GetBlock(2));
          intermediateTrees[i]
            = makeTree<dataType>(treesNodes[i], treesArcs[i]);
        } else {
          treesNodes[i]
            = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(0));
          vtkUnstructuredGrid *persistenceDiagram
            = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(0));
          intermediateTrees[i] = makeBDTreeFromPDGrid<dataType>(
            persistenceDiagram, useSadMaxPairs[i]);
          isPersistenceDiagram = true;
        }
      }
      return isPersistenceDiagram;
    }

    template <class dataType>
    bool constructTrees(
      std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
      std::vector<MergeTree<dataType>> &intermediateTrees,
      std::vector<vtkUnstructuredGrid *> &treesNodes,
      std::vector<vtkUnstructuredGrid *> &treesArcs,
      std::vector<vtkDataSet *> &treesSegmentation,
      bool useSadMaxPairs = true) {
      std::vector<bool> const useSadMaxPairsVec(
        inputTrees.size(), useSadMaxPairs);
      return constructTrees(inputTrees, intermediateTrees, treesNodes,
                            treesArcs, treesSegmentation, useSadMaxPairsVec);
    }

    template <class dataType>
    bool constructTrees(
      std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
      std::vector<MergeTree<dataType>> &intermediateTrees,
      bool useSadMaxPairs = true) {
      std::vector<vtkUnstructuredGrid *> treesNodes;
      std::vector<vtkUnstructuredGrid *> treesArcs;
      std::vector<vtkDataSet *> treesSegmentation;
      return constructTrees(inputTrees, intermediateTrees, treesNodes,
                            treesArcs, treesSegmentation, useSadMaxPairs);
    }
  } // namespace ftm
} // namespace ttk
