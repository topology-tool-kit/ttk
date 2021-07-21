/// \ingroup vtk
/// \class ttk::ttkFTMTreeUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#ifndef _TTKFTMTREEUTILS_H
#define _TTKFTMTREEUTILS_H

#pragma once

#include <FTMTree.h>

#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

namespace ttk {
  namespace ftm {

    void manageInconsistentArcsMultiParent(FTMTree_MT *tree) {
      ftm::idNode treeRoot = 0;
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i) {
        if(tree->getNode(i)->getNumberOfDownSuperArcs() != 0
           and tree->getNode(i)->getNumberOfUpSuperArcs() == 0)
          treeRoot = i;
      }

      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
        if(tree->getNode(i)->getNumberOfUpSuperArcs() > 1) {
          ftm::idNode lowestParent = std::numeric_limits<ftm::idNode>::max();
          for(long unsigned int j = 0;
              j < tree->getNode(i)->getNumberOfUpSuperArcs(); ++j) {
            auto tParent
              = tree->getSuperArc(tree->getNode(i)->getUpSuperArcId(j))
                  ->getUpNodeId();
            lowestParent = (lowestParent > tParent) ? tParent : lowestParent;
          }

          for(long unsigned int j = 0;
              j < tree->getNode(i)->getNumberOfUpSuperArcs(); ++j) {
            ftm::idSuperArc nodeArcId = tree->getNode(i)->getUpSuperArcId(j);
            auto tParent = tree->getSuperArc(nodeArcId)->getUpNodeId();
            if(tParent != lowestParent) {

              if(tParent == treeRoot) {
                for(long unsigned int k = 0;
                    k < tree->getNode(i)->getNumberOfDownSuperArcs(); ++k) {
                  ftm::idSuperArc nodeArcId2
                    = tree->getNode(i)->getDownSuperArcId(k);
                  auto tChildren
                    = tree->getSuperArc(nodeArcId2)->getDownNodeId();
                  if(tChildren > i) {
                    tree->getNode(i)->removeDownSuperArc(nodeArcId2);
                    tree->getNode(tChildren)->removeUpSuperArc(nodeArcId2);
                    tree->makeSuperArc(tChildren, treeRoot);
                    break;
                  }
                }
              }

              // Delete down arc from old parent
              tree->getNode(tParent)->removeDownSuperArc(nodeArcId);
              // Delete up arc from node
              tree->getNode(i)->removeUpSuperArc(nodeArcId);
            }
          }
        }
    }

    void removeSelfLink(FTMTree_MT *tree) {
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i) {
        for(unsigned int j = 0; j < tree->getNode(i)->getNumberOfUpSuperArcs();
            ++j) {
          ftm::idSuperArc nodeArcId = tree->getNode(i)->getUpSuperArcId(j);
          auto tParent = tree->getSuperArc(nodeArcId)->getUpNodeId();
          if(tParent == i) {
            // Delete down arc
            tree->getNode(i)->removeDownSuperArc(nodeArcId);
            // Delete up arc
            tree->getNode(i)->removeUpSuperArc(nodeArcId);
          }
        }
      }
    }

    template <class dataType>
    MergeTree<dataType> makeTree(vtkUnstructuredGrid *treeNodes,
                                 vtkUnstructuredGrid *treeArcs) {
      // Init Scalars
      // std::cout << "// Init Scalars" << std::endl;
      Scalars scalars;
      vtkSmartPointer<vtkDataArray> nodesScalar
        = treeNodes->GetPointData()->GetArray("Scalar"); // 1: Scalar
      scalars.size = nodesScalar->GetNumberOfTuples();
      std::vector<dataType> scalarsValues;
      for(int i = 0; i < nodesScalar->GetNumberOfTuples(); ++i)
        scalarsValues.push_back(nodesScalar->GetTuple1(i));
      // scalars.values = ttkUtils::GetVoidPointer(nodesScalar);
      scalars.values = (void *)(scalarsValues.data());

      // Init Tree
      // std::cout << "// Init tree" << std::endl;
      Params params;
      params.treeType = Join_Split;
      // MergeTree<dataType> mergeTree(scalars, params);
      MergeTree<dataType> mergeTree(scalars, scalarsValues, params);

      // Add Nodes
      // std::cout << "// Add nodes" << std::endl;
      vtkSmartPointer<vtkDataArray> nodesId
        = treeNodes->GetPointData()->GetArray("NodeId"); // 0: NodeId
      vtkIdType nodesNumTuples = nodesId->GetNumberOfTuples();
      for(vtkIdType i = 0; i < nodesNumTuples; ++i) {
        mergeTree.tree.makeNode(i);
      }

      // Add Arcs
      // std::cout << "// Add arcs" << std::endl;
      vtkSmartPointer<vtkDataArray> arcsUp
        = treeArcs->GetCellData()->GetArray("upNodeId"); // 1: upNodeId
      vtkSmartPointer<vtkDataArray> arcsDown
        = treeArcs->GetCellData()->GetArray("downNodeId"); // 2: downNodeId
      vtkIdType arcsNumTuples = arcsUp->GetNumberOfTuples();
      vtkSmartPointer<vtkDataArray> dummyArcArray
        = treeArcs->GetCellData()->GetArray("isDummyArc");
      std::set<std::tuple<double, double>> added_arcs; // Avoid duplicates
      for(vtkIdType i = 0; i < arcsNumTuples; ++i) {
        if(dummyArcArray != nullptr and dummyArcArray->GetTuple1(i) == 1)
          continue;
        double downId = arcsDown->GetTuple1(i);
        double upId = arcsUp->GetTuple1(i);
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

    void loadBlocks(std::vector<vtkMultiBlockDataSet *> &inputTrees,
                    vtkMultiBlockDataSet *blocks) {
      if(blocks != nullptr) {
        inputTrees.resize(blocks->GetNumberOfBlocks());
        for(size_t i = 0; i < inputTrees.size(); ++i) {
          inputTrees[i]
            = vtkMultiBlockDataSet::SafeDownCast(blocks->GetBlock(i));
        }
      }
    }

    template <class dataType>
    void constructTrees(std::vector<vtkMultiBlockDataSet *> &inputTrees,
                        std::vector<MergeTree<dataType>> &intermediateTrees,
                        std::vector<vtkUnstructuredGrid *> &treesNodes,
                        std::vector<vtkUnstructuredGrid *> &treesArcs,
                        std::vector<vtkDataSet *> &treesSegmentation) {
      const int numInputs = inputTrees.size();
      intermediateTrees = std::vector<MergeTree<dataType>>(numInputs);
      treesNodes = std::vector<vtkUnstructuredGrid *>(numInputs);
      treesArcs = std::vector<vtkUnstructuredGrid *>(numInputs);
      treesSegmentation = std::vector<vtkDataSet *>(numInputs);
      for(int i = 0; i < numInputs; i++) {
        treesNodes[i]
          = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(0));
        treesArcs[i]
          = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(1));
        if(inputTrees[i]->GetNumberOfBlocks() > 2)
          treesSegmentation[i]
            = vtkDataSet::SafeDownCast(inputTrees[i]->GetBlock(2));
        intermediateTrees[i] = makeTree<dataType>(treesNodes[i], treesArcs[i]);
      }
    }

    template <class dataType>
    void constructTrees(std::vector<vtkMultiBlockDataSet *> &inputTrees,
                        std::vector<MergeTree<dataType>> &intermediateTrees) {
      std::vector<vtkUnstructuredGrid *> treesNodes;
      std::vector<vtkUnstructuredGrid *> treesArcs;
      std::vector<vtkDataSet *> treesSegmentation;
      constructTrees(inputTrees, intermediateTrees, treesNodes, treesArcs,
                     treesSegmentation);
    }
  } // namespace ftm
} // namespace ttk

#endif
