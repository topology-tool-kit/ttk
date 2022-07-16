/// \ingroup vtk
/// \class ttk::ttkFTMTreeUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#pragma once

#include <FTMTree.h>
#include <FTMTreeUtils.h>

#include <ttkUtils.h>

#include <vtkCellData.h>
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
      Scalars scalars;
      vtkSmartPointer<vtkDataArray> nodesScalar
        = treeNodes->GetPointData()->GetArray("Scalar"); // 1: Scalar
      scalars.size = nodesScalar->GetNumberOfTuples();
      std::vector<dataType> scalarsValues(nodesScalar->GetNumberOfTuples());
      for(int i = 0; i < nodesScalar->GetNumberOfTuples(); ++i) {
        int index = (treeNodeIdArray ? treeNodeIdArray->GetTuple1(i) : i);
        scalarsValues[index] = nodesScalar->GetTuple1(i);
      }
      scalars.values = (void *)(scalarsValues.data());

      // Init Tree
      Params params;
      params.treeType = Join_Split;
      MergeTree<dataType> mergeTree(scalars, scalarsValues, params);

      // Add Nodes
      vtkSmartPointer<vtkDataArray> nodesId
        = treeNodes->GetPointData()->GetArray("NodeId"); // 0: NodeId
      vtkIdType nodesNumTuples = nodesId->GetNumberOfTuples();
      for(vtkIdType i = 0; i < nodesNumTuples; ++i) {
        mergeTree.tree.makeNode(i);
      }

      // Add Arcs
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

    inline void
      loadBlocks(std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
                 vtkMultiBlockDataSet *blocks) {
      if(blocks != nullptr) {
        inputTrees.resize(
          vtkMultiBlockDataSet::SafeDownCast(blocks->GetBlock(0))
            ->GetNumberOfBlocks());
        for(size_t i = 0; i < inputTrees.size(); ++i) {
          vtkSmartPointer<vtkMultiBlockDataSet> vtkBlock
            = vtkSmartPointer<vtkMultiBlockDataSet>::New();
          vtkBlock->SetNumberOfBlocks(blocks->GetNumberOfBlocks());
          for(unsigned int j = 0; j < blocks->GetNumberOfBlocks(); ++j)
            vtkBlock->SetBlock(
              j, vtkMultiBlockDataSet::SafeDownCast(blocks->GetBlock(j))
                   ->GetBlock(i));
          inputTrees[i] = vtkBlock;
        }
      }
    }

    template <class dataType>
    void constructTrees(
      std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
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
    void constructTrees(
      std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
      std::vector<MergeTree<dataType>> &intermediateTrees) {
      std::vector<vtkUnstructuredGrid *> treesNodes;
      std::vector<vtkUnstructuredGrid *> treesArcs;
      std::vector<vtkDataSet *> treesSegmentation;
      constructTrees(inputTrees, intermediateTrees, treesNodes, treesArcs,
                     treesSegmentation);
    }
  } // namespace ftm
} // namespace ttk
