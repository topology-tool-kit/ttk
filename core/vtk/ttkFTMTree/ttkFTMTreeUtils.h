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
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

using namespace ttk;
using namespace ftm;

template <typename dataType>
struct MergeTree {
  ftm::Scalars scalars;
  std::vector<dataType> scalarsValues;
  ftm::Params params;
  ftm::FTMTree_MT tree;

  ftm::Scalars emptyScalars() {
    ftm::Scalars scalarsT;
    scalarsT.size = 0;
    dataType *scalarsValuesT = nullptr;
    scalarsT.values = (void *)scalarsValuesT;
    return scalarsT;
  }

  ftm::Params emptyParams() {
    ftm::Params paramsT;
    paramsT.treeType = ftm::Join_Split;
    return paramsT;
  }

  MergeTree() : MergeTree(emptyScalars(), emptyParams()) {
  }

  MergeTree(ftm::Scalars scalarsT, ftm::Params paramsT)
    : scalars(scalarsT), params(paramsT),
      tree(&params, &scalars, params.treeType) {
    tree.makeAlloc();
  }

  MergeTree(ftm::Scalars scalarsT,
            std::vector<dataType> scalarValuesT,
            ftm::Params paramsT)
    : scalars(scalarsT), scalarsValues(scalarValuesT), params(paramsT),
      tree(&params, &scalars, params.treeType) {
    scalars.values = (void *)(scalarsValues.data());
    tree.makeAlloc();
  }
};

template <class dataType>
MergeTree<dataType> makeTree(vtkUnstructuredGrid *treeNodes,
                             vtkUnstructuredGrid *treeArcs) {
  // Init Scalars
  Scalars scalars;
  vtkSmartPointer<vtkDataArray> nodesScalar
    = treeNodes->GetPointData()->GetArray("Scalar"); // 1: Scalar
  scalars.size = nodesScalar->GetNumberOfTuples();
  scalars.values = ttkUtils::GetVoidPointer(nodesScalar);

  // Init Tree
  Params params;
  params.treeType = Join_Split;
  MergeTree<dataType> mergeTree(scalars, params);
  mergeTree.tree.makeAlloc();

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
  std::set<std::tuple<double, double>> added_arcs; // Avoid duplicates
  for(vtkIdType i = 0; i < arcsNumTuples; ++i) {
    double downId = arcsDown->GetTuple1(i);
    double upId = arcsUp->GetTuple1(i);
    auto it = added_arcs.find(std::make_tuple(downId, upId));
    if(it == added_arcs.end()) { // arc not added yet
      mergeTree.tree.makeSuperArc(downId, upId); // (down, Up)
      added_arcs.insert(std::make_tuple(downId, upId));
    }
  }

  // Manage inconsistent arcs
  // manageInconsistentArcsMultiParent(tree);;

  // Remove self link
  // removeSelfLink(tree);

  // tree->printTree2();

  return mergeTree;
}

#endif
