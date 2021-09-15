/// \ingroup base
/// \class vtk::ttkMergeTreeTemporalReductionUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///

#ifndef _TTKMERGETREETEMPORALREDUCTIONUTILS_H
#define _TTKMERGETREETEMPORALREDUCTIONUTILS_H

#include <FTMTreeUtils.h>

#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

using namespace ttk;
using namespace ftm;

// -------------------------------------------------------------------------------------
// Temporal Subsampling
// -------------------------------------------------------------------------------------
template <class dataType>
void makeTemporalSubsamplingOutput(
  std::vector<MergeTree<dataType>> &intermediateMTrees,
  std::vector<std::vector<double>> &embedding,
  std::vector<MergeTree<dataType>> &allMT,
  std::vector<int> removed,
  vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1,
  vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1,
  bool displayRemoved = false) {
  auto fullSize = intermediateMTrees.size() + removed.size();
  std::vector<SimplexId> nodeSimplexId(fullSize);
  vtkNew<vtkIntArray> treeId{};
  treeId->SetName("TreeId");
  vtkNew<vtkIntArray> nodePathId{};
  nodePathId->SetName("NodePathId");
  vtkNew<vtkIntArray> pathId{};
  pathId->SetName("PathId");
  vtkNew<vtkIntArray> nodeRemoved{};
  nodeRemoved->SetName("isNodeRemoved");
  vtkNew<vtkIntArray> nodeId{};
  nodeId->SetName("NodeId");
  vtkNew<vtkIntArray> arcId{};
  arcId->SetName("ArcId");

  vtkSmartPointer<vtkUnstructuredGrid> vtkArcs
    = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for(unsigned int i = 0; i < embedding[0].size(); ++i) {
    if(not displayRemoved and i >= intermediateMTrees.size())
      continue;

    // Next point
    double point[3] = {embedding[0][i], embedding[1][i], 0};
    nodeSimplexId[i] = points->InsertNextPoint(point);

    // Add TreeID field
    int index = i;
    if(i >= intermediateMTrees.size())
      index = removed[i - intermediateMTrees.size()];
    treeId->InsertNextTuple1(index);

    // Add NodePathId field
    int thisPathId = (i < intermediateMTrees.size() ? 0 : 1);
    nodePathId->InsertNextTuple1(thisPathId);

    // Add NodeId
    nodeId->InsertNextTuple1(i);

    // Add cell
    if(i > 0 and i < intermediateMTrees.size()) {
      vtkIdType pointIds[2];
      pointIds[0] = nodeSimplexId[i - 1];
      pointIds[1] = nodeSimplexId[i];
      vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);

      // Add PathId field
      pathId->InsertNextTuple1(0);

      // Add Arc Id
      arcId->InsertNextTuple1(i);
    }
  }

  if(displayRemoved) {
    for(unsigned int i = 0; i < removed.size(); ++i) {
      int index = removed[i];
      // First cell
      if(i == 0) {
        vtkIdType pointIds[2];
        pointIds[0] = nodeSimplexId[index - 1];
        pointIds[1] = nodeSimplexId[intermediateMTrees.size() + i];
        vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
        pathId->InsertNextTuple1(1);

        if(removed[i + 1] != index + 1) {
          vtkIdType pointIds2[2];
          pointIds2[0] = nodeSimplexId[intermediateMTrees.size() + i];
          pointIds2[1] = nodeSimplexId[index + 1];
          vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds2);
          pathId->InsertNextTuple1(1);
        }
      }
      // Between cell
      if(i > 0 and i < removed.size() - 1) {
        vtkIdType pointIds[2];
        int prevIndex
          = (removed[i - 1] == index - 1 ? intermediateMTrees.size() + i - 1
                                         : index - 1);
        pointIds[0] = nodeSimplexId[prevIndex];
        pointIds[1] = nodeSimplexId[intermediateMTrees.size() + i];
        vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
        pathId->InsertNextTuple1(1);

        vtkIdType pointIds2[2];
        int nextvIndex
          = (removed[i + 1] == index + 1 ? intermediateMTrees.size() + i + 1
                                         : index + 1);
        pointIds2[0] = nodeSimplexId[intermediateMTrees.size() + i];
        pointIds2[1] = nodeSimplexId[nextvIndex];
        vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds2);
        pathId->InsertNextTuple1(1);
      }
      // End cell
      if(i == removed.size() - 1) {
        vtkIdType pointIds[2];
        pointIds[0] = nodeSimplexId[intermediateMTrees.size() + i];
        pointIds[1] = nodeSimplexId[index + 1];
        vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
        pathId->InsertNextTuple1(1);

        if(removed[i - 1] != index - 1) {
          vtkIdType pointIds2[2];
          pointIds2[0] = nodeSimplexId[index - 1];
          pointIds2[1] = nodeSimplexId[intermediateMTrees.size() + i];
          vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds2);
          pathId->InsertNextTuple1(1);
        }
      }
    }
  } else {
    for(unsigned int i = 0; i < removed.size(); ++i) {
      int before = removed[i] - 1;
      int index = i - 1;
      if(i > 0)
        while(before == removed[index]) {
          before = removed[index] - 1;
          --index;
        }

      int after = removed[i] + 1;
      index = i + 1;
      if(i < removed.size() - 1)
        while(after == removed[index]) {
          after = removed[index] + 1;
          ++index;
        }

      vtkIdType pointIds[2];
      pointIds[0] = nodeSimplexId[before];
      pointIds[1] = nodeSimplexId[after];
      vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
      pathId->InsertNextTuple1(1);
      arcId->InsertNextTuple1(intermediateMTrees.size() + i);
    }
  }

  for(unsigned int i = 0; i < fullSize; ++i)
    nodeRemoved->InsertNextTuple1(0);
  for(unsigned int i = 0; i < removed.size(); ++i) {
    nodeRemoved->SetTuple1(removed[i], 1);
    nodeRemoved->SetTuple1(intermediateMTrees.size() + i, 1);
  }

  vtkOutputNode1->SetPoints(points);
  vtkOutputNode1->GetPointData()->AddArray(treeId);
  vtkOutputNode1->GetPointData()->AddArray(nodePathId);
  vtkOutputNode1->GetPointData()->AddArray(nodeRemoved);
  vtkOutputNode1->GetPointData()->AddArray(nodeId);
  vtkArcs->SetPoints(points);
  vtkArcs->GetCellData()->AddArray(pathId);
  vtkArcs->GetCellData()->AddArray(arcId);
  vtkOutputArc1->ShallowCopy(vtkArcs);
}

template <class dataType>
void makeTemporalSubsamplingETDOutput(
  std::vector<MergeTree<dataType>> &intermediateMTrees,
  std::vector<double> &emptyTreeDistances,
  std::vector<MergeTree<dataType>> &allMT,
  std::vector<int> removed,
  vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1,
  vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1,
  double DistanceAxisStretch) {
  std::vector<std::vector<double>> embedding(2);
  for(int i = 0; i < (int)intermediateMTrees.size() + (int)removed.size();
      ++i) {
    double y = emptyTreeDistances[i] * DistanceAxisStretch;
    double x = i;
    if(i >= (int)intermediateMTrees.size())
      x = removed[i - (int)intermediateMTrees.size()];

    embedding[0].push_back(x);
    embedding[1].push_back(y);
  }

  // Create output
  makeTemporalSubsamplingOutput<dataType>(intermediateMTrees, embedding, allMT,
                                          removed, vtkOutputNode1,
                                          vtkOutputArc1);
}

#endif
