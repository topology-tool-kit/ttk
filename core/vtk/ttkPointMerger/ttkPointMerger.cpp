#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include <Geometry.h>
#include <Triangulation.h>
#include <ttkPointMerger.h>

#include <array>

vtkStandardNewMacro(ttkPointMerger);

ttkPointMerger::ttkPointMerger() {
  this->setDebugMsgPrefix("PointMerger");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkPointMerger::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkPointMerger::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkPointMerger::RequestData(vtkInformation *ttkNotUsed(request),
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {

  using ttk::SimplexId;
  ttk::Timer t;

  auto input = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  auto triangulation = ttkAlgorithm::GetTriangulation(input);

  if(triangulation == nullptr) {
    return -1;
  }

  if(BoundaryOnly) {
    triangulation->preconditionBoundaryVertices();
  }

  SimplexId vertexNumber = input->GetNumberOfPoints();
  std::vector<SimplexId> candidateVertices;

  if(BoundaryOnly) {
    for(SimplexId i = 0; i < vertexNumber; i++) {
      if(triangulation->isVertexOnBoundary(i)) {
        candidateVertices.push_back(i);
      }
    }
  } else {
    candidateVertices.resize(vertexNumber);
    for(SimplexId i = 0; i < vertexNumber; i++)
      candidateVertices[i] = i;
  }

  std::vector<std::vector<SimplexId>> closePoints(vertexNumber);

  this->printMsg("Computing pointwise distances with "
                 + std::to_string(candidateVertices.size()) + " candidates");

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)candidateVertices.size(); i++) {
    double distance = -1;
    std::array<double, 3> p0{};

    SimplexId vertexId0 = candidateVertices[i];

    input->GetPoint(vertexId0, p0.data());
    for(SimplexId j = 0; j < (SimplexId)candidateVertices.size(); j++) {
      if(i != j) {
        std::array<double, 3> p1{};
        SimplexId vertexId1 = candidateVertices[j];
        input->GetPoint(vertexId1, p1.data());
        distance = ttk::Geometry::distance(p0.data(), p1.data());
        if(distance < DistanceThreshold) {
          closePoints[vertexId0].push_back(vertexId1);
        }
      }
    }
  }

  this->printMsg("Merging points...");

  std::vector<double> minMergeDistance(vertexNumber, -1);
  std::vector<double> maxMergeDistance(vertexNumber, -1);
  std::vector<SimplexId> mergeCount(vertexNumber, 0);
  std::vector<SimplexId> mergeMap(vertexNumber, -1);

  for(SimplexId i = 0; i < vertexNumber; i++) {
    if(closePoints[i].size()) {

      SimplexId targetVertexId = i;
      do {
        SimplexId nextTargetVertexId = targetVertexId;
        for(SimplexId j = 0; j < (SimplexId)closePoints[targetVertexId].size();
            j++) {
          if(closePoints[targetVertexId][j] < nextTargetVertexId)
            nextTargetVertexId = closePoints[targetVertexId][j];
        }
        targetVertexId = nextTargetVertexId;
        mergeMap[i] = targetVertexId;
      } while(mergeMap[targetVertexId] != targetVertexId);
      mergeCount[targetVertexId]++;

      if(i != targetVertexId) {
        std::array<double, 3> p0{}, p1{};
        input->GetPoint(i, p0.data());
        input->GetPoint(targetVertexId, p1.data());
        double distance = ttk::Geometry::distance(p0.data(), p1.data());
        if((minMergeDistance[targetVertexId] == -1)
           || (distance < minMergeDistance[targetVertexId])) {
          minMergeDistance[targetVertexId] = distance;
        }
        if((maxMergeDistance[targetVertexId] == -1)
           || (distance > maxMergeDistance[targetVertexId])) {
          maxMergeDistance[targetVertexId] = distance;
        }
      }
    } else {
      mergeMap[i] = i;
    }
  }

  SimplexId vertexIdGen = 0;
  std::vector<SimplexId> old2new(vertexNumber, 0);
  std::vector<SimplexId> new2old;
  for(SimplexId i = 0; i < vertexNumber; i++) {
    if(mergeMap[i] == i) {
      old2new[i] = vertexIdGen;
      vertexIdGen++;
      new2old.push_back(i);
    } else {
      old2new[i] = old2new[mergeMap[i]];
    }
  }

  // now create the output
  vtkNew<vtkPoints> pointSet{};
  pointSet->SetNumberOfPoints(vertexIdGen);

  vtkNew<vtkIdTypeArray> mergeCountArray{};
  mergeCountArray->SetNumberOfTuples(vertexIdGen);
  mergeCountArray->SetName("VertexMergeCount");

  vtkNew<vtkDoubleArray> minDistanceArray{};
  minDistanceArray->SetNumberOfTuples(vertexIdGen);
  minDistanceArray->SetName("MinMergeDistance");

  vtkNew<vtkDoubleArray> maxDistanceArray{};
  maxDistanceArray->SetNumberOfTuples(vertexIdGen);
  maxDistanceArray->SetName("MaxMergeDistance");

  std::vector<vtkNew<vtkDoubleArray>> pointData{};
  pointData.resize(input->GetPointData()->GetNumberOfArrays());
  for(size_t i = 0; i < pointData.size(); i++) {
    pointData[i]->SetName(input->GetPointData()->GetArray(i)->GetName());
    pointData[i]->SetNumberOfComponents(
      input->GetPointData()->GetArray(i)->GetNumberOfComponents());
    pointData[i]->SetNumberOfTuples(vertexIdGen);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < vertexIdGen; i++) {
    std::array<double, 3> p{};
    input->GetPoint(new2old[i], p.data());
    pointSet->SetPoint(i, p.data());

    mergeCountArray->SetTuple1(i, mergeCount[new2old[i]]);
    minDistanceArray->SetTuple1(i, minMergeDistance[new2old[i]]);
    maxDistanceArray->SetTuple1(i, maxMergeDistance[new2old[i]]);

    for(size_t j = 0; j < pointData.size(); j++) {
      std::vector<double> data(pointData[j]->GetNumberOfComponents());
      input->GetPointData()->GetArray(j)->GetTuple(new2old[i], data.data());
      pointData[j]->SetTuple(i, data.data());
    }
  }

  output->SetPoints(pointSet);
  for(size_t i = 0; i < pointData.size(); i++) {
    output->GetPointData()->AddArray(pointData[i]);
  }
  output->GetPointData()->AddArray(mergeCountArray);
  output->GetPointData()->AddArray(minDistanceArray);
  output->GetPointData()->AddArray(maxDistanceArray);

  // now do the cells
  output->Allocate();
  std::vector<vtkNew<vtkDoubleArray>> cellData{};
  cellData.resize(input->GetCellData()->GetNumberOfArrays());
  for(size_t i = 0; i < cellData.size(); i++) {
    cellData[i]->SetName(input->GetCellData()->GetArray(i)->GetName());
    cellData[i]->SetNumberOfComponents(
      input->GetCellData()->GetArray(i)->GetNumberOfComponents());
  }

  for(SimplexId i = 0; i < input->GetNumberOfCells(); i++) {
    vtkNew<vtkGenericCell> c{};
    input->GetCell(i, c);

    std::vector<SimplexId> newVertexIds;
    for(int j = 0; j < c->GetNumberOfPoints(); j++) {
      SimplexId vertexId = old2new[c->GetPointId(j)];
      bool isIn = false;
      for(SimplexId k = 0; k < (SimplexId)newVertexIds.size(); k++) {
        if(newVertexIds[k] == vertexId) {
          isIn = true;
          break;
        }
      }
      if(!isIn) {
        newVertexIds.push_back(vertexId);
      }
    }

    vtkNew<vtkIdList> idList{};
    idList->SetNumberOfIds(newVertexIds.size());
    for(size_t j = 0; j < newVertexIds.size(); j++) {
      idList->SetId(j, newVertexIds[j]);
    }

    vtkIdType cellId = -1;

    if((c->GetCellDimension() == 1) && (newVertexIds.size() == 2)) {
      cellId = output->InsertNextCell(VTK_LINE, idList);
    }
    if(c->GetCellDimension() == 2) {
      if(newVertexIds.size() > 2) {
        cellId = output->InsertNextCell(VTK_POLYGON, idList);
      }
    }
    if(c->GetCellDimension() == 3) {
      if(newVertexIds.size() == 4) {
        cellId = output->InsertNextCell(VTK_TETRA, idList);
      }
      if(newVertexIds.size() == 8) {
        cellId = output->InsertNextCell(VTK_HEXAHEDRON, idList);
      } else {
        this->printWrn("Ill-defined cell type for cell #" + std::to_string(i)
                       + " (" + std::to_string(newVertexIds.size())
                       + " vertices)");
      }
    }

    if(cellId != -1) {
      // insert the cell data
      for(size_t j = 0; j < cellData.size(); j++) {
        std::vector<double> data(cellData[j]->GetNumberOfComponents());
        input->GetCellData()->GetArray(j)->GetTuple(i, data.data());
        cellData[j]->InsertNextTuple(data.data());
      }
    }
  }
  for(size_t i = 0; i < cellData.size(); i++)
    output->GetCellData()->AddArray(cellData[i]);

  this->printMsg(
    "Merge performed", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 1;
}
