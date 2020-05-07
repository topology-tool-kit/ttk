#include <ttkPointMerger.h>

#include <Geometry.h>

#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkGenericCell.h>


vtkStandardNewMacro(ttkPointMerger);

ttkPointMerger::ttkPointMerger() {
  this->setDebugMsgPrefix("PointMerger");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkPointMerger::~ttkPointMerger() {
}

int ttkPointMerger::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkPointMerger::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkPointMerger::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  auto input = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  ttk::Timer t;

  auto triangulation = ttkAlgorithm::GetTriangulation(input);

  if(!triangulation)
    return -1;

  if(BoundaryOnly)
    triangulation->preconditionBoundaryVertices();

  ttk::SimplexId vertexNumber = input->GetNumberOfPoints();
  std::vector<ttk::SimplexId> candidateVertices;

  if(BoundaryOnly) {
    for(ttk::SimplexId i = 0; i < vertexNumber; i++) {
      if(triangulation->isVertexOnBoundary(i)) {
        candidateVertices.push_back(i);
      }
    }
  } else {
    candidateVertices.resize(vertexNumber);
    for(ttk::SimplexId i = 0; i < vertexNumber; i++)
      candidateVertices[i] = i;
  }

  std::vector<std::vector<ttk::SimplexId>> closePoints(vertexNumber);

  this->printMsg(
    "Computing pointwise distances (" + std::to_string(candidateVertices.size()) +  " candidates)",
    0, 0, this->threadNumber_,
    ttk::debug::LineMode::REPLACE
  );

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(ttk::SimplexId i = 0; i < (ttk::SimplexId)candidateVertices.size(); i++) {
    double distance = -1;
    std::vector<double> p0(3);

    ttk::SimplexId vertexId0 = candidateVertices[i];

    input->GetPoint(vertexId0, p0.data());
    for(ttk::SimplexId j = 0; j < (ttk::SimplexId)candidateVertices.size(); j++) {
      if(i != j) {
        std::vector<double> p1(3);
        ttk::SimplexId vertexId1 = candidateVertices[j];
        input->GetPoint(vertexId1, p1.data());
        distance = ttk::Geometry::distance(p0.data(), p1.data());
        if(distance < DistanceThreshold) {
          closePoints[vertexId0].push_back(vertexId1);
        }
      }
    }
  }

  this->printMsg(
    "Computing pointwise distances (" + std::to_string(candidateVertices.size()) +  " candidates)",
    1, t.getElapsedTime(), this->threadNumber_
  );

  t.reStart();
  this->printMsg(
    "Merging Points",
    0, 0, this->threadNumber_,
    ttk::debug::LineMode::REPLACE
  );

  std::vector<double> minMergeDistance(vertexNumber, -1);
  std::vector<double> maxMergeDistance(vertexNumber, -1);
  std::vector<ttk::SimplexId> mergeCount(vertexNumber, 0);
  std::vector<ttk::SimplexId> mergeMap(vertexNumber, -1);

  for(ttk::SimplexId i = 0; i < vertexNumber; i++) {
    if(closePoints[i].size()) {

      ttk::SimplexId targetVertexId = i;
      do {
        ttk::SimplexId nextTargetVertexId = targetVertexId;
        for(ttk::SimplexId j = 0; j < (ttk::SimplexId)closePoints[targetVertexId].size();
            j++) {
          if(closePoints[targetVertexId][j] < nextTargetVertexId)
            nextTargetVertexId = closePoints[targetVertexId][j];
        }
        targetVertexId = nextTargetVertexId;
        mergeMap[i] = targetVertexId;
      } while(mergeMap[targetVertexId] != targetVertexId);
      mergeCount[targetVertexId]++;

      if(i != targetVertexId) {
        std::vector<double> p0(3), p1(3);
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

  ttk::SimplexId vertexIdGen = 0;
  std::vector<ttk::SimplexId> old2new(vertexNumber, 0);
  std::vector<ttk::SimplexId> new2old;
  for(ttk::SimplexId i = 0; i < vertexNumber; i++) {
    if(mergeMap[i] == i) {
      old2new[i] = vertexIdGen;
      vertexIdGen++;
      new2old.push_back(i);
    } else {
      old2new[i] = old2new[mergeMap[i]];
    }
  }

  // now create the output
  vtkSmartPointer<vtkPoints> pointSet = vtkSmartPointer<vtkPoints>::New();
  pointSet->SetNumberOfPoints(vertexIdGen);

  vtkSmartPointer<vtkIdTypeArray> mergeCountArray
    = vtkSmartPointer<vtkIdTypeArray>::New();
  mergeCountArray->SetNumberOfTuples(vertexIdGen);
  mergeCountArray->SetName("VertexMergeCount");

  vtkSmartPointer<vtkDoubleArray> minDistanceArray
    = vtkSmartPointer<vtkDoubleArray>::New();
  minDistanceArray->SetNumberOfTuples(vertexIdGen);
  minDistanceArray->SetName("MinMergeDistance");

  vtkSmartPointer<vtkDoubleArray> maxDistanceArray
    = vtkSmartPointer<vtkDoubleArray>::New();
  maxDistanceArray->SetNumberOfTuples(vertexIdGen);
  maxDistanceArray->SetName("MaxMergeDistance");

  std::vector<vtkSmartPointer<vtkDoubleArray>> pointData;
  pointData.resize(input->GetPointData()->GetNumberOfArrays());
  for(int i = 0; i < (int)pointData.size(); i++) {
    pointData[i] = vtkSmartPointer<vtkDoubleArray>::New();
    pointData[i]->SetName(input->GetPointData()->GetArray(i)->GetName());
    pointData[i]->SetNumberOfComponents(
      input->GetPointData()->GetArray(i)->GetNumberOfComponents());
    pointData[i]->SetNumberOfTuples(vertexIdGen);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(ttk::SimplexId i = 0; i < vertexIdGen; i++) {
    std::vector<double> p(3);
    input->GetPoint(new2old[i], p.data());
    pointSet->SetPoint(i, p.data());

    mergeCountArray->SetTuple1(i, mergeCount[new2old[i]]);
    minDistanceArray->SetTuple1(i, minMergeDistance[new2old[i]]);
    maxDistanceArray->SetTuple1(i, maxMergeDistance[new2old[i]]);

    for(int j = 0; j < (int)pointData.size(); j++) {
      std::vector<double> data(pointData[j]->GetNumberOfComponents());
      input->GetPointData()->GetArray(j)->GetTuple(new2old[i], data.data());
      pointData[j]->SetTuple(i, data.data());
    }
  }

  output->SetPoints(pointSet);
  for(int i = 0; i < (int)pointData.size(); i++) {
    output->GetPointData()->AddArray(pointData[i]);
  }
  output->GetPointData()->AddArray(mergeCountArray);
  output->GetPointData()->AddArray(minDistanceArray);
  output->GetPointData()->AddArray(maxDistanceArray);

  // now do the cells
  output->Allocate();
  std::vector<vtkSmartPointer<vtkDoubleArray>> cellData;
  cellData.resize(input->GetCellData()->GetNumberOfArrays());
  for(int i = 0; i < (int)cellData.size(); i++) {
    cellData[i] = vtkSmartPointer<vtkDoubleArray>::New();
    cellData[i]->SetName(input->GetCellData()->GetArray(i)->GetName());
    cellData[i]->SetNumberOfComponents(
      input->GetCellData()->GetArray(i)->GetNumberOfComponents());
  }
  for(ttk::SimplexId i = 0; i < input->GetNumberOfCells(); i++) {

    vtkSmartPointer<vtkGenericCell> c = vtkSmartPointer<vtkGenericCell>::New();
    input->GetCell(i, c);

    std::vector<ttk::SimplexId> newVertexIds;
    for(int j = 0; j < c->GetNumberOfPoints(); j++) {
      ttk::SimplexId vertexId = old2new[c->GetPointId(j)];
      bool isIn = false;
      for(ttk::SimplexId k = 0; k < (ttk::SimplexId)newVertexIds.size(); k++) {
        if(newVertexIds[k] == vertexId) {
          isIn = true;
          break;
        }
      }
      if(!isIn) {
        newVertexIds.push_back(vertexId);
      }
    }

    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    idList->SetNumberOfIds(newVertexIds.size());
    for(ttk::SimplexId j = 0; j < (ttk::SimplexId)newVertexIds.size(); j++) {
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
        this->printWrn("Ill-defined cell type for cell #" + std::to_string(i) + " (" + std::to_string(newVertexIds.size()) + " vertices)");
      }
    }

    if(cellId != -1) {
      // insert the cell data
      for(int j = 0; j < (int)cellData.size(); j++) {
        std::vector<double> data(cellData[j]->GetNumberOfComponents());
        input->GetCellData()->GetArray(j)->GetTuple(i, data.data());
        cellData[j]->InsertNextTuple(data.data());
      }
    }
  }
  for(int i = 0; i < (int)cellData.size(); i++)
    output->GetCellData()->AddArray(cellData[i]);

  this->printMsg(
    "Merging Points",
    1, t.getElapsedTime(), this->threadNumber_
  );

  return 1;
}
