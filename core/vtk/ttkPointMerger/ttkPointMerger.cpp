#include <ttkPointMerger.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPointMerger)

  int ttkPointMerger::doIt(vector<vtkDataSet *> &inputs,
                           vector<vtkDataSet *> &outputs) {

  Timer t;
  Memory m;

  string dataType = inputs[0]->GetClassName();
  if(dataType != "vtkUnstructuredGrid") {
    stringstream msg;
    msg << "[ttkPointMerger] The input is not a 'vtkUnstructuredGrid'!" << endl;
    msg << "[ttkPointMerger] Doing nothing..." << endl;
    dMsg(cerr, msg.str(), Debug::infoMsg);

    outputs[0]->ShallowCopy(inputs[0]);

    return 0;
  }

  vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

  if(!triangulation)
    return -1;

  triangulation->setWrapper(this);
  if(BoundaryOnly)
    triangulation->preprocessBoundaryVertices();

  SimplexId vertexNumber = input->GetNumberOfPoints();
  vector<SimplexId> candidateVertices;

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

  vector<vector<SimplexId>> closePoints(vertexNumber);

  {
    stringstream msg;
    msg << "[ttkPointMerger] Computing pointwise distances ("
        << candidateVertices.size() << " candidates)..." << endl;
    dMsg(cout, msg.str(), Debug::timeMsg);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)candidateVertices.size(); i++) {
    double distance = -1;
    vector<double> p0(3);

    SimplexId vertexId0 = candidateVertices[i];

    input->GetPoint(vertexId0, p0.data());
    for(SimplexId j = 0; j < (SimplexId)candidateVertices.size(); j++) {
      if(i != j) {
        vector<double> p1(3);
        SimplexId vertexId1 = candidateVertices[j];
        input->GetPoint(vertexId1, p1.data());
        distance = Geometry::distance(p0.data(), p1.data());
        if(distance < DistanceThreshold) {
          closePoints[vertexId0].push_back(vertexId1);
        }
      }
    }
  }

  {
    stringstream msg;
    msg << "[ttkPointMerger] Merging points..." << endl;
    dMsg(cout, msg.str(), Debug::timeMsg);
  }

  vector<double> minMergeDistance(vertexNumber, -1);
  vector<double> maxMergeDistance(vertexNumber, -1);
  vector<SimplexId> mergeCount(vertexNumber, 0);
  vector<SimplexId> mergeMap(vertexNumber, -1);

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
        vector<double> p0(3), p1(3);
        input->GetPoint(i, p0.data());
        input->GetPoint(targetVertexId, p1.data());
        double distance = Geometry::distance(p0.data(), p1.data());
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
  vector<SimplexId> old2new(vertexNumber, 0);
  vector<SimplexId> new2old;
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

  vector<vtkSmartPointer<vtkDoubleArray>> pointData;
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
  for(SimplexId i = 0; i < vertexIdGen; i++) {
    vector<double> p(3);
    input->GetPoint(new2old[i], p.data());
    pointSet->SetPoint(i, p.data());

    mergeCountArray->SetTuple1(i, mergeCount[new2old[i]]);
    minDistanceArray->SetTuple1(i, minMergeDistance[new2old[i]]);
    maxDistanceArray->SetTuple1(i, maxMergeDistance[new2old[i]]);

    for(int j = 0; j < (int)pointData.size(); j++) {
      vector<double> data(pointData[j]->GetNumberOfComponents());
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
  vector<vtkSmartPointer<vtkDoubleArray>> cellData;
  cellData.resize(input->GetCellData()->GetNumberOfArrays());
  for(int i = 0; i < (int)cellData.size(); i++) {
    cellData[i] = vtkSmartPointer<vtkDoubleArray>::New();
    cellData[i]->SetName(input->GetCellData()->GetArray(i)->GetName());
    cellData[i]->SetNumberOfComponents(
      input->GetCellData()->GetArray(i)->GetNumberOfComponents());
  }
  for(SimplexId i = 0; i < input->GetNumberOfCells(); i++) {

    vtkSmartPointer<vtkGenericCell> c = vtkSmartPointer<vtkGenericCell>::New();
    input->GetCell(i, c);

    vector<SimplexId> newVertexIds;
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

    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    idList->SetNumberOfIds(newVertexIds.size());
    for(SimplexId j = 0; j < (SimplexId)newVertexIds.size(); j++) {
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
        stringstream msg;
        msg << "[ttkPointMerger] Ill-defined cell type for cell #" << i << " ("
            << newVertexIds.size() << " vertices)" << endl;
        dMsg(cerr, msg.str(), Debug::infoMsg);
      }
    }

    if(cellId != -1) {
      // insert the cell data
      for(int j = 0; j < (int)cellData.size(); j++) {
        vector<double> data(cellData[j]->GetNumberOfComponents());
        input->GetCellData()->GetArray(j)->GetTuple(i, data.data());
        cellData[j]->InsertNextTuple(data.data());
      }
    }
  }
  for(int i = 0; i < (int)cellData.size(); i++)
    output->GetCellData()->AddArray(cellData[i]);

  {
    stringstream msg;
    msg << "[ttkPointMerger] Merge performed in " << t.getElapsedTime() << " s."
        << endl;
    msg << "[ttkPointMerger] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
