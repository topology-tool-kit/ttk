#include <ttkMeshSubdivision.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkMeshSubdivision)

  ttkMeshSubdivision::ttkMeshSubdivision() {

  // init
  IterationNumber = 1;
  UseAllCores = true;
}

ttkMeshSubdivision::~ttkMeshSubdivision() {
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkMeshSubdivision::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkMeshSubdivision::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkMeshSubdivision] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkMeshSubdivision::doIt(vtkUnstructuredGrid *input,
                             vtkUnstructuredGrid *output) {

  Timer t;

  vtkSmartPointer<vtkUnstructuredGrid> tmpGrid
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  output->DeepCopy(input);
  tmpGrid->DeepCopy(input);

  for(int i = 0; i < IterationNumber; i++) {

    vector<vector<vtkSmartPointer<vtkIdList>>> newCells(
      tmpGrid->GetNumberOfCells());
    vector<vector<vector<double>>> newPoints(tmpGrid->GetNumberOfCells());

    // for each cell, for each point, several scalar fields
    vector<vector<vector<double>>> newPointData(tmpGrid->GetNumberOfCells());
    vector<vector<vector<double>>> newCellData(tmpGrid->GetNumberOfCells());

    // make the call thread-safe
    vtkSmartPointer<vtkGenericCell> threadCell
      = vtkSmartPointer<vtkGenericCell>::New();
    tmpGrid->GetCell(0, threadCell);
    double value = 0;
    for(int j = 0; j < tmpGrid->GetPointData()->GetNumberOfArrays(); j++) {
      if(tmpGrid->GetPointData()->GetArray(j)->GetNumberOfComponents() == 1) {
        tmpGrid->GetPointData()->GetArray(j)->GetTuple(0, &value);
      }
    }
    for(int j = 0; j < tmpGrid->GetCellData()->GetNumberOfArrays(); j++) {
      if(tmpGrid->GetCellData()->GetArray(j)->GetNumberOfComponents() == 1) {
        tmpGrid->GetCellData()->GetArray(j)->GetTuple(0, &value);
      }
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId j = 0; j < (SimplexId)tmpGrid->GetNumberOfCells(); j++) {

      vtkSmartPointer<vtkGenericCell> cell
        = vtkSmartPointer<vtkGenericCell>::New();
      tmpGrid->GetCell(j, cell);

      vector<double> cellValues(tmpGrid->GetCellData()->GetNumberOfArrays());
      for(int k = 0; k < (int)cellValues.size(); k++) {
        if(tmpGrid->GetCellData()->GetArray(k)->GetNumberOfComponents() == 1) {
          tmpGrid->GetCellData()->GetArray(k)->GetTuple(j, &(cellValues[k]));
        }
      }

      newCells[j].resize(cell->GetNumberOfPoints());
      for(int k = 0; k < (int)newCells[j].size(); k++) {
        newCells[j][k] = vtkSmartPointer<vtkIdList>::New();
      }
      newCellData[j].resize(cell->GetNumberOfPoints());
      for(int k = 0; k < (int)newCellData[j].size(); k++) {
        newCellData[j][k] = cellValues;
      }

      newPoints[j].resize(cell->GetNumberOfPoints() + cell->GetNumberOfEdges()
                          + cell->GetNumberOfFaces() + 1);
      newPointData[j].resize(cell->GetNumberOfPoints()
                             + cell->GetNumberOfEdges()
                             + cell->GetNumberOfFaces() + 1);
      for(int k = 0; k < (int)newPointData[j].size(); k++) {
        newPointData[j][k].resize(tmpGrid->GetPointData()->GetNumberOfArrays());
      }

      SimplexId pointCounter = 0;
      SimplexId globalPointCounter
        = j
          * (cell->GetNumberOfPoints() + cell->GetNumberOfEdges()
             + cell->GetNumberOfFaces() + 1);

      // 0) add the vertices themselves
      vector<SimplexId> vertexMap(cell->GetNumberOfPoints());
      for(int k = 0; k < (int)cell->GetNumberOfPoints(); k++) {
        double p[3];
        tmpGrid->GetPoint(cell->GetPointId(k), p);
        newPoints[j][pointCounter].resize(3, 0);
        for(int l = 0; l < 3; l++)
          newPoints[j][pointCounter][l] = p[l];

        // copy the point data
        double val;
        for(int l = 0; l < tmpGrid->GetPointData()->GetNumberOfArrays(); l++) {
          if(tmpGrid->GetPointData()->GetArray(l)->GetNumberOfComponents()
             == 1) {
            tmpGrid->GetPointData()->GetArray(l)->GetTuple(
              cell->GetPointId(k), &val);
            newPointData[j][pointCounter][l] = val;
          }
        }

        vertexMap[k] = globalPointCounter;
        pointCounter++;
        globalPointCounter++;
      }

      vector<SimplexId> edgeMap(cell->GetNumberOfEdges());
      // 1) create the edge list and create one new vertex per edge
      for(int k = 0; k < cell->GetNumberOfEdges(); k++) {

        vtkCell *edge = cell->GetEdge(k);
        double p0[3], p1[3];
        tmpGrid->GetPoint(edge->GetPointId(0), p0);
        tmpGrid->GetPoint(edge->GetPointId(1), p1);

        edgeMap[k] = globalPointCounter;

        newPoints[j][pointCounter].resize(3);
        for(int l = 0; l < 3; l++)
          newPoints[j][pointCounter][l] = 0.5 * (p0[l] + p1[l]);

        // take care of the point data
        double value0 = 0, value1 = 0;
        for(int l = 0; l < tmpGrid->GetPointData()->GetNumberOfArrays(); l++) {
          if(tmpGrid->GetPointData()->GetArray(l)->GetNumberOfComponents()
             == 1) {
            tmpGrid->GetPointData()->GetArray(l)->GetTuple(
              edge->GetPointId(0), &value0);
            tmpGrid->GetPointData()->GetArray(l)->GetTuple(
              edge->GetPointId(1), &value1);
            newPointData[j][pointCounter][l] = 0.5 * value0 + 0.5 * value1;
          }
        }

        pointCounter++;
        globalPointCounter++;
      }

      // 2) create the face list and create one new vertex per face
      vector<SimplexId> faceMap(cell->GetNumberOfFaces());
      for(int k = 0; k < cell->GetNumberOfFaces(); k++) {

        vtkCell *face = cell->GetFace(k);
        newPoints[j][pointCounter].resize(3, 0);

        for(int l = 0; l < face->GetNumberOfPoints(); l++) {
          double p[3];
          tmpGrid->GetPoint(face->GetPointId(l), p);
          for(int m = 0; m < 3; m++) {
            newPoints[j][pointCounter][m] += p[m];
          }

          // take care of the point data
          double val = 0;
          for(int m = 0; m < tmpGrid->GetPointData()->GetNumberOfArrays();
              m++) {

            if(tmpGrid->GetPointData()->GetArray(m)->GetNumberOfComponents()
               == 1) {
              tmpGrid->GetPointData()->GetArray(m)->GetTuple(
                face->GetPointId(l), &val);
              newPointData[j][pointCounter][m]
                += (1 / ((double)face->GetNumberOfPoints())) * val;
            }
          }
        }

        faceMap[k] = globalPointCounter;

        for(int l = 0; l < 3; l++)
          newPoints[j][pointCounter][l] /= (double)face->GetNumberOfPoints();

        pointCounter++;
        globalPointCounter++;
      }

      // 3) create one new vertex per cell
      for(int k = 0; k < (int)cell->GetNumberOfPoints(); k++) {
        double p[3];
        tmpGrid->GetPoint(cell->GetPointId(k), p);
        newPoints[j][pointCounter].resize(3, 0);
        for(int l = 0; l < 3; l++)
          newPoints[j][pointCounter][l] += p[l];

        // take care of the point data
        double val = 0;
        for(int l = 0; l < tmpGrid->GetPointData()->GetNumberOfArrays(); l++) {

          if(tmpGrid->GetPointData()->GetArray(l)->GetNumberOfComponents()
             == 1) {
            tmpGrid->GetPointData()->GetArray(l)->GetTuple(
              cell->GetPointId(k), &val);
            newPointData[j][pointCounter][l]
              += (1 / ((double)cell->GetNumberOfPoints())) * val;
          }
        }
      }
      for(int k = 0; k < 3; k++)
        newPoints[j][pointCounter][k] /= (double)cell->GetNumberOfPoints();

      // 4) create the new cells:
      // for each cell, take each vertex:
      //   grab all incoming edges and faces + barycenter + vertex
      //   make a cell out of that
      for(int k = 0; k < (int)cell->GetNumberOfPoints(); k++) {

        SimplexId vertexId = cell->GetPointId(k);

        if(cell->GetCellDimension() == 2) {
          // insert the id of the vertex
          newCells[j][k]->InsertNextId(vertexMap[k]);

          // take care of the edges
          SimplexId firstEdge = -1;
          for(SimplexId l = 0; l < (SimplexId)edgeMap.size(); l++) {

            vtkCell *edge = cell->GetEdge(l);
            SimplexId vertexId0 = edge->GetPointId(0);
            SimplexId vertexId1 = edge->GetPointId(1);

            if((vertexId == vertexId0) || (vertexId == vertexId1)) {
              // add the point id to the cell
              newCells[j][k]->InsertNextId(edgeMap[l]);
              firstEdge = l;
              break;
            }
          }

          // insert the id of the barycenter of the cell
          newCells[j][k]->InsertNextId(globalPointCounter);

          // take care of the second edge
          for(SimplexId l = 0; l < (SimplexId)edgeMap.size(); l++) {

            vtkCell *edge = cell->GetEdge(l);
            SimplexId vertexId0 = edge->GetPointId(0);
            SimplexId vertexId1 = edge->GetPointId(1);

            if((l != firstEdge)
               && ((vertexId == vertexId0) || (vertexId == vertexId1))) {
              // add the point id to the cell
              newCells[j][k]->InsertNextId(edgeMap[l]);
              break;
            }
          }
        } else if(cell->GetCellDimension() == 3) {

          // insert the id of the vertex
          newCells[j][k]->InsertNextId(vertexMap[k]);

          // take care of the edges
          SimplexId firstEdge = -1;
          pair<SimplexId, SimplexId> firstEdgeVertices;
          for(SimplexId l = 0; l < (SimplexId)edgeMap.size(); l++) {

            vtkCell *edge = cell->GetEdge(l);
            SimplexId vertexId0 = edge->GetPointId(0);
            SimplexId vertexId1 = edge->GetPointId(1);

            if((vertexId == vertexId0) || (vertexId == vertexId1)) {
              // add the point id to the cell
              newCells[j][k]->InsertNextId(edgeMap[l]);
              firstEdgeVertices.first = vertexId0;
              firstEdgeVertices.second = vertexId1;
              firstEdge = l;
              break;
            }
          }

          vector<SimplexId> firstFaceVertices;
          SimplexId firstFace = -1;
          // take care of the faces
          for(SimplexId l = 0; l < (SimplexId)faceMap.size(); l++) {
            vtkCell *face = cell->GetFace(l);

            // test if this face contains the first edge
            bool vertex0In = false;
            bool vertex1In = false;
            for(int m = 0; m < (int)face->GetNumberOfPoints(); m++) {
              if(face->GetPointId(m) == firstEdgeVertices.first) {
                vertex0In = true;
              }
              if(face->GetPointId(m) == firstEdgeVertices.second) {
                vertex1In = true;
              }
            }

            if((vertex0In) && (vertex1In)) {
              newCells[j][k]->InsertNextId(faceMap[l]);

              firstFace = l;
              for(int n = 0; n < face->GetNumberOfPoints(); n++) {
                firstFaceVertices.push_back(face->GetPointId(n));
              }

              break;
            }
            if(firstFace != -1)
              break;
          }

          // take care of the second edge
          SimplexId secondEdge = -1;
          for(SimplexId l = 0; l < (SimplexId)edgeMap.size(); l++) {
            // make sure the found edge belongs to the face identified above
            // cell->GetEdge(l)->GetPointId(0)
            // cell->GetEdge(l)->GetPointId(1)
            // both should be in

            vtkCell *edge = cell->GetEdge(l);
            SimplexId vertexId0 = edge->GetPointId(0);
            SimplexId vertexId1 = edge->GetPointId(1);

            // check if this edge belongs to the face identified right before
            bool vertex0In = false;
            bool vertex1In = false;
            for(SimplexId m = 0; m < (SimplexId)firstFaceVertices.size(); m++) {
              if(firstFaceVertices[m] == vertexId0)
                vertex0In = true;
              if(firstFaceVertices[m] == vertexId1)
                vertex1In = true;
            }

            if((l != firstEdge) && (vertex0In) && (vertex1In)
               && ((vertexId == vertexId0) || (vertexId == vertexId1))) {
              // add the point id to the cell
              newCells[j][k]->InsertNextId(edgeMap[l]);

              secondEdge = l;
              break;
            }
          }

          // front face done

          // take care of the third edge
          pair<SimplexId, SimplexId> thirdEdgeVertices;
          for(SimplexId l = 0; l < (SimplexId)edgeMap.size(); l++) {

            vtkCell *edge = cell->GetEdge(l);
            SimplexId vertexId0 = edge->GetPointId(0);
            SimplexId vertexId1 = edge->GetPointId(1);

            if((l != firstEdge) && (l != secondEdge)
               && ((vertexId == vertexId0) || (vertexId == vertexId1))) {
              // add the point id to the cell
              newCells[j][k]->InsertNextId(edgeMap[l]);
              thirdEdgeVertices.first = vertexId0;
              thirdEdgeVertices.second = vertexId1;
              break;
            }
          }

          // take care of the second face
          SimplexId secondFace = -1;
          // take care of the faces
          for(SimplexId l = (SimplexId)faceMap.size() - 1; l >= 0; l--) {
            vtkCell *face = cell->GetFace(l);

            if(l != firstFace) {

              bool vertex0In = false;
              bool vertex1In = false;
              for(SimplexId m = 0; m < (SimplexId)face->GetNumberOfPoints();
                  m++) {
                if(face->GetPointId(m) == thirdEdgeVertices.first) {
                  vertex0In = true;
                }
                if(face->GetPointId(m) == thirdEdgeVertices.second) {
                  vertex1In = true;
                }
              }

              if((vertex0In) && (vertex1In)) {
                newCells[j][k]->InsertNextId(faceMap[l]);

                secondFace = l;
                break;
              }
            }
          }

          // insert the id of the barycenter of the cell
          newCells[j][k]->InsertNextId(globalPointCounter);

          // take care of the last face
          SimplexId thirdFace = -1;
          for(SimplexId l = 0; l < (SimplexId)faceMap.size(); l++) {
            vtkCell *face = cell->GetFace(l);

            if((l != firstFace) && (l != secondFace)) {
              for(int m = 0; m < (int)face->GetNumberOfPoints(); m++) {
                if(face->GetPointId(m) == vertexId) {
                  newCells[j][k]->InsertNextId(faceMap[l]);

                  thirdFace = l;
                  break;
                }
              }
              if(thirdFace != -1)
                break;
            }
          }
        }
      }
    }

    // now merge things
    vtkSmartPointer<vtkPoints> pointSet = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellArray
      = vtkSmartPointer<vtkCellArray>::New();

    vector<vtkSmartPointer<vtkDoubleArray>> pointData;
    pointData.resize(tmpGrid->GetPointData()->GetNumberOfArrays());
    for(int j = 0; j < (int)pointData.size(); j++) {
      pointData[j] = vtkSmartPointer<vtkDoubleArray>::New();
      pointData[j]->SetName(tmpGrid->GetPointData()->GetArray(j)->GetName());
    }

    vector<vtkSmartPointer<vtkDoubleArray>> cellData;
    cellData.resize(tmpGrid->GetCellData()->GetNumberOfArrays());
    for(int j = 0; j < (int)cellData.size(); j++) {
      cellData[j] = vtkSmartPointer<vtkDoubleArray>::New();
      cellData[j]->SetName(tmpGrid->GetCellData()->GetArray(j)->GetName());
    }

    // order is really important here
    for(SimplexId j = 0; j < (SimplexId)newPoints.size(); j++) {
      for(int k = 0; k < (int)newPoints[j].size(); k++) {
        pointSet->InsertNextPoint(newPoints[j][k].data());
      }
      for(int k = 0; k < (int)newPointData[j].size(); k++) {
        for(int l = 0; l < (int)newPointData[j][k].size(); l++) {
          pointData[l]->InsertNextTuple1(newPointData[j][k][l]);
        }
      }
    }
    output->SetPoints(pointSet);
    for(int j = 0; j < (int)pointData.size(); j++) {
      output->GetPointData()->AddArray(pointData[j]);
    }

    for(SimplexId j = 0; j < (SimplexId)newCells.size(); j++) {
      for(int k = 0; k < (int)newCells[j].size(); k++) {
        cellArray->InsertNextCell(newCells[j][k]);
      }

      for(int k = 0; k < (int)newCellData[j].size(); k++) {
        for(int l = 0; l < (int)newCellData[j][k].size(); l++) {
          cellData[l]->InsertNextTuple1(newCellData[j][k][l]);
        }
      }
    }
    if(tmpGrid->GetCell(0)->GetCellDimension() == 3) {
      output->SetCells(VTK_HEXAHEDRON, cellArray);
    } else if(tmpGrid->GetCell(0)->GetCellDimension() == 2) {
      output->SetCells(VTK_POLYGON, cellArray);
    }
    for(int j = 0; j < (int)cellData.size(); j++) {
      output->GetCellData()->AddArray(cellData[j]);
    }

    if(i != IterationNumber - 1) {
      tmpGrid->DeepCopy(output);
    }
  }

  {
    stringstream msg;
    msg << "[ttkMeshSubdivision] Subdivision computed in " << t.getElapsedTime()
        << " s. (" << output->GetNumberOfCells() << " cell(s))" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkMeshSubdivision::RequestData(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {

  Memory m;

  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkUnstructuredGrid *input = vtkUnstructuredGrid::GetData(inputVector[0]);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkMeshSubdivision] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
