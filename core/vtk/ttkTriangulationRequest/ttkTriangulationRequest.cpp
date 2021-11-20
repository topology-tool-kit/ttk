#include <ttkTriangulationRequest.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkSignedCharArray.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkTriangulationRequest);

ttkTriangulationRequest::ttkTriangulationRequest() {
  this->setDebugMsgPrefix("TriangulationRequest");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkTriangulationRequest::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTriangulationRequest::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkTriangulationRequest::RequestData(vtkInformation *ttkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  ttk::Timer timer;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector);

  const auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  const int dimensionality = triangulation->getDimensionality();

  vtkNew<vtkPoints> points{};
  vtkNew<vtkUnstructuredGrid> cells{};
  vtkNew<ttkSimplexIdTypeArray> cellIds{};
  cellIds->SetNumberOfComponents(1);
  cellIds->SetName("CellId");
  vtkNew<vtkSignedCharArray> cellDims{};
  cellDims->SetNumberOfComponents(1);
  cellDims->SetName("CellDimension");

  using ttk::SimplexId;

  std::vector<SimplexId> vertices;
  const SimplexId numberOfVertices = triangulation->getNumberOfVertices();
  std::vector<SimplexId> isVisited(numberOfVertices, -1);

  this->printMsg(ttk::debug::Separator::L1);
  this->printMsg({
    {"#Threads", std::to_string(this->threadNumber_)},
    {"#Vertices", std::to_string(numberOfVertices)},
  });

  auto addVertex = [&](const SimplexId vertexId) {
    if(vertexId == -1)
      return vtkIdType(-1);

    std::array<float, 3> p{};
    triangulation->getVertexPoint(vertexId, p[0], p[1], p[2]);
    vertices.push_back(vertexId);
    return points->InsertNextPoint(p.data());
  };

  auto addEdge = [&](const SimplexId edgeId) {
    std::array<vtkIdType, 2> pointIds{};
    std::array<SimplexId, 2> vertexIds{};

    for(int i = 0; i < 2; ++i) {
      triangulation->getEdgeVertex(edgeId, i, vertexIds[i]);

      const SimplexId vertexId = vertexIds[i];
      if(vertexId == -1)
        return -1;

      if(isVisited[vertexId] == -1) {
        pointIds[i] = addVertex(vertexId);
        isVisited[vertexId] = pointIds[i];
      } else
        pointIds[i] = isVisited[vertexId];
    }

    cells->InsertNextCell(VTK_LINE, 2, pointIds.data());
    cellIds->InsertNextTuple1(edgeId);
    cellDims->InsertNextTuple1(1);

    return 0;
  };

  auto addTriangle = [&](const SimplexId triangleId) {
    std::array<vtkIdType, 3> pointIds{};
    std::array<SimplexId, 3> vertexIds{};

    for(int i = 0; i < 3; ++i) {
      if(dimensionality == 3)
        triangulation->getTriangleVertex(triangleId, i, vertexIds[i]);
      else
        triangulation->getCellVertex(triangleId, i, vertexIds[i]);

      const SimplexId vertexId = vertexIds[i];
      if(vertexId == -1)
        return -1;

      if(isVisited[vertexId] == -1) {
        pointIds[i] = addVertex(vertexId);
        isVisited[vertexId] = pointIds[i];
      } else
        pointIds[i] = isVisited[vertexId];
    }

    cells->InsertNextCell(VTK_TRIANGLE, 3, pointIds.data());
    cellIds->InsertNextTuple1(triangleId);
    cellDims->InsertNextTuple1(2);

    return 0;
  };

  auto addTetra = [&](const SimplexId tetraId) {
    std::array<vtkIdType, 4> pointIds{};
    std::array<SimplexId, 4> vertexIds{};

    for(int i = 0; i < 4; ++i) {
      triangulation->getCellVertex(tetraId, i, vertexIds[i]);

      const SimplexId vertexId = vertexIds[i];
      if(vertexId == -1)
        return -1;

      if(isVisited[vertexId] == -1) {
        pointIds[i] = addVertex(vertexId);
        isVisited[vertexId] = pointIds[i];
      } else
        pointIds[i] = isVisited[vertexId];
    }

    cells->InsertNextCell(VTK_TETRA, 4, pointIds.data());
    cellIds->InsertNextTuple1(tetraId);
    cellDims->InsertNextTuple1(3);

    return 0;
  };

  auto addStar = [&](const SimplexId starId) {
    if(dimensionality == 2)
      addTriangle(starId);
    else if(dimensionality == 3)
      addTetra(starId);
  };

  // parse SimplexIdentifier into a vector of identifiers
  std::vector<SimplexId> ids{};
  std::istringstream iss(this->SimplexIdentifier);
  for(SimplexId i; iss >> i;) {
    ids.emplace_back(i);
    if(iss.peek() == ',') {
      iss.ignore();
    }
  }

  // do minimum preprocess and put watchdog on SimplexIdentifier
  for(const auto si : ids) {
    switch(this->SimplexType) {
      case SIMPLEX::VERTEX:
        if(si < 0 or si >= numberOfVertices) {
          this->printErr("Vertex ID beyond the range.");
          return 0;
        }
        break;

      case SIMPLEX::EDGE:
        triangulation->preconditionEdges();
        if(si < 0 or si >= triangulation->getNumberOfEdges()) {
          this->printErr("Edge ID beyond the range.");
          return 0;
        }
        break;

      case SIMPLEX::TRIANGLE:
        if(dimensionality == 2) {
          if(si < 0 or si >= triangulation->getNumberOfCells()) {
            this->printErr("Triangle ID beyond the range.");
            return 0;
          }
        } else if(dimensionality == 3) {
          triangulation->preconditionTriangles();
          if(si < 0 or si >= triangulation->getNumberOfTriangles()) {
            this->printErr("Triangle ID beyond the range.");
            return 0;
          }
        }
        break;

      case SIMPLEX::TETRA:
        if(dimensionality == 3)
          if(si < 0 or si >= triangulation->getNumberOfCells()) {
            this->printErr("Tetrahedron ID beyond the range.");
            return 0;
          }
        break;
    }

    switch(this->RequestType) {
      case REQUEST::COMPUTE_SIMPLEX:
        switch(this->SimplexType) {
          case SIMPLEX::VERTEX: {
            const auto vid = addVertex(si);
            cells->InsertNextCell(VTK_VERTEX, 1, &vid);
            cellIds->InsertNextTuple1(vid);
            cellDims->InsertNextTuple1(0);
          } break;

          case SIMPLEX::EDGE:
            addEdge(si);
            break;

          case SIMPLEX::TRIANGLE:
            addTriangle(si);
            break;

          case SIMPLEX::TETRA:
            if(dimensionality == 3)
              addTetra(si);
            break;
        }
        break;

      case REQUEST::COMPUTE_FACET:
        switch(this->SimplexType) {
          case SIMPLEX::VERTEX:
            break;

          case SIMPLEX::EDGE:
            for(int i = 0; i < 2; ++i) {
              SimplexId vertexId;
              triangulation->getEdgeVertex(si, i, vertexId);
              addVertex(vertexId);
            }
            break;

          case SIMPLEX::TRIANGLE:
            if(dimensionality == 2) {
              triangulation->preconditionCellEdges();
              for(int i = 0; i < 3; ++i) {
                SimplexId edgeId;
                triangulation->getCellEdge(si, i, edgeId);
                addEdge(edgeId);
              }
            } else if(dimensionality == 3) {
              triangulation->preconditionTriangleEdges();
              for(int i = 0; i < 3; ++i) {
                SimplexId edgeId;
                triangulation->getTriangleEdge(si, i, edgeId);
                addEdge(edgeId);
              }
            }
            break;

          case SIMPLEX::TETRA:
            if(dimensionality == 3) {
              triangulation->preconditionCellTriangles();
              for(int i = 0; i < 4; ++i) {
                SimplexId triangleId;
                triangulation->getCellTriangle(si, i, triangleId);
                addTriangle(triangleId);
              }
            }
            break;
        }
        break;

      case REQUEST::COMPUTE_COFACET:
        switch(this->SimplexType) {
          case SIMPLEX::VERTEX:
            triangulation->preconditionVertexNeighbors();
            triangulation->preconditionVertexEdges();
            {
              const SimplexId edgeNumber
                = triangulation->getVertexEdgeNumber(si);
              for(SimplexId i = 0; i < edgeNumber; ++i) {
                SimplexId edgeId;
                triangulation->getVertexEdge(si, i, edgeId);
                addEdge(edgeId);
              }
            }
            break;

          case SIMPLEX::EDGE:
            if(dimensionality == 2) {
              triangulation->preconditionEdgeStars();
              const SimplexId starNumber = triangulation->getEdgeStarNumber(si);
              for(SimplexId i = 0; i < starNumber; ++i) {
                SimplexId starId;
                triangulation->getEdgeStar(si, i, starId);
                addStar(starId);
              }
            } else if(dimensionality == 3) {
              triangulation->preconditionEdgeTriangles();
              const SimplexId triangleNumber
                = triangulation->getEdgeTriangleNumber(si);
              for(SimplexId i = 0; i < triangleNumber; ++i) {
                SimplexId triangleId;
                triangulation->getEdgeTriangle(si, i, triangleId);
                addTriangle(triangleId);
              }
            }
            break;

          case SIMPLEX::TRIANGLE:
            if(dimensionality == 3) {
              triangulation->preconditionTriangleStars();
              const SimplexId starNumber
                = triangulation->getTriangleStarNumber(si);
              for(SimplexId i = 0; i < starNumber; ++i) {
                SimplexId starId;
                triangulation->getTriangleStar(si, i, starId);
                addStar(starId);
              }
            }
            break;

          case SIMPLEX::TETRA:
            break;
        }
        break;

      case REQUEST::COMPUTE_STAR:
        switch(this->SimplexType) {
          case SIMPLEX::VERTEX:
            triangulation->preconditionVertexStars();
            {
              const SimplexId starNumber
                = triangulation->getVertexStarNumber(si);
              for(SimplexId i = 0; i < starNumber; ++i) {
                SimplexId starId;
                triangulation->getVertexStar(si, i, starId);
                addStar(starId);
              }
            }
            break;

          case SIMPLEX::EDGE:
            triangulation->preconditionEdgeStars();
            {
              const SimplexId starNumber = triangulation->getEdgeStarNumber(si);
              for(SimplexId i = 0; i < starNumber; ++i) {
                SimplexId starId;
                triangulation->getEdgeStar(si, i, starId);
                addStar(starId);
              }
            }
            break;

          case SIMPLEX::TRIANGLE:
            if(dimensionality == 3) {
              triangulation->preconditionTriangleStars();
              const SimplexId starNumber
                = triangulation->getTriangleStarNumber(si);
              for(SimplexId i = 0; i < starNumber; ++i) {
                SimplexId starId;
                triangulation->getTriangleStar(si, i, starId);
                addStar(starId);
              }
            }
            break;

          case SIMPLEX::TETRA:
            break;
        }
        break;

      case REQUEST::COMPUTE_LINK:
        switch(this->SimplexType) {
          case SIMPLEX::VERTEX:
            triangulation->preconditionVertexLinks();
            if(dimensionality == 2) {
              triangulation->preconditionEdges();
              const SimplexId linkNumber
                = triangulation->getVertexLinkNumber(si);
              for(SimplexId i = 0; i < linkNumber; ++i) {
                SimplexId linkId;
                triangulation->getVertexLink(si, i, linkId);
                addEdge(linkId);
              }
            } else if(dimensionality == 3) {
              triangulation->preconditionVertexTriangles();
              const SimplexId linkNumber
                = triangulation->getVertexLinkNumber(si);
              for(SimplexId i = 0; i < linkNumber; ++i) {
                SimplexId linkId;
                triangulation->getVertexLink(si, i, linkId);
                addTriangle(linkId);
              }
            }
            break;

          case SIMPLEX::EDGE:
            triangulation->preconditionEdgeLinks();
            if(dimensionality == 2) {
              const SimplexId linkNumber = triangulation->getEdgeLinkNumber(si);
              for(SimplexId i = 0; i < linkNumber; ++i) {
                SimplexId linkId;
                triangulation->getEdgeLink(si, i, linkId);
                addVertex(linkId);
              }
            } else if(dimensionality == 3) {
              const SimplexId linkNumber = triangulation->getEdgeLinkNumber(si);
              for(SimplexId i = 0; i < linkNumber; ++i) {
                SimplexId linkId;
                triangulation->getEdgeLink(si, i, linkId);
                addEdge(linkId);
              }
            }
            break;

          case SIMPLEX::TRIANGLE:
            if(dimensionality == 3) {
              triangulation->preconditionTriangleLinks();
              const SimplexId linkNumber
                = triangulation->getTriangleLinkNumber(si);
              for(SimplexId i = 0; i < linkNumber; ++i) {
                SimplexId linkId;
                triangulation->getTriangleLink(si, i, linkId);
                addVertex(linkId);
              }
            }
            break;

          case SIMPLEX::TETRA:
            break;
        }
        break;
    }
  }

  cells->SetPoints(points);
  cells->GetCellData()->AddArray(cellIds);
  cells->GetCellData()->AddArray(cellDims);

  output->ShallowCopy(cells);

  if(KeepAllDataArrays) {
    vtkPointData *inputPointData = input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!inputPointData) {
      this->printErr("Input has no point data.");
      return 0;
    }
#endif
    const int numberOfInputArrays = inputPointData->GetNumberOfArrays();

    vtkPointData *outputPointData = output->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!outputPointData) {
      this->printErr("Output has no point data.");
      return 0;
    }
#endif

    for(int i = 0; i < numberOfInputArrays; ++i) {
      vtkDataArray *arr = inputPointData->GetArray(i);

      if(arr and arr->GetNumberOfComponents() == 1) {
        vtkDataArray *newArr = arr->NewInstance();
        if(newArr == nullptr) {
          continue;
        }
        newArr->SetName(arr->GetName());
        newArr->SetNumberOfComponents(1);

        for(SimplexId v : vertices)
          newArr->InsertNextTuple1(arr->GetTuple1(v));

        outputPointData->AddArray(newArr);
      }
    }
  }

  {
    this->printMsg("Complete", 1, timer.getElapsedTime());
    this->printMsg(ttk::debug::Separator::L1);
  }
  return 1;
}
