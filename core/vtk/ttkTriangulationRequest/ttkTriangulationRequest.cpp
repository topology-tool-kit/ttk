#include <ttkMacros.h>
#include <ttkTriangulationRequest.h>
#include <ttkUtils.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkTriangulationRequest);

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

int ttkTriangulationRequest::RequestData(vtkInformation *request,
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  ttk::Timer timer;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector);

  const auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  const int dimensionality = triangulation->getDimensionality();
  const Request requestType = static_cast<Request>(RequestType);
  const Simplex simplexType = static_cast<Simplex>(SimplexType);

  vtkNew<vtkPoints> points{};
  vtkNew<vtkUnstructuredGrid> cells{};

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

    float p[3];
    triangulation->getVertexPoint(vertexId, p[0], p[1], p[2]);
    vertices.push_back(vertexId);
    return points->InsertNextPoint(p);
  };

  auto addEdge = [&](const SimplexId edgeId) {
    vtkIdType pointIds[2];
    SimplexId vertexIds[2];

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

    cells->InsertNextCell(VTK_LINE, 2, pointIds);

    return 0;
  };

  auto addTriangle = [&](const SimplexId triangleId) {
    vtkIdType pointIds[3];
    SimplexId vertexIds[3];

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

    cells->InsertNextCell(VTK_TRIANGLE, 3, pointIds);

    return 0;
  };

  auto addTetra = [&](const SimplexId tetraId) {
    vtkIdType pointIds[4];
    SimplexId vertexIds[4];

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

    cells->InsertNextCell(VTK_TETRA, 4, pointIds);

    return 0;
  };

  auto addStar = [&](const SimplexId starId) {
    if(dimensionality == 2)
      addTriangle(starId);
    else if(dimensionality == 3)
      addTetra(starId);
  };

  // do minimum preprocess and put watchdog on SimplexIdentifier
  switch(simplexType) {
    case Vertex:
      if(SimplexIdentifier < 0 or SimplexIdentifier >= numberOfVertices) {
        this->printErr("Vertex ID beyond the range.");
        return 0;
      }
      break;

    case Edge:
      triangulation->preconditionEdges();
      if(SimplexIdentifier < 0
         or SimplexIdentifier >= triangulation->getNumberOfEdges()) {
        this->printErr("Edge ID beyond the range.");
        return 0;
      }
      break;

    case Triangle:
      if(dimensionality == 2) {
        if(SimplexIdentifier < 0
           or SimplexIdentifier >= triangulation->getNumberOfCells()) {
          this->printErr("Triangle ID beyond the range.");
          return 0;
        }
      } else if(dimensionality == 3) {
        triangulation->preconditionTriangles();
        if(SimplexIdentifier < 0
           or SimplexIdentifier >= triangulation->getNumberOfTriangles()) {
          this->printErr("Triangle ID beyond the range.");
          return 0;
        }
      }
      break;

    case Tetra:
      if(dimensionality == 3)
        if(SimplexIdentifier < 0
           or SimplexIdentifier >= triangulation->getNumberOfCells()) {
          this->printErr("Tetrahedron ID beyond the range.");
          return 0;
        }
      break;
  }

  switch(requestType) {
    case ComputeSimplex:
      switch(simplexType) {
        case Vertex: {
          const auto vid = addVertex(SimplexIdentifier);
          cells->InsertNextCell(VTK_VERTEX, 1, &vid);
        } break;

        case Edge:
          addEdge(SimplexIdentifier);
          break;

        case Triangle:
          addTriangle(SimplexIdentifier);
          break;

        case Tetra:
          if(dimensionality == 3)
            addTetra(SimplexIdentifier);
          break;
      }
      break;

    case ComputeFacet:
      switch(simplexType) {
        case Vertex:
          break;

        case Edge:
          for(int i = 0; i < 2; ++i) {
            SimplexId vertexId;
            triangulation->getEdgeVertex(SimplexIdentifier, i, vertexId);
            addVertex(vertexId);
          }
          break;

        case Triangle:
          if(dimensionality == 2) {
            triangulation->preconditionCellEdges();
            for(int i = 0; i < 3; ++i) {
              SimplexId edgeId;
              triangulation->getCellEdge(SimplexIdentifier, i, edgeId);
              addEdge(edgeId);
            }
          } else if(dimensionality == 3) {
            triangulation->preconditionTriangleEdges();
            for(int i = 0; i < 3; ++i) {
              SimplexId edgeId;
              triangulation->getTriangleEdge(SimplexIdentifier, i, edgeId);
              addEdge(edgeId);
            }
          }
          break;

        case Tetra:
          if(dimensionality == 3) {
            triangulation->preconditionCellTriangles();
            for(int i = 0; i < 4; ++i) {
              SimplexId triangleId;
              triangulation->getCellTriangle(SimplexIdentifier, i, triangleId);
              addTriangle(triangleId);
            }
          }
          break;
      }
      break;

    case ComputeCofacet:
      switch(simplexType) {
        case Vertex:
          triangulation->preconditionVertexNeighbors();
          triangulation->preconditionVertexEdges();
          {
            const SimplexId edgeNumber
              = triangulation->getVertexEdgeNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < edgeNumber; ++i) {
              SimplexId edgeId;
              triangulation->getVertexEdge(SimplexIdentifier, i, edgeId);
              addEdge(edgeId);
            }
          }
          break;

        case Edge:
          if(dimensionality == 2) {
            triangulation->preconditionEdgeStars();
            const SimplexId starNumber
              = triangulation->getEdgeStarNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < starNumber; ++i) {
              SimplexId starId;
              triangulation->getEdgeStar(SimplexIdentifier, i, starId);
              addStar(starId);
            }
          } else if(dimensionality == 3) {
            triangulation->preconditionEdgeTriangles();
            const SimplexId triangleNumber
              = triangulation->getEdgeTriangleNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < triangleNumber; ++i) {
              SimplexId triangleId;
              triangulation->getEdgeTriangle(SimplexIdentifier, i, triangleId);
              addTriangle(triangleId);
            }
          }
          break;

        case Triangle:
          if(dimensionality == 3) {
            triangulation->preconditionTriangleStars();
            const SimplexId starNumber
              = triangulation->getTriangleStarNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < starNumber; ++i) {
              SimplexId starId;
              triangulation->getTriangleStar(SimplexIdentifier, i, starId);
              addStar(starId);
            }
          }
          break;

        case Tetra:
          break;
      }
      break;

    case ComputeStar:
      switch(simplexType) {
        case Vertex:
          triangulation->preconditionVertexStars();
          {
            const SimplexId starNumber
              = triangulation->getVertexStarNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < starNumber; ++i) {
              SimplexId starId;
              triangulation->getVertexStar(SimplexIdentifier, i, starId);
              addStar(starId);
            }
          }
          break;

        case Edge:
          triangulation->preconditionEdgeStars();
          {
            const SimplexId starNumber
              = triangulation->getEdgeStarNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < starNumber; ++i) {
              SimplexId starId;
              triangulation->getEdgeStar(SimplexIdentifier, i, starId);
              addStar(starId);
            }
          }
          break;

        case Triangle:
          if(dimensionality == 3) {
            triangulation->preconditionTriangleStars();
            const SimplexId starNumber
              = triangulation->getTriangleStarNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < starNumber; ++i) {
              SimplexId starId;
              triangulation->getTriangleStar(SimplexIdentifier, i, starId);
              addStar(starId);
            }
          }
          break;

        case Tetra:
          break;
      }
      break;

    case ComputeLink:
      switch(simplexType) {
        case Vertex:
          triangulation->preconditionVertexLinks();
          if(dimensionality == 2) {
            triangulation->preconditionEdges();
            const SimplexId linkNumber
              = triangulation->getVertexLinkNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < linkNumber; ++i) {
              SimplexId linkId;
              triangulation->getVertexLink(SimplexIdentifier, i, linkId);
              addEdge(linkId);
            }
          } else if(dimensionality == 3) {
            triangulation->preconditionVertexTriangles();
            const SimplexId linkNumber
              = triangulation->getVertexLinkNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < linkNumber; ++i) {
              SimplexId linkId;
              triangulation->getVertexLink(SimplexIdentifier, i, linkId);
              addTriangle(linkId);
            }
          }
          break;

        case Edge:
          triangulation->preconditionEdgeLinks();
          if(dimensionality == 2) {
            const SimplexId linkNumber
              = triangulation->getEdgeLinkNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < linkNumber; ++i) {
              SimplexId linkId;
              triangulation->getEdgeLink(SimplexIdentifier, i, linkId);
              addVertex(linkId);
            }
          } else if(dimensionality == 3) {
            const SimplexId linkNumber
              = triangulation->getEdgeLinkNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < linkNumber; ++i) {
              SimplexId linkId;
              triangulation->getEdgeLink(SimplexIdentifier, i, linkId);
              addEdge(linkId);
            }
          }
          break;

        case Triangle:
          if(dimensionality == 3) {
            triangulation->preconditionTriangleLinks();
            const SimplexId linkNumber
              = triangulation->getTriangleLinkNumber(SimplexIdentifier);
            for(SimplexId i = 0; i < linkNumber; ++i) {
              SimplexId linkId;
              triangulation->getTriangleLink(SimplexIdentifier, i, linkId);
              addVertex(linkId);
            }
          }
          break;

        case Tetra:
          break;
      }
      break;
  }

  cells->SetPoints(points);

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
