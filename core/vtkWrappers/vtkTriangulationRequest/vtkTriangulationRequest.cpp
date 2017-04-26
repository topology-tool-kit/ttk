#include                  <vtkTriangulationRequest.h>

vtkStandardNewMacro(vtkTriangulationRequest)

  int vtkTriangulationRequest::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

    Memory m;

    vtkDataSet* input=inputs[0];
    vtkUnstructuredGrid* output=vtkUnstructuredGrid::SafeDownCast(outputs[0]);

    Triangulation *triangulation=vtkTriangulation::getTriangulation(input);

    if(!triangulation)
      return -1;

    triangulation->setWrapper(this);
    const int dimensionality=triangulation->getDimensionality();
    const Request requestType=static_cast<Request>(RequestType);
    const Simplex simplexType=static_cast<Simplex>(SimplexType);

    vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnstructuredGrid> cells=vtkSmartPointer<vtkUnstructuredGrid>::New();

    const int numberOfVertices=triangulation->getNumberOfVertices();
    vector<vtkIdType> isVisited(numberOfVertices, -1);

    cells->Allocate();

    auto addVertex=[&](const int vertexId){
      if(vertexId==-1) return vtkIdType(-1);

      float p[3];
      triangulation->getVertexPoint(vertexId, p[0], p[1], p[2]);
      return points->InsertNextPoint(p);
    };

    auto addEdge=[&](const int edgeId){
      vtkIdType pointIds[2];
      int vertexIds[2];

      for(int i=0; i<2; ++i){
        triangulation->getEdgeVertex(edgeId, i, vertexIds[i]);

        const int vertexId=vertexIds[i];
        if(vertexId==-1) return -1;

        if(isVisited[vertexId]==-1){
          pointIds[i]=addVertex(vertexId);
          isVisited[vertexId]=pointIds[i];
        }
        else
          pointIds[i]=isVisited[vertexId];
      }

      cells->InsertNextCell(VTK_LINE, 2, pointIds);

      return 0;
    };

    auto addTriangle=[&](const int triangleId){
      vtkIdType pointIds[3];
      int vertexIds[3];

      for(int i=0; i<3; ++i){
        if(dimensionality==3)
          triangulation->getTriangleVertex(triangleId, i, vertexIds[i]);
        else
          triangulation->getCellVertex(triangleId, i, vertexIds[i]);

        const int vertexId=vertexIds[i];
        if(vertexId==-1) return -1;

        if(isVisited[vertexId]==-1){
          pointIds[i]=addVertex(vertexId);
          isVisited[vertexId]=pointIds[i];
        }
        else
          pointIds[i]=isVisited[vertexId];
      }

      cells->InsertNextCell(VTK_TRIANGLE, 3, pointIds);

      return 0;
    };

    auto addTetra=[&](const int tetraId){
      vtkIdType pointIds[4];
      int vertexIds[4];

      for(int i=0; i<4; ++i){
        triangulation->getCellVertex(tetraId, i, vertexIds[i]);

        const int vertexId=vertexIds[i];
        if(vertexId==-1) return -1;

        if(isVisited[vertexId]==-1){
          pointIds[i]=addVertex(vertexId);
          isVisited[vertexId]=pointIds[i];
        }
        else
          pointIds[i]=isVisited[vertexId];
      }

      cells->InsertNextCell(VTK_TETRA, 4, pointIds);

      return 0;
    };

    auto addStar=[&](const int starId){
      if(dimensionality==2)
        addTriangle(starId);
      else if(dimensionality==3)
        addTetra(starId);
    };

    // do minimum preprocess and put watchdog on SimplexId
    switch(simplexType){
      case Vertex:
        if(SimplexId<0 or SimplexId>=numberOfVertices) return -1;
        break;

      case Edge:
        triangulation->preprocessEdges();
        if(SimplexId<0 or SimplexId>=triangulation->getNumberOfEdges()) return -1;
        break;

      case Triangle:
        if(dimensionality==2){
          if(SimplexId<0 or SimplexId>=triangulation->getNumberOfCells()) return -1;
        }
        else if(dimensionality==3){
          triangulation->preprocessTriangles();
          if(SimplexId<0 or SimplexId>=triangulation->getNumberOfTriangles()) return -1;
        }
        break;

      case Tetra:
        if(dimensionality==3)
          if(SimplexId<0 or SimplexId>=triangulation->getNumberOfCells()) return -1;
        break;
    }

    switch(requestType){
      case ComputeSimplex:
        switch(simplexType){
          case Vertex:
            addVertex(SimplexId);
            break;

          case Edge:
            addEdge(SimplexId);
            break;

          case Triangle:
            addTriangle(SimplexId);
            break;

          case Tetra:
            if(dimensionality==3)
              addTetra(SimplexId);
            break;
        }
        break;

      case ComputeFacet:
        switch(simplexType){
          case Vertex:
            break;

          case Edge:
            for(int i=0; i<2; ++i){
              int vertexId;
              triangulation->getEdgeVertex(SimplexId, i, vertexId);
              addVertex(vertexId);
            }
            break;

          case Triangle:
            if(dimensionality==2){
              triangulation->preprocessCellEdges();
              for(int i=0; i<3; ++i){
                int edgeId;
                triangulation->getCellEdge(SimplexId, i, edgeId);
                addEdge(edgeId);
              }
            }
            else if(dimensionality==3){
              triangulation->preprocessTriangleEdges();
              for(int i=0; i<3; ++i){
                int edgeId;
                triangulation->getEdgeVertex(SimplexId, i, edgeId);
                addEdge(edgeId);
              }
            }
            break;

          case Tetra:
            if(dimensionality==3){
              triangulation->preprocessCellTriangles();
              for(int i=0; i<4; ++i){
                int triangleId;
                triangulation->getCellTriangle(SimplexId, i, triangleId);
                addTriangle(triangleId);
              }
            }
            break;
        }
        break;

      case ComputeCofacet:
        switch(simplexType){
          case Vertex:
            triangulation->preprocessVertexEdges();
            {
              const int edgeNumber=triangulation->getVertexEdgeNumber(SimplexId);
              for(int i=0; i<edgeNumber; ++i){
                int edgeId;
                triangulation->getVertexEdge(SimplexId, i, edgeId);
                addEdge(edgeId);
              }
            }
            break;

          case Edge:
            if(dimensionality==2){
              triangulation->preprocessEdgeStars();
              const int starNumber=triangulation->getEdgeStarNumber(SimplexId);
              for(int i=0; i<starNumber; ++i){
                int starId;
                triangulation->getEdgeStar(SimplexId, i, starId);
                addStar(starId);
              }
            }
            else if(dimensionality==3){
              triangulation->preprocessEdgeTriangles();
              const int triangleNumber=triangulation->getEdgeTriangleNumber(SimplexId);
              for(int i=0; i<triangleNumber; ++i){
                int triangleId;
                triangulation->getEdgeTriangle(SimplexId, i, triangleId);
                addTriangle(triangleId);
              }
            }
            break;

          case Triangle:
            if(dimensionality==3){
              triangulation->preprocessTriangleStars();
              const int starNumber=triangulation->getTriangleStarNumber(SimplexId);
              for(int i=0; i<starNumber; ++i){
                int starId;
                triangulation->getTriangleStar(SimplexId, i, starId);
                addStar(starId);
              }
            }
            break;

          case Tetra:
            break;
        }
        break;

      case ComputeStar:
        switch(simplexType){
          case Vertex:
            triangulation->preprocessVertexStars();
            {
              const int starNumber=triangulation->getVertexStarNumber(SimplexId);
              for(int i=0; i<starNumber; ++i){
                int starId;
                triangulation->getVertexStar(SimplexId, i, starId);
                addStar(starId);
              }
            }
            break;

          case Edge:
            triangulation->preprocessEdgeStars();
            {
              const int starNumber=triangulation->getEdgeStarNumber(SimplexId);
              for(int i=0; i<starNumber; ++i){
                int starId;
                triangulation->getEdgeStar(SimplexId, i, starId);
                addStar(starId);
              }
            }
            break;

          case Triangle:
            if(dimensionality==3){
              triangulation->preprocessTriangleStars();
              const int starNumber=triangulation->getTriangleStarNumber(SimplexId);
              for(int i=0; i<starNumber; ++i){
                int starId;
                triangulation->getTriangleStar(SimplexId, i, starId);
                addStar(starId);
              }
            }
            break;

          case Tetra:
            break;
        }
        break;

      case ComputeLink:
        switch(simplexType){
          case Vertex:
            triangulation->preprocessVertexLinks();
            if(dimensionality==2){
              const int linkNumber=triangulation->getVertexLinkNumber(SimplexId);
              for(int i=0; i<linkNumber; ++i){
                int linkId;
                triangulation->getVertexLink(SimplexId, i, linkId);
                addEdge(linkId);
              }
            }
            else if(dimensionality==3){
              const int linkNumber=triangulation->getVertexLinkNumber(SimplexId);
              for(int i=0; i<linkNumber; ++i){
                int linkId;
                triangulation->getVertexLink(SimplexId, i, linkId);
                addTriangle(linkId);
              }
            }
            break;

          case Edge:
            triangulation->preprocessEdgeLinks();
            if(dimensionality==2){
              const int linkNumber=triangulation->getEdgeLinkNumber(SimplexId);
              for(int i=0; i<linkNumber; ++i){
                int linkId;
                triangulation->getEdgeLink(SimplexId, i, linkId);
                addVertex(linkId);
              }
            }
            else if(dimensionality==3){
              const int linkNumber=triangulation->getEdgeLinkNumber(SimplexId);
              for(int i=0; i<linkNumber; ++i){
                int linkId;
                triangulation->getEdgeLink(SimplexId, i, linkId);
                addEdge(linkId);
              }
            }
            break;

          case Triangle:
            if(dimensionality==3){
              triangulation->preprocessTriangleLinks();
              const int linkNumber=triangulation->getTriangleLinkNumber(SimplexId);
              for(int i=0; i<linkNumber; ++i){
                int linkId;
                triangulation->getTriangleLink(SimplexId, i, linkId);
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

    {
      stringstream msg;
      msg << "[vtkTriangulationRequest] Memory usage: " << m.getElapsedUsage() 
        << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
    }

    return 0;
  }
