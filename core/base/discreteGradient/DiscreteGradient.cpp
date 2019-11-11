#include <DiscreteGradient.h>

using namespace std;
using namespace ttk;
using namespace dcg;

DiscreteGradient::DiscreteGradient()
  : IterationThreshold{-1}, ReverseSaddleMaximumConnection{false},
    ReverseSaddleSaddleConnection{false}, CollectPersistencePairs{false},
    ReturnSaddleConnectors{false}, SaddleConnectorsPersistenceThreshold{0},

    dimensionality_{-1}, numberOfVertices_{0}, gradient_{}, dmtMax2PL_{},

    inputScalarField_{nullptr}, inputOffsets_{nullptr},
    inputTriangulation_{nullptr},

    outputCriticalPoints_numberOfPoints_{nullptr},
    outputCriticalPoints_points_{nullptr},
    outputCriticalPoints_points_cellDimensions_{nullptr},
    outputCriticalPoints_points_cellIds_{nullptr},
    outputCriticalPoints_points_cellScalars_{nullptr},
    outputCriticalPoints_points_isOnBoundary_{nullptr},
    outputCriticalPoints_points_PLVertexIdentifiers_{nullptr},
    outputCriticalPoints_points_manifoldSize_{nullptr},

    outputGradientGlyphs_numberOfPoints_{nullptr},
    outputGradientGlyphs_points_{nullptr},
    outputGradientGlyphs_points_pairOrigins_{nullptr},
    outputGradientGlyphs_numberOfCells_{nullptr},
    outputGradientGlyphs_cells_{nullptr},
    outputGradientGlyphs_cells_pairTypes_{nullptr},

    outputPersistencePairs_{nullptr} {
}

DiscreteGradient::~DiscreteGradient() {
}

int DiscreteGradient::getDimensionality() const {
  return dimensionality_;
}

int DiscreteGradient::getNumberOfDimensions() const {
  return dimensionality_ + 1;
}

SimplexId DiscreteGradient::getNumberOfCells(const int dimension) const {
  if(dimensionality_ == 2) {
    switch(dimension) {
      case 0:
        return inputTriangulation_->getNumberOfVertices();
        break;

      case 1:
        return inputTriangulation_->getNumberOfEdges();
        break;

      case 2:
        return inputTriangulation_->getNumberOfCells();
        break;
    }
  } else if(dimensionality_ == 3) {
    switch(dimension) {
      case 0:
        return inputTriangulation_->getNumberOfVertices();
        break;

      case 1:
        return inputTriangulation_->getNumberOfEdges();
        break;

      case 2:
        return inputTriangulation_->getNumberOfTriangles();
        break;

      case 3:
        return inputTriangulation_->getNumberOfCells();
        break;
    }
  }

  return -1;
}

bool DiscreteGradient::isMinimum(const Cell &cell) const {
  if(cell.dim_ == 0)
    return (gradient_[0][0][cell.id_] == -1);

  return false;
}

bool DiscreteGradient::isSaddle1(const Cell &cell) const {
  if(cell.dim_ == 1)
    return (gradient_[0][1][cell.id_] == -1
            and gradient_[1][1][cell.id_] == -1);

  return false;
}

bool DiscreteGradient::isSaddle2(const Cell &cell) const {
  if(dimensionality_ == 3 and cell.dim_ == 2)
    return (gradient_[1][2][cell.id_] == -1
            and gradient_[2][2][cell.id_] == -1);

  return false;
}

bool DiscreteGradient::isMaximum(const Cell &cell) const {
  if(dimensionality_ == 2 and cell.dim_ == 2)
    return (gradient_[1][2][cell.id_] == -1);

  if(dimensionality_ == 3 and cell.dim_ == 3)
    return (gradient_[2][3][cell.id_] == -1);

  return false;
}

bool DiscreteGradient::isCellCritical(const int cellDim,
                                      const SimplexId cellId) const {
  if(dimensionality_ == 2) {
    switch(cellDim) {
      case 0:
        return (gradient_[0][0][cellId] == -1);
        break;

      case 1:
        return (gradient_[0][1][cellId] == -1
                and gradient_[1][1][cellId] == -1);
        break;

      case 2:
        return (gradient_[1][2][cellId] == -1);
        break;
    }
  } else if(dimensionality_ == 3) {
    switch(cellDim) {
      case 0:
        return (gradient_[0][0][cellId] == -1);
        break;

      case 1:
        return (gradient_[0][1][cellId] == -1
                and gradient_[1][1][cellId] == -1);
        break;

      case 2:
        return (gradient_[1][2][cellId] == -1
                and gradient_[2][2][cellId] == -1);
        break;

      case 3:
        return (gradient_[2][3][cellId] == -1);
        break;
    }
  }
  return false;
}

bool DiscreteGradient::isCellCritical(const Cell &cell) const {
  return isCellCritical(cell.dim_, cell.id_);
}

bool DiscreteGradient::isBoundary(const Cell &cell) const {
  const int cellDim = cell.dim_;
  const SimplexId cellId = cell.id_;

  if(dimensionality_ == 2) {
    switch(cellDim) {
      case 0:
        return inputTriangulation_->isVertexOnBoundary(cellId);

      case 1:
        return inputTriangulation_->isEdgeOnBoundary(cellId);

      case 2:
        for(int i = 0; i < 3; ++i) {
          SimplexId edgeId;
          inputTriangulation_->getCellEdge(cellId, i, edgeId);
          if(inputTriangulation_->isEdgeOnBoundary(edgeId))
            return true;
        }
        break;
    }
  } else if(dimensionality_ == 3) {
    switch(cellDim) {
      case 0:
        return inputTriangulation_->isVertexOnBoundary(cellId);

      case 1:
        return inputTriangulation_->isEdgeOnBoundary(cellId);

      case 2:
        return inputTriangulation_->isTriangleOnBoundary(cellId);

      case 3:
        for(int i = 0; i < 4; ++i) {
          SimplexId triangleId;
          inputTriangulation_->getCellTriangle(cellId, i, triangleId);
          if(inputTriangulation_->isTriangleOnBoundary(triangleId))
            return true;
        }
        break;
    }
  }
  return false;
}

SimplexId DiscreteGradient::getPairedCell(const Cell &cell,
                                          bool isReverse) const {
  SimplexId id{-1};
  if(dimensionality_ == 2) {
    switch(cell.dim_) {
      case 0:
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        inputTriangulation_->getVertexEdge(
          cell.id_, gradient_[0][0][cell.id_], id);
#else
        return gradient_[0][0][cell.id_];
#endif
        break;

      case 1:
        if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
          inputTriangulation_->getEdgeVertex(
            cell.id_, gradient_[0][1][cell.id_], id);
          return id;
#else
          return gradient_[0][1][cell.id_];
#endif
        }

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        inputTriangulation_->getEdgeStar(
          cell.id_, gradient_[1][1][cell.id_], id);
#else
        return gradient_[1][1][cell.id_];
#endif
        break;

      case 2:
        if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
          inputTriangulation_->getCellEdge(
            cell.id_, gradient_[1][2][cell.id_], id);
          return id;
#else
          return gradient_[1][2][cell.id_];
#endif
        }
        break;
    }
  } else if(dimensionality_ == 3) {
    switch(cell.dim_) {
      case 0:
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        inputTriangulation_->getVertexEdge(
          cell.id_, gradient_[0][0][cell.id_], id);
#else
        return gradient_[0][0][cell.id_];
#endif
        break;

      case 1:
        if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
          inputTriangulation_->getEdgeVertex(
            cell.id_, gradient_[0][1][cell.id_], id);
          return id;
#else
          return gradient_[0][1][cell.id_];
#endif
        }

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        inputTriangulation_->getEdgeTriangle(
          cell.id_, gradient_[1][1][cell.id_], id);
#else
        return gradient_[1][1][cell.id_];
#endif
        break;

      case 2:
        if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
          inputTriangulation_->getTriangleEdge(
            cell.id_, gradient_[1][2][cell.id_], id);
          return id;
#else
          return gradient_[1][2][cell.id_];
#endif
        }

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        inputTriangulation_->getTriangleStar(
          cell.id_, gradient_[2][2][cell.id_], id);
#else
        return gradient_[2][2][cell.id_];
#endif
        break;

      case 3:
        if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
          inputTriangulation_->getCellTriangle(
            cell.id_, gradient_[2][3][cell.id_], id);
          return id;
#else
          return gradient_[2][3][cell.id_];
#endif
        }
        break;
    }
  }

  return id;
}

int DiscreteGradient::getCriticalPoints(vector<Cell> &criticalPoints) const {
  // foreach dimension
  const int numberOfDimensions = getNumberOfDimensions();
  for(int i = 0; i < numberOfDimensions; ++i) {

    // foreach cell of that dimension
    const SimplexId numberOfCells = getNumberOfCells(i);
    for(SimplexId j = 0; j < numberOfCells; ++j) {
      const Cell cell(i, j);

      if(isCellCritical(cell))
        criticalPoints.push_back(cell);
    }
  }

  return 0;
}

int DiscreteGradient::getDescendingPath(const Cell &cell,
                                        vector<Cell> &vpath) const {
  if(dimensionality_ == 2) {
    if(cell.dim_ == 0) {
      // assume that cellId is a vertex
      SimplexId currentId = cell.id_;
      SimplexId connectedEdgeId;
      do {
        // add a vertex
        const Cell vertex(0, currentId);
        vpath.push_back(std::move(vertex));

        if(isCellCritical(vertex))
          break;

        connectedEdgeId = getPairedCell(vertex);
        if(connectedEdgeId == -1)
          break;

        // add an edge
        const Cell edge(1, connectedEdgeId);
        vpath.push_back(std::move(edge));

        if(isCellCritical(edge))
          break;

        for(int i = 0; i < 2; ++i) {
          SimplexId vertexId;
          inputTriangulation_->getEdgeVertex(connectedEdgeId, i, vertexId);

          if(vertexId != currentId) {
            currentId = vertexId;
            break;
          }
        }

      } while(connectedEdgeId != -1);
    }
  } else if(dimensionality_ == 3) {
    if(cell.dim_ == 0) {
      // assume that cellId is a vertex
      SimplexId currentId = cell.id_;
      SimplexId connectedEdgeId;
      do {
        // add a vertex
        const Cell vertex(0, currentId);
        vpath.push_back(std::move(vertex));

        if(isCellCritical(vertex))
          break;

        connectedEdgeId = getPairedCell(vertex);
        if(connectedEdgeId == -1)
          break;

        // add an edge
        const Cell edge(1, connectedEdgeId);
        vpath.push_back(std::move(edge));

        if(isCellCritical(edge))
          break;

        for(int i = 0; i < 2; ++i) {
          SimplexId vertexId;
          inputTriangulation_->getEdgeVertex(connectedEdgeId, i, vertexId);

          if(vertexId != currentId) {
            currentId = vertexId;
            break;
          }
        }

      } while(connectedEdgeId != -1);
    }
  }

  return 0;
}

int DiscreteGradient::getDescendingPathThroughWall(
  const wallId_t wallId,
  const Cell &saddle2,
  const Cell &saddle1,
  const vector<wallId_t> &isVisited,
  vector<Cell> *const vpath,
  const bool enableCycleDetector) const {
  // debug
  const SimplexId numberOfEdges = inputTriangulation_->getNumberOfEdges();
  vector<char> isCycle;
  if(enableCycleDetector)
    isCycle.resize(numberOfEdges, false);

  if(dimensionality_ == 3) {
    // add the 2-saddle to the path
    if(vpath)
      vpath->push_back(saddle2);

    SimplexId currentId = -1;
    {
      int nconnections = 0;
      for(int i = 0; i < 3; ++i) {
        SimplexId edgeId;
        inputTriangulation_->getTriangleEdge(saddle2.id_, i, edgeId);
        if(isVisited[edgeId] == wallId) {
          // saddle2 can be adjacent to saddle1 on the wall
          if(isSaddle1(Cell(1, edgeId))) {
            if(vpath)
              vpath->push_back(Cell(1, edgeId));
            return false;
          }

          currentId = edgeId;
          ++nconnections;
        }
      }
      if(nconnections > 1)
        return true;
    }

    int oldId;
    do {

      // debug
      if(enableCycleDetector) {
        if(!isCycle[currentId]) {
          isCycle[currentId] = true;
        } else {
          cout << "[DiscreteGradient] Error : cycle detected on the wall of "
                  "1-saddleId="
               << saddle1.id_ << endl;
          break;
        }
      }

      oldId = currentId;

      // add an edge
      const Cell edge(1, currentId);
      if(vpath)
        vpath->push_back(edge);

      if(isCellCritical(edge))
        break;

      const SimplexId connectedTriangleId = getPairedCell(edge);

      // add a triangle
      const Cell triangle(2, connectedTriangleId);
      if(vpath)
        vpath->push_back(triangle);

      if(isCellCritical(triangle))
        break;

      int nconnections = 0;
      for(int i = 0; i < 3; ++i) {
        SimplexId edgeId;
        inputTriangulation_->getTriangleEdge(connectedTriangleId, i, edgeId);

        if(isVisited[edgeId] == wallId and edgeId != oldId) {
          currentId = edgeId;
          ++nconnections;
        }
      }
      if(nconnections > 1)
        return true;

      // stop at convergence caused by boundary effect
    } while(currentId != oldId);
  }

  return false;
}

int DiscreteGradient::getAscendingPath(const Cell &cell,
                                       vector<Cell> &vpath,
                                       const bool enableCycleDetector) const {

  const SimplexId numberOfCells = inputTriangulation_->getNumberOfCells();
  vector<char> isCycle;
  if(enableCycleDetector)
    isCycle.resize(numberOfCells, false);

  if(dimensionality_ == 2) {
    if(cell.dim_ == 2) {
      // assume that cellId is a triangle
      SimplexId currentId = cell.id_;
      SimplexId oldId;
      do {
        oldId = currentId;

        // add a triangle
        const Cell triangle(2, currentId);
        vpath.push_back(triangle);

        if(isCellCritical(triangle))
          break;

        const SimplexId connectedEdgeId = getPairedCell(triangle, true);
        if(connectedEdgeId == -1)
          break;

        // add an edge
        const Cell edge(1, connectedEdgeId);
        vpath.push_back(edge);

        if(isCellCritical(edge))
          break;

        const SimplexId starNumber
          = inputTriangulation_->getEdgeStarNumber(connectedEdgeId);
        for(SimplexId i = 0; i < starNumber; ++i) {
          SimplexId starId;
          inputTriangulation_->getEdgeStar(connectedEdgeId, i, starId);

          if(starId != currentId) {
            currentId = starId;
            break;
          }
        }

        // stop at convergence caused by boundary effect
      } while(currentId != oldId);
    }
  } else if(dimensionality_ == 3) {
    if(cell.dim_ == 3) {
      // assume that cellId is a tetra
      SimplexId currentId = cell.id_;
      SimplexId oldId;
      do {

        // debug
        if(enableCycleDetector) {
          if(!isCycle[currentId]) {
            isCycle[currentId] = true;
          } else {
            cout << "[DiscreteGradient] Error : cycle detected in the path "
                    "from tetraId="
                 << cell.id_ << endl;
            break;
          }
        }

        oldId = currentId;

        // add a tetra
        const Cell tetra(3, currentId);
        vpath.push_back(tetra);

        if(isCellCritical(tetra))
          break;

        const SimplexId connectedTriangleId = getPairedCell(tetra, true);
        if(connectedTriangleId == -1)
          break;

        // add a triangle
        const Cell triangle(2, connectedTriangleId);
        vpath.push_back(triangle);

        if(isCellCritical(triangle))
          break;

        const SimplexId starNumber
          = inputTriangulation_->getTriangleStarNumber(connectedTriangleId);
        for(SimplexId i = 0; i < starNumber; ++i) {
          SimplexId starId;
          inputTriangulation_->getTriangleStar(connectedTriangleId, i, starId);

          if(starId != currentId) {
            currentId = starId;
            break;
          }
        }

        // stop at convergence caused by boundary effect
      } while(currentId != oldId);
    }
  }

  return 0;
}

bool DiscreteGradient::getAscendingPathThroughWall(
  const wallId_t wallId,
  const Cell &saddle1,
  const Cell &saddle2,
  const vector<wallId_t> &isVisited,
  vector<Cell> *const vpath,
  const bool enableCycleDetector) const {
  // debug
  const SimplexId numberOfTriangles
    = inputTriangulation_->getNumberOfTriangles();
  vector<char> isCycle;
  if(enableCycleDetector)
    isCycle.resize(numberOfTriangles, false);

  if(dimensionality_ == 3) {
    // add the 1-saddle to the path
    if(vpath)
      vpath->push_back(saddle1);

    SimplexId currentId = -1;
    {
      int nconnections = 0;
      const SimplexId triangleNumber
        = inputTriangulation_->getEdgeTriangleNumber(saddle1.id_);
      for(SimplexId i = 0; i < triangleNumber; ++i) {
        SimplexId triangleId;
        inputTriangulation_->getEdgeTriangle(saddle1.id_, i, triangleId);
        if(isVisited[triangleId] == wallId) {
          // saddle1 can be adjacent to saddle2 on the wall
          if(isSaddle2(Cell(2, triangleId))) {
            if(vpath)
              vpath->push_back(Cell(2, triangleId));
            return false;
          }

          currentId = triangleId;
          ++nconnections;
        }
      }
      if(nconnections > 1)
        return true;
    }

    SimplexId oldId;
    do {

      // debug
      if(enableCycleDetector) {
        if(!isCycle[currentId]) {
          isCycle[currentId] = true;
        } else {
          cout << "[DiscreteGradient] Error : cycle detected on the wall of "
                  "2-saddleId="
               << saddle2.id_ << endl;
          break;
        }
      }

      oldId = currentId;

      // add a triangle
      const Cell triangle(2, currentId);
      if(vpath)
        vpath->push_back(triangle);

      if(isCellCritical(triangle))
        break;

      const SimplexId connectedEdgeId = getPairedCell(triangle, true);

      // add an edge
      const Cell edge(1, connectedEdgeId);
      if(vpath)
        vpath->push_back(edge);

      if(isCellCritical(edge))
        break;

      int nconnections = 0;
      const SimplexId triangleNumber
        = inputTriangulation_->getEdgeTriangleNumber(connectedEdgeId);
      for(SimplexId i = 0; i < triangleNumber; ++i) {
        SimplexId triangleId;
        inputTriangulation_->getEdgeTriangle(connectedEdgeId, i, triangleId);

        if(isVisited[triangleId] == wallId and triangleId != oldId) {
          currentId = triangleId;
          ++nconnections;
        }
      }
      if(nconnections > 1)
        return true;

      // stop at convergence caused by boundary effect
    } while(currentId != oldId);
  }

  return false;
}

int DiscreteGradient::getDescendingWall(const wallId_t wallId,
                                        const Cell &cell,
                                        vector<wallId_t> &isVisited,
                                        vector<Cell> *const wall,
                                        set<SimplexId> *const saddles) const {
  if(dimensionality_ == 3) {
    if(cell.dim_ == 2) {
      // assume that cellId is a triangle
      const SimplexId originId = cell.id_;

      queue<SimplexId> bfs;
      bfs.push(originId);

      // BFS traversal
      while(!bfs.empty()) {
        const SimplexId triangleId = bfs.front();
        bfs.pop();

        if(isVisited[triangleId] != wallId) {
          isVisited[triangleId] = wallId;

          // add the triangle
          if(wall)
            wall->push_back(Cell(2, triangleId));

          for(int j = 0; j < 3; ++j) {
            SimplexId edgeId;
            inputTriangulation_->getTriangleEdge(triangleId, j, edgeId);

            if(saddles and isSaddle1(Cell(1, edgeId)))
              saddles->insert(edgeId);

            const SimplexId pairedCellId = getPairedCell(Cell(1, edgeId));

            if(pairedCellId != -1 and pairedCellId != triangleId)
              bfs.push(pairedCellId);
          }
        }
      }
    }
  }

  return 0;
}

int DiscreteGradient::getAscendingWall(const wallId_t wallId,
                                       const Cell &cell,
                                       vector<wallId_t> &isVisited,
                                       vector<Cell> *const wall,
                                       set<SimplexId> *const saddles) const {
  if(dimensionality_ == 3) {
    if(cell.dim_ == 1) {
      // assume that cellId is an edge
      const SimplexId originId = cell.id_;

      queue<SimplexId> bfs;
      bfs.push(originId);

      // BFS traversal
      while(!bfs.empty()) {
        const SimplexId edgeId = bfs.front();
        bfs.pop();

        if(isVisited[edgeId] != wallId) {
          isVisited[edgeId] = wallId;

          // add the edge
          if(wall)
            wall->push_back(Cell(1, edgeId));

          const SimplexId triangleNumber
            = inputTriangulation_->getEdgeTriangleNumber(edgeId);
          for(SimplexId j = 0; j < triangleNumber; ++j) {
            SimplexId triangleId;
            inputTriangulation_->getEdgeTriangle(edgeId, j, triangleId);

            if(saddles and isSaddle2(Cell(2, triangleId)))
              saddles->insert(triangleId);

            const SimplexId pairedCellId
              = getPairedCell(Cell(2, triangleId), true);

            if(pairedCellId != -1 and pairedCellId != edgeId)
              bfs.push(pairedCellId);
          }
        }
      }
    }
  }

  return 0;
}

int DiscreteGradient::reverseAscendingPath(const vector<Cell> &vpath) {
  if(dimensionality_ == 2) {
    // assume that the first cell is an edge
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId edgeId = vpath[i].id_;
      const SimplexId triangleId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 3; ++k) {
        SimplexId tmp;
        inputTriangulation_->getCellEdge(triangleId, k, tmp);
        if(tmp == edgeId) {
          gradient_[1][2][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < inputTriangulation_->getEdgeStarNumber(edgeId); ++k) {
        SimplexId tmp;
        inputTriangulation_->getEdgeStar(edgeId, k, tmp);
        if(tmp == triangleId) {
          gradient_[1][1][edgeId] = k;
          break;
        }
      }
#else
      gradient_[1][2][triangleId] = edgeId;
      gradient_[1][1][edgeId] = triangleId;
#endif
    }
  } else if(dimensionality_ == 3) {
    // assume that the first cell is a triangle
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId triangleId = vpath[i].id_;
      const SimplexId tetraId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 4; ++k) {
        SimplexId tmp;
        inputTriangulation_->getCellTriangle(tetraId, k, tmp);
        if(tmp == triangleId) {
          gradient_[2][3][tetraId] = k;
          break;
        }
      }
      for(int k = 0; k < inputTriangulation_->getTriangleStarNumber(triangleId);
          ++k) {
        SimplexId tmp;
        inputTriangulation_->getTriangleStar(triangleId, k, tmp);
        if(tmp == tetraId) {
          gradient_[2][2][triangleId] = k;
          break;
        }
      }
#else
      gradient_[2][3][tetraId] = triangleId;
      gradient_[2][2][triangleId] = tetraId;
#endif
    }
  }

  return 0;
}

int DiscreteGradient::reverseAscendingPathOnWall(const vector<Cell> &vpath) {
  if(dimensionality_ == 3) {
    // assume that the first cell is an edge
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId edgeId = vpath[i].id_;
      const SimplexId triangleId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 3; ++k) {
        SimplexId tmp;
        inputTriangulation_->getTriangleEdge(triangleId, k, tmp);
        if(tmp == edgeId) {
          gradient_[1][2][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < inputTriangulation_->getEdgeTriangleNumber(edgeId);
          ++k) {
        SimplexId tmp;
        inputTriangulation_->getEdgeTriangle(edgeId, k, tmp);
        if(tmp == triangleId) {
          gradient_[1][1][edgeId] = k;
          break;
        }
      }
#else
      gradient_[1][2][triangleId] = edgeId;
      gradient_[1][1][edgeId] = triangleId;
#endif
    }
  }

  return 0;
}

int DiscreteGradient::reverseDescendingPathOnWall(const vector<Cell> &vpath) {
  if(dimensionality_ == 3) {
    // assume that the first cell is a triangle
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId triangleId = vpath[i].id_;
      const SimplexId edgeId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 3; ++k) {
        SimplexId tmp;
        inputTriangulation_->getTriangleEdge(triangleId, k, tmp);
        if(tmp == edgeId) {
          gradient_[1][2][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < inputTriangulation_->getEdgeTriangleNumber(edgeId);
          ++k) {
        SimplexId tmp;
        inputTriangulation_->getEdgeTriangle(edgeId, k, tmp);
        if(tmp == triangleId) {
          gradient_[1][1][edgeId] = k;
          break;
        }
      }
#else
      gradient_[1][1][edgeId] = triangleId;
      gradient_[1][2][triangleId] = edgeId;
#endif
    }
  }

  return 0;
}

int DiscreteGradient::getEdgeIncenter(const SimplexId edgeId,
                                      float incenter[3]) const {
  return inputTriangulation_->getEdgeIncenter(edgeId, incenter);
}

int DiscreteGradient::getTriangleIncenter(const SimplexId triangleId,
                                          float incenter[3]) const {
  return inputTriangulation_->getTriangleIncenter(triangleId, incenter);
}

int DiscreteGradient::getTetraIncenter(const SimplexId tetraId,
                                       float incenter[3]) const {
  return inputTriangulation_->getTetraIncenter(tetraId, incenter);
}

int DiscreteGradient::getCellIncenter(const Cell &cell,
                                      float incenter[3]) const {
  return inputTriangulation_->getCellIncenter(cell.id_, cell.dim_, incenter);
}

int DiscreteGradient::getCriticalPointMap(
  const vector<pair<SimplexId, char>> &criticalPoints, vector<char> &isPL) {
  isPL.resize(numberOfVertices_);
  std::fill(isPL.begin(), isPL.end(), 0);
  for(pair<SimplexId, char> criticalPoint : criticalPoints) {
    const SimplexId criticalPointId = criticalPoint.first;
    const char criticalPointType = criticalPoint.second;

    isPL[criticalPointId] = criticalPointType;
  }

  return 0;
}

int DiscreteGradient::setGradientGlyphs() const {
  (*outputGradientGlyphs_numberOfPoints_) = 0;
  (*outputGradientGlyphs_numberOfCells_) = 0;

  SimplexId pointId{};
  SimplexId cellId{};

  // foreach dimension
  const int numberOfDimensions = getNumberOfDimensions();
  for(int i = 0; i < numberOfDimensions - 1; ++i) {
    // foreach cell of that dimension
    const SimplexId numberOfCells = getNumberOfCells(i);
    for(SimplexId j = 0; j < numberOfCells; ++j) {
      const Cell cell(i, j);

      const SimplexId pairedCellId = getPairedCell(cell);
      if(pairedCellId != -1) {
        // get gradient pair
        const int pairedCellDim = i + 1;

        const Cell pairedCell(pairedCellDim, pairedCellId);

        float p0[3];
        getCellIncenter(cell, p0);

        outputGradientGlyphs_points_->push_back(p0[0]);
        outputGradientGlyphs_points_->push_back(p0[1]);
        outputGradientGlyphs_points_->push_back(p0[2]);

        outputGradientGlyphs_points_pairOrigins_->push_back(false);

        float p1[3];
        getCellIncenter(pairedCell, p1);

        outputGradientGlyphs_points_->push_back(p1[0]);
        outputGradientGlyphs_points_->push_back(p1[1]);
        outputGradientGlyphs_points_->push_back(p1[2]);

        outputGradientGlyphs_points_pairOrigins_->push_back(true);

        outputGradientGlyphs_cells_->push_back(2);
        outputGradientGlyphs_cells_->push_back(pointId);
        outputGradientGlyphs_cells_->push_back(pointId + 1);

        outputGradientGlyphs_cells_pairTypes_->push_back(i);

        pointId += 2;
        cellId += 1;
      }
    }
  }

  (*outputGradientGlyphs_numberOfPoints_) = pointId;
  (*outputGradientGlyphs_numberOfCells_) = cellId;

  return 0;
}
