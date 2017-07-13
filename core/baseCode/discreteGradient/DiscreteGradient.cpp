#include                  <DiscreteGradient.h>

DiscreteGradient::DiscreteGradient():
  IterationThreshold{-1},
  ReverseSaddleMaximumConnection{},
  ReverseSaddleSaddleConnection{},
  CollectPersistencePairs{},

  dimensionality_{-1},
  gradient_{},

  inputScalarField_{},
  inputOffsets_{},
  inputTriangulation_{},

  outputCriticalPoints_numberOfPoints_{},
  outputCriticalPoints_points_{},
  outputCriticalPoints_points_cellDimensions_{},
  outputCriticalPoints_points_cellIds_{},
  outputCriticalPoints_points_cellScalars_{},
  outputCriticalPoints_points_isOnBoundary_{},

  outputGradientGlyphs_numberOfPoints_{},
  outputGradientGlyphs_points_{},
  outputGradientGlyphs_points_pairOrigins_{},
  outputGradientGlyphs_numberOfCells_{},
  outputGradientGlyphs_cells_{},
  outputGradientGlyphs_cells_pairTypes_{},

  outputPersistencePairs_{}
{}

DiscreteGradient::~DiscreteGradient(){
}

int DiscreteGradient::getDimensionality() const{
  return dimensionality_;
}

int DiscreteGradient::getNumberOfDimensions() const{
  return dimensionality_+1;
}

int DiscreteGradient::getNumberOfCells(const int dimension) const{
  if(dimensionality_==2){
    switch(dimension){
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
  }
  else if(dimensionality_==3){
    switch(dimension){
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

bool DiscreteGradient::isMinimum(const Cell& cell) const{
  if(cell.dim_==0)
    return (gradient_[0][0][cell.id_]==-1);

  return false;
}

bool DiscreteGradient::isSaddle1(const Cell& cell) const{
  if(cell.dim_==1)
    return (gradient_[0][1][cell.id_]==-1 and gradient_[1][1][cell.id_]==-1);

  return false;
}

bool DiscreteGradient::isSaddle2(const Cell& cell) const{
  if(dimensionality_==3 and cell.dim_==2)
    return (gradient_[1][2][cell.id_]==-1 and gradient_[2][2][cell.id_]==-1);

  return false;
}

bool DiscreteGradient::isMaximum(const Cell& cell) const{
  if(dimensionality_==2 and cell.dim_==2)
    return (gradient_[1][2][cell.id_]==-1);

  if(dimensionality_==3 and cell.dim_==3)
    return (gradient_[2][3][cell.id_]==-1);

  return false;
}

bool DiscreteGradient::isCellCritical(const Cell& cell) const{
  if(dimensionality_==2){
    switch(cell.dim_){
      case 0:
        return (gradient_[0][0][cell.id_]==-1);
        break;

      case 1:
        return (gradient_[0][1][cell.id_]==-1 and gradient_[1][1][cell.id_]==-1);
        break;

      case 2:
        return (gradient_[1][2][cell.id_]==-1);
        break;
    }
  }
  else if(dimensionality_==3){
    switch(cell.dim_){
      case 0:
        return (gradient_[0][0][cell.id_]==-1);
        break;

      case 1:
        return (gradient_[0][1][cell.id_]==-1 and gradient_[1][1][cell.id_]==-1);
        break;

      case 2:
        return (gradient_[1][2][cell.id_]==-1 and gradient_[2][2][cell.id_]==-1);
        break;

      case 3:
        return (gradient_[2][3][cell.id_]==-1);
        break;
    }
  }

  return false;
}

bool DiscreteGradient::isBoundary(const Cell& cell) const{
  const int cellDim=cell.dim_;
  const int cellId=cell.id_;

  if(dimensionality_==2){
    switch(cellDim){
      case 0:
        return inputTriangulation_->isVertexOnBoundary(cellId);

      case 1:
        return inputTriangulation_->isEdgeOnBoundary(cellId);

      case 2:
        for(int i=0; i<3; ++i){
          int edgeId;
          inputTriangulation_->getCellEdge(cellId, i, edgeId);
          if(inputTriangulation_->isEdgeOnBoundary(edgeId))
            return true;
        }
        break;
    }
  }
  else if(dimensionality_==3){
    switch(cellDim){
      case 0:
        return inputTriangulation_->isVertexOnBoundary(cellId);

      case 1:
        return inputTriangulation_->isEdgeOnBoundary(cellId);

      case 2:
        return inputTriangulation_->isTriangleOnBoundary(cellId);

      case 3:
        for(int i=0; i<4; ++i){
          int triangleId;
          inputTriangulation_->getCellTriangle(cellId, i, triangleId);
          if(inputTriangulation_->isTriangleOnBoundary(triangleId))
            return true;
        }
        break;
    }
  }
  return false;
}

int DiscreteGradient::getPairedCell(const Cell& cell, bool isReverse) const{
  if(dimensionality_==2){
    switch(cell.dim_){
      case 0:
        return gradient_[0][0][cell.id_];
        break;

      case 1:
        if(isReverse)
          return gradient_[0][1][cell.id_];

        return gradient_[1][1][cell.id_];
        break;

      case 2:
        if(isReverse)
          return gradient_[1][2][cell.id_];
        break;
    }
  }
  else if(dimensionality_==3){
    switch(cell.dim_){
      case 0:
        return gradient_[0][0][cell.id_];
        break;

      case 1:
        if(isReverse)
          return gradient_[0][1][cell.id_];

        return gradient_[1][1][cell.id_];
        break;

      case 2:
        if(isReverse)
          return gradient_[1][2][cell.id_];

        return gradient_[2][2][cell.id_];
        break;

      case 3:
        if(isReverse)
          return gradient_[2][3][cell.id_];
        break;
    }
  }

  return -1;
}

int DiscreteGradient::getCriticalPoints(vector<Cell>& criticalPoints) const{

  // foreach dimension
  const int numberOfDimensions=getNumberOfDimensions();
  for(int i=0; i<numberOfDimensions; ++i){

    // foreach cell of that dimension
    const int numberOfCells=getNumberOfCells(i);
    for(int j=0; j<numberOfCells; ++j){
      const Cell cell(i,j);

      if(isCellCritical(cell))
        criticalPoints.push_back(std::move(cell));
    }
  }

  return 0;
}

int DiscreteGradient::getDescendingPath(const Cell& cell, vector<Cell>& vpath) const{
  if(dimensionality_==2){
    if(cell.dim_==0){
      // assume that cellId is a vertex
      int currentId=cell.id_;
      int connectedEdgeId;
      do{
        // add a vertex
        const Cell vertex(0,currentId);
        vpath.push_back(std::move(vertex));

        if(isCellCritical(vertex)) break;

        connectedEdgeId=getPairedCell(vertex);
        if(connectedEdgeId==-1) break;

        // add an edge
        const Cell edge(1,connectedEdgeId);
        vpath.push_back(std::move(edge));

        if(isCellCritical(edge)) break;

        for(int i=0; i<2; ++i){
          int vertexId;
          inputTriangulation_->getEdgeVertex(connectedEdgeId, i, vertexId);

          if(vertexId!=currentId){
            currentId=vertexId;
            break;
          }
        }

      }while(connectedEdgeId!=-1);
    }
  }
  else if(dimensionality_==3){
    if(cell.dim_==0){
      // assume that cellId is a vertex
      int currentId=cell.id_;
      int connectedEdgeId;
      do{
        // add a vertex
        const Cell vertex(0,currentId);
        vpath.push_back(std::move(vertex));

        if(isCellCritical(vertex)) break;

        connectedEdgeId=getPairedCell(vertex);
        if(connectedEdgeId==-1) break;

        // add an edge
        const Cell edge(1,connectedEdgeId);
        vpath.push_back(std::move(edge));

        if(isCellCritical(edge)) break;

        for(int i=0; i<2; ++i){
          int vertexId;
          inputTriangulation_->getEdgeVertex(connectedEdgeId, i, vertexId);

          if(vertexId!=currentId){
            currentId=vertexId;
            break;
          }
        }

      }while(connectedEdgeId!=-1);
    }
  }

  return 0;
}

int DiscreteGradient::getDescendingPathThroughWall(const wallId_t wallId,
    const Cell& saddle2,
    const Cell& saddle1,
    const vector<wallId_t>& isVisited,
    vector<Cell>* const vpath,
    const bool enableCycleDetector) const{
  // debug
  const int numberOfEdges=inputTriangulation_->getNumberOfEdges();
  vector<char> isCycle;
  if(enableCycleDetector)
    isCycle.resize(numberOfEdges, false);

  if(dimensionality_==3){
    // add the 2-saddle to the path
    if(vpath)
      vpath->push_back(saddle2);

    int currentId=-1;
    {
      int nconnections=0;
      for(int i=0; i<3; ++i){
        int edgeId;
        inputTriangulation_->getTriangleEdge(saddle2.id_, i, edgeId);
        if(isVisited[edgeId]==wallId){
          // saddle2 can be adjacent to saddle1 on the wall
          if(isSaddle1(Cell(1,edgeId))){
            if(vpath)
              vpath->push_back(Cell(1,edgeId));
            return false;
          }

          currentId=edgeId;
          ++nconnections;
        }
      }
      if(nconnections>1) return true;
    }

    int oldId;
    do{

      // debug
      if(enableCycleDetector){
        if(!isCycle[currentId]){
          isCycle[currentId]=true;
        }
        else{
          cout << "[DiscreteGradient] Error : cycle detected on the wall of 1-saddleId=" << saddle1.id_ << endl;
          break;
        }
      }

      oldId=currentId;

      // add an edge
      const Cell edge(1,currentId);
      if(vpath)
        vpath->push_back(edge);

      if(isCellCritical(edge)) break;

      const int connectedTriangleId=getPairedCell(edge);

      // add a triangle
      const Cell triangle(2,connectedTriangleId);
      if(vpath)
        vpath->push_back(triangle);

      if(isCellCritical(triangle)) break;

      int nconnections=0;
      for(int i=0; i<3; ++i){
        int edgeId;
        inputTriangulation_->getTriangleEdge(connectedTriangleId, i, edgeId);

        if(isVisited[edgeId]==wallId and edgeId!=oldId){
          currentId=edgeId;
          ++nconnections;
        }
      }
      if(nconnections>1) return true;

      // stop at convergence caused by boundary effect
    }while(currentId!=oldId);
  }

  return false;
}

int DiscreteGradient::getAscendingPath(const Cell& cell,
    vector<Cell>& vpath,
    const bool enableCycleDetector) const{

  const int numberOfCells=inputTriangulation_->getNumberOfCells();
  vector<char> isCycle;
  if(enableCycleDetector)
    isCycle.resize(numberOfCells, false);

  if(dimensionality_==2){
    if(cell.dim_==2){
      // assume that cellId is a triangle
      int currentId=cell.id_;
      int oldId;
      do{
        oldId=currentId;

        // add a triangle
        const Cell triangle(2, currentId);
        vpath.push_back(triangle);

        if(isCellCritical(triangle)) break;

        const int connectedEdgeId=getPairedCell(triangle, true);
        if(connectedEdgeId==-1) break;

        // add an edge
        const Cell edge(1, connectedEdgeId);
        vpath.push_back(edge);

        if(isCellCritical(edge)) break;

        const int starNumber=inputTriangulation_->getEdgeStarNumber(connectedEdgeId);
        for(int i=0; i<starNumber; ++i){
          int starId;
          inputTriangulation_->getEdgeStar(connectedEdgeId, i, starId);

          if(starId!=currentId){
            currentId=starId;
            break;
          }
        }

        // stop at convergence caused by boundary effect
      }while(currentId!=oldId);
    }
  }
  else if(dimensionality_==3){
    if(cell.dim_==3){
      // assume that cellId is a tetra
      int currentId=cell.id_;
      int oldId;
      do{

        // debug
        if(enableCycleDetector){
          if(!isCycle[currentId]){
            isCycle[currentId]=true;
          }
          else{
            cout << "[DiscreteGradient] Error : cycle detected in the path from tetraId=" << cell.id_ << endl;
            break;
          }
        }

        oldId=currentId;

        // add a tetra
        const Cell tetra(3,currentId);
        vpath.push_back(tetra);

        if(isCellCritical(tetra)) break;

        const int connectedTriangleId=getPairedCell(tetra, true);
        if(connectedTriangleId==-1) break;

        // add a triangle
        const Cell triangle(2,connectedTriangleId);
        vpath.push_back(triangle);

        if(isCellCritical(triangle)) break;

        const int starNumber=inputTriangulation_->getTriangleStarNumber(connectedTriangleId);
        for(int i=0; i<starNumber; ++i){
          int starId;
          inputTriangulation_->getTriangleStar(connectedTriangleId, i, starId);

          if(starId!=currentId){
            currentId=starId;
            break;
          }
        }

        // stop at convergence caused by boundary effect
      }while(currentId!=oldId);
    }
  }

  return 0;
}

bool DiscreteGradient::getAscendingPathThroughWall(const wallId_t wallId,
    const Cell& saddle1,
    const Cell& saddle2,
    const vector<wallId_t>& isVisited,
    vector<Cell>* const vpath,
    const bool enableCycleDetector) const{
  // debug
  const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
  vector<char> isCycle;
  if(enableCycleDetector)
    isCycle.resize(numberOfTriangles, false);

  if(dimensionality_==3){
    // add the 1-saddle to the path
    if(vpath)
      vpath->push_back(saddle1);

    int currentId=-1;
    {
      int nconnections=0;
      const int triangleNumber=inputTriangulation_->getEdgeTriangleNumber(saddle1.id_);
      for(int i=0; i<triangleNumber; ++i){
        int triangleId;
        inputTriangulation_->getEdgeTriangle(saddle1.id_, i, triangleId);
        if(isVisited[triangleId]==wallId){
          // saddle1 can be adjacent to saddle2 on the wall
          if(isSaddle2(Cell(2,triangleId))){
            if(vpath)
              vpath->push_back(Cell(2,triangleId));
            return false;
          }

          currentId=triangleId;
          ++nconnections;
        }
      }
      if(nconnections>1) return true;
    }

    int oldId;
    do{

      // debug
      if(enableCycleDetector){
        if(!isCycle[currentId]){
          isCycle[currentId]=true;
        }
        else{
          cout << "[DiscreteGradient] Error : cycle detected on the wall of 2-saddleId=" << saddle2.id_ << endl;
          break;
        }
      }

      oldId=currentId;

      // add a triangle
      const Cell triangle(2,currentId);
      if(vpath)
        vpath->push_back(triangle);

      if(isCellCritical(triangle)) break;

      const int connectedEdgeId=getPairedCell(triangle, true);

      // add an edge
      const Cell edge(1,connectedEdgeId);
      if(vpath)
        vpath->push_back(edge);

      if(isCellCritical(edge)) break;

      int nconnections=0;
      const int triangleNumber=inputTriangulation_->getEdgeTriangleNumber(connectedEdgeId);
      for(int i=0; i<triangleNumber; ++i){
        int triangleId;
        inputTriangulation_->getEdgeTriangle(connectedEdgeId, i, triangleId);

        if(isVisited[triangleId]==wallId and triangleId!=oldId){
          currentId=triangleId;
          ++nconnections;
        }
      }
      if(nconnections>1) return true;

      // stop at convergence caused by boundary effect
    }while(currentId!=oldId);
  }

  return false;
}

int DiscreteGradient::getDescendingWall(const wallId_t wallId,
    const Cell& cell,
    vector<wallId_t>& isVisited,
    vector<Cell>* const wall,
    set<int>* const saddles) const{
  if(dimensionality_==3){
    if(cell.dim_==2){
      // assume that cellId is a triangle
      const int originId=cell.id_;

      queue<int> bfs;
      bfs.push(originId);

      // BFS traversal
      while(!bfs.empty()){
        const int triangleId=bfs.front();
        bfs.pop();

        if(isVisited[triangleId]!=wallId){
          isVisited[triangleId]=wallId;

          // add the triangle
          if(wall)
            wall->push_back(Cell(2,triangleId));

          for(int j=0; j<3; ++j){
            int edgeId;
            inputTriangulation_->getTriangleEdge(triangleId, j, edgeId);

            if(saddles and isSaddle1(Cell(1,edgeId)))
              saddles->insert(edgeId);

            const int pairedCellId=getPairedCell(Cell(1, edgeId));

            if(pairedCellId!=-1 and pairedCellId!=triangleId)
              bfs.push(pairedCellId);
          }
        }
      }
    }
  }

  return 0;
}

int DiscreteGradient::getAscendingWall(const wallId_t wallId,
    const Cell& cell,
    vector<wallId_t>& isVisited,
    vector<Cell>* const wall,
    set<int>* const saddles) const{
  if(dimensionality_==3){
    if(cell.dim_==1){
      // assume that cellId is an edge
      const int originId=cell.id_;

      queue<int> bfs;
      bfs.push(originId);

      // BFS traversal
      while(!bfs.empty()){
        const int edgeId=bfs.front();
        bfs.pop();

        if(isVisited[edgeId]!=wallId){
          isVisited[edgeId]=wallId;

          // add the edge
          if(wall)
            wall->push_back(Cell(1,edgeId));

          const int triangleNumber=inputTriangulation_->getEdgeTriangleNumber(edgeId);
          for(int j=0; j<triangleNumber; ++j){
            int triangleId;
            inputTriangulation_->getEdgeTriangle(edgeId, j, triangleId);

            if(saddles and isSaddle2(Cell(2,triangleId)))
              saddles->insert(triangleId);

            const int pairedCellId=getPairedCell(Cell(2, triangleId), true);

            if(pairedCellId!=-1 and pairedCellId!=edgeId)
              bfs.push(pairedCellId);
          }
        }
      }
    }
  }

  return 0;
}

int DiscreteGradient::reverseAscendingPath(const vector<Cell>& vpath){
  if(dimensionality_==2){
    // assume that the first cell is an edge
    const int numberOfCellsInPath=vpath.size();
    for(int i=0; i<numberOfCellsInPath; i+=2){
      const int edgeId=vpath[i].id_;
      const int triangleId=vpath[i+1].id_;

      gradient_[1][2][triangleId]=edgeId;
      gradient_[1][1][edgeId]=triangleId;
    }
  }
  else if(dimensionality_==3){
    // assume that the first cell is a triangle
    const int numberOfCellsInPath=vpath.size();
    for(int i=0; i<numberOfCellsInPath; i+=2){
      const int triangleId=vpath[i].id_;
      const int tetraId=vpath[i+1].id_;

      gradient_[2][3][tetraId]=triangleId;
      gradient_[2][2][triangleId]=tetraId;
    }
  }

  return 0;
}

int DiscreteGradient::reverseAscendingPathOnWall(const vector<Cell>& vpath){
  if(dimensionality_==3){
    // assume that the first cell is an edge
    const int numberOfCellsInPath=vpath.size();
    for(int i=0; i<numberOfCellsInPath; i+=2){
      const int edgeId=vpath[i].id_;
      const int triangleId=vpath[i+1].id_;

      gradient_[1][2][triangleId]=edgeId;
      gradient_[1][1][edgeId]=triangleId;
    }
  }

  return 0;
}

int DiscreteGradient::reverseDescendingPathOnWall(const vector<Cell>& vpath){
  if(dimensionality_==3){
    // assume that the first cell is a triangle
    const int numberOfCellsInPath=vpath.size();
    for(int i=0; i<numberOfCellsInPath; i+=2){
      const int triangleId=vpath[i].id_;
      const int edgeId=vpath[i+1].id_;

      gradient_[1][1][edgeId]=triangleId;
      gradient_[1][2][triangleId]=edgeId;
    }
  }

  return 0;
}

int DiscreteGradient::getEdgeIncenter(const int edgeId, float incenter[3]) const{
  int vertexId[2];
  inputTriangulation_->getEdgeVertex(edgeId, 0, vertexId[0]);
  inputTriangulation_->getEdgeVertex(edgeId, 1, vertexId[1]);

  float p[6];
  inputTriangulation_->getVertexPoint(vertexId[0],p[0],p[1],p[2]);
  inputTriangulation_->getVertexPoint(vertexId[1],p[3],p[4],p[5]);

  incenter[0]=0.5*p[0]+0.5*p[3];
  incenter[1]=0.5*p[1]+0.5*p[4];
  incenter[2]=0.5*p[2]+0.5*p[5];

  return 0;
}

int DiscreteGradient::getTriangleIncenter(const int triangleId, float incenter[3]) const{
  int vertexId[3];
  if(dimensionality_==2){
    inputTriangulation_->getCellVertex(triangleId, 0, vertexId[0]);
    inputTriangulation_->getCellVertex(triangleId, 1, vertexId[1]);
    inputTriangulation_->getCellVertex(triangleId, 2, vertexId[2]);
  }
  else if(dimensionality_==3){
    inputTriangulation_->getTriangleVertex(triangleId, 0, vertexId[0]);
    inputTriangulation_->getTriangleVertex(triangleId, 1, vertexId[1]);
    inputTriangulation_->getTriangleVertex(triangleId, 2, vertexId[2]);
  }

  float p[9];
  inputTriangulation_->getVertexPoint(vertexId[0], p[0], p[1], p[2]);
  inputTriangulation_->getVertexPoint(vertexId[1], p[3], p[4], p[5]);
  inputTriangulation_->getVertexPoint(vertexId[2], p[6], p[7], p[8]);

  float d[3];
  d[0]=Geometry::distance(p+3,p+6);
  d[1]=Geometry::distance(p,p+6);
  d[2]=Geometry::distance(p,p+3);
  const float sum=d[0]+d[1]+d[2];

  d[0]=d[0]/sum;
  d[1]=d[1]/sum;
  d[2]=d[2]/sum;

  incenter[0]=d[0]*p[0]+d[1]*p[3]+d[2]*p[6];
  incenter[1]=d[0]*p[1]+d[1]*p[4]+d[2]*p[7];
  incenter[2]=d[0]*p[2]+d[1]*p[5]+d[2]*p[8];

  return 0;
}

int DiscreteGradient::getTetraIncenter(const int tetraId, float incenter[3]) const{
  incenter[0]=0.0f;
  incenter[1]=0.0f;
  incenter[2]=0.0f;

  float p[3];
  for(int i=0; i<4; ++i){
    int triangleId;
    inputTriangulation_->getCellTriangle(tetraId, i, triangleId);

    getTriangleIncenter(triangleId, p);
    incenter[0]+=p[0];
    incenter[1]+=p[1];
    incenter[2]+=p[2];
  }

  incenter[0]/=4.0f;
  incenter[1]/=4.0f;
  incenter[2]/=4.0f;

  return 0;
}

int DiscreteGradient::getCellIncenter(const Cell& cell, float incenter[3]) const{
  switch(cell.dim_){
    case 0:
      inputTriangulation_->getVertexPoint(cell.id_, incenter[0], incenter[1], incenter[2]);
      break;

    case 1:
      getEdgeIncenter(cell.id_, incenter);
      break;

    case 2:
      getTriangleIncenter(cell.id_, incenter);
      break;

    case 3:
      getTetraIncenter(cell.id_, incenter);
      break;
  }

  return  0;
}

int DiscreteGradient::setGradientGlyphs() const{
  (*outputGradientGlyphs_numberOfPoints_)=0;
  (*outputGradientGlyphs_numberOfCells_)=0;

  int pointId{};
  int cellId{};

  // foreach dimension
  const int numberOfDimensions=getNumberOfDimensions();
  for(int i=0; i<numberOfDimensions-1; ++i){
    // foreach cell of that dimension
    const int numberOfCells=getNumberOfCells(i);
    for(int j=0; j<numberOfCells; ++j){
      const Cell cell(i,j);

      const int pairedCellId=getPairedCell(cell);
      if(pairedCellId!=-1){
        // get gradient pair
        const int pairedCellDim=i+1;

        const Cell pairedCell(pairedCellDim,pairedCellId);

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
        outputGradientGlyphs_cells_->push_back(pointId+1);

        outputGradientGlyphs_cells_pairTypes_->push_back(i);

        pointId+=2;
        cellId+=1;
      }
    }
  }

  (*outputGradientGlyphs_numberOfPoints_)=pointId;
  (*outputGradientGlyphs_numberOfCells_)=cellId;

  return 0;
}
