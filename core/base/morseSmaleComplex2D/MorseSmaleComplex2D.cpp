#include<MorseSmaleComplex2D.h>

MorseSmaleComplex2D::MorseSmaleComplex2D():
  AbstractMorseSmaleComplex()
{
}

MorseSmaleComplex2D::~MorseSmaleComplex2D(){
}

int MorseSmaleComplex2D::getSeparatrices(const vector<Cell>& criticalPoints,
    vector<Separatrix>& separatrices,
    vector<vector<Cell>>& separatricesGeometry) const{

  vector<int> saddleIndexes;
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==1)
      saddleIndexes.push_back(i);
  }
  const int numberOfSaddles=saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfAscendingPaths=2, numberOfDescendingPaths=2
  const int numberOfSeparatrices=4*numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSaddles; ++i){
    const int saddleIndex=saddleIndexes[i];
    const Cell& saddle=criticalPoints[saddleIndex];

    // add ascending vpaths
    if(ComputeAscendingSeparatrices1){
      const int starNumber=inputTriangulation_->getEdgeStarNumber(saddle.id_);
      for(int j=0; j<starNumber; ++j){
        const int shift=j;

        int triangleId;
        inputTriangulation_->getEdgeStar(saddle.id_, j, triangleId);

        vector<Cell> vpath{saddle};
        discreteGradient_.getAscendingPath(Cell(2,triangleId), vpath);

        const Cell& lastCell=vpath.back();
        if(lastCell.dim_==2 and discreteGradient_.isCellCritical(lastCell)){
          const int separatrixIndex=4*i+shift;

          separatricesGeometry[separatrixIndex]=std::move(vpath);
          separatrices[separatrixIndex]=std::move(Separatrix(true,saddle,lastCell,false,separatrixIndex));
        }
      }
    }

    // add descending vpaths
    if(ComputeDescendingSeparatrices1){
      for(int j=0; j<2; ++j){
        const int shift=j+2;

        int vertexId;
        inputTriangulation_->getEdgeVertex(saddle.id_, j, vertexId);

        vector<Cell> vpath{saddle};
        discreteGradient_.getDescendingPath(Cell(0,vertexId), vpath);

        const Cell& lastCell=vpath.back();
        if(lastCell.dim_==0 and discreteGradient_.isCellCritical(lastCell)){
          const int separatrixIndex=4*i+shift;

          separatricesGeometry[separatrixIndex]=std::move(vpath);
          separatrices[separatrixIndex]=std::move(Separatrix(true,saddle,lastCell,false,separatrixIndex));
        }
      }
    }
  }

  return 0;
}

int MorseSmaleComplex2D::setAscendingSegmentation(const vector<Cell>& criticalPoints,
    vector<int>& maxSeeds,
    int* const morseSmaleManifold,
    int& numberOfMaxima) const{
  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
  std::fill(morseSmaleManifold,morseSmaleManifold+numberOfVertices, -1);

  const int numberOfCells=inputTriangulation_->getNumberOfCells();
  vector<int> morseSmaleManifoldOnCells(numberOfCells, -1);

  // get the seeds : maxima
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==2)
      maxSeeds.push_back(criticalPoint.id_);
  }
  const int numberOfSeeds=maxSeeds.size();
  numberOfMaxima=numberOfSeeds;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSeeds; ++i){
    queue<int> bfs;

    // push the seed
    {
      const int seedId=maxSeeds[i];
      bfs.push(seedId);
    }

    // BFS traversal
    while(!bfs.empty()){
      const int triangleId=bfs.front();
      bfs.pop();

      if(morseSmaleManifoldOnCells[triangleId]==-1){
        morseSmaleManifoldOnCells[triangleId]=i;

        for(int j=0; j<3; ++j){
          int edgeId;
          inputTriangulation_->getCellEdge(triangleId, j, edgeId);

          const int starNumber=inputTriangulation_->getEdgeStarNumber(edgeId);
          for(int k=0; k<starNumber; ++k){
            int neighborId;
            inputTriangulation_->getEdgeStar(edgeId, k, neighborId);

            const int pairedCellId=discreteGradient_.getPairedCell(Cell(2, neighborId), true);

            if(pairedCellId==edgeId)
              bfs.push(neighborId);
          }
        }
      }
    }
  }

  // put segmentation infos from cells to points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfVertices; ++i){
    int starId;
    inputTriangulation_->getVertexStar(i, 0, starId);
    morseSmaleManifold[i]=morseSmaleManifoldOnCells[starId];
  }

  return 0;
}

int MorseSmaleComplex2D::setDescendingSegmentation(const vector<Cell>& criticalPoints,
    int* const morseSmaleManifold,
    int& numberOfMinima) const{
  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
  std::fill(morseSmaleManifold,morseSmaleManifold+numberOfVertices, -1);

  // get the seeds : minima
  vector<int> seeds;
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==0)
      seeds.push_back(criticalPoint.id_);
  }
  const int numberOfSeeds=seeds.size();
  numberOfMinima=numberOfSeeds;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSeeds; ++i){
    queue<int> bfs;

    // push the seed
    {
      const int seedId=seeds[i];
      bfs.push(seedId);
    }

    // BFS traversal
    while(!bfs.empty()){
      const int vertexId=bfs.front();
      bfs.pop();

      if(morseSmaleManifold[vertexId]==-1){
        morseSmaleManifold[vertexId]=i;

        const int edgeNumber=inputTriangulation_->getVertexEdgeNumber(vertexId);
        for(int j=0; j<edgeNumber; ++j){
          int edgeId;
          inputTriangulation_->getVertexEdge(vertexId, j, edgeId);

          for(int k=0; k<2; ++k){
            int neighborId;
            inputTriangulation_->getEdgeVertex(edgeId, k, neighborId);

            const int pairedCellId=discreteGradient_.getPairedCell(Cell(0, neighborId));

            if(pairedCellId==edgeId)
              bfs.push(neighborId);
          }
        }
      }
    }
  }

  return 0;
}

int MorseSmaleComplex2D::setFinalSegmentation(const int numberOfMaxima,
    const int numberOfMinima,
    const int* const ascendingManifold,
    const int* const descendingManifold,
    int* const morseSmaleManifold) const{
  vector<vector<pair<int,int>>> minTable(numberOfMinima);

  int id{};
  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
  for(int i=0; i<numberOfVertices; ++i){
    const int d=ascendingManifold[i];
    const int a=descendingManifold[i];

    if(a==-1 or d==-1){
      morseSmaleManifold[i]=-1;
      continue;
    }

    vector<pair<int,int>>& table=minTable[a];
    int foundId=-1;
    for(const pair<int,int>& p : table){
      if(p.first == d)
        foundId=p.second;
    }

    // add new association (a,d)
    if(foundId==-1){
      table.push_back(make_pair(d,id));
      morseSmaleManifold[i]=id;
      ++id;
    }
    // update to saved associationId
    else
      morseSmaleManifold[i]=foundId;

  }

  return 0;
}
