#include<MorseSmaleComplex2D.h>

MorseSmaleComplex2D::MorseSmaleComplex2D():
  AbstractMorseSmaleComplex()
{
}

MorseSmaleComplex2D::~MorseSmaleComplex2D(){
}

int MorseSmaleComplex2D::getAscendingSeparatrices1(const vector<Cell>& criticalPoints,
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
    {
      const int starNumber=inputTriangulation_->getEdgeStarNumber(saddle.id_);
      for(int j=0; j<starNumber; ++j){
        const int shift=j;

        int triangleId;
        inputTriangulation_->getEdgeStar(saddle.id_, j, triangleId);

        vector<Cell> vpath;
        vpath.push_back(saddle);
        discreteGradient_.getAscendingPath(Cell(2,triangleId), vpath);

        const Cell& lastCell=vpath.back();
        if(lastCell.dim_==2 and discreteGradient_.isCellCritical(lastCell)){
          const int separatrixIndex=4*i+shift;

          separatricesGeometry[separatrixIndex]=std::move(vpath);
          separatrices[separatrixIndex]=std::move(Separatrix(true,saddle,lastCell,false,separatrixIndex));
        }
      }
    }
  }

  return 0;
}
