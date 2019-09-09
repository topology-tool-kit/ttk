#include <MorseSmaleComplex2D.h>

using namespace std;
using namespace ttk;
using namespace dcg;

MorseSmaleComplex2D::MorseSmaleComplex2D() : AbstractMorseSmaleComplex() {
}

MorseSmaleComplex2D::~MorseSmaleComplex2D() {
}

int MorseSmaleComplex2D::getAscendingSeparatrices1(
  const vector<Cell> &criticalPoints,
  vector<Separatrix> &separatrices,
  vector<vector<Cell>> &separatricesGeometry) const {

  vector<SimplexId> saddleIndexes;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 1)
      saddleIndexes.push_back(i);
  }
  const SimplexId numberOfSaddles = saddleIndexes.size();

  // estimation of the number of separatrices, apriori :
  // numberOfAscendingPaths=2, numberOfDescendingPaths=2
  const SimplexId numberOfSeparatrices = 4 * numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const SimplexId saddleIndex = saddleIndexes[i];
    const Cell &saddle = criticalPoints[saddleIndex];

    // add ascending vpaths
    {
      const SimplexId starNumber
        = inputTriangulation_->getEdgeStarNumber(saddle.id_);
      for(SimplexId j = 0; j < starNumber; ++j) {
        const SimplexId shift = j;

        SimplexId triangleId;
        inputTriangulation_->getEdgeStar(saddle.id_, j, triangleId);

        vector<Cell> vpath;
        vpath.push_back(saddle);
        discreteGradient_.getAscendingPath(Cell(2, triangleId), vpath);

        const Cell &lastCell = vpath.back();
        if(lastCell.dim_ == 2 and discreteGradient_.isCellCritical(lastCell)) {
          const SimplexId separatrixIndex = 4 * i + shift;

          separatricesGeometry[separatrixIndex] = std::move(vpath);
          separatrices[separatrixIndex]
            = Separatrix(true, saddle, lastCell, false, separatrixIndex);
        }
      }
    }
  }

  return 0;
}
