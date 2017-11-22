#ifndef _BOTTLENECKDISTANCEIMPL_H
#define _BOTTLENECKDISTANCEIMPL_H

using namespace std;

template <typename dataType>
int BottleneckDistance::execute(
  const bool usePersistenceMetric,
  const double alpha)
{
  Timer t;

  this->computeBottleneck(
    static_cast<const vector<diagramTuple>*> (outputCT1_),
    static_cast<const vector<diagramTuple>*> (outputCT2_),
    static_cast<vector<matchingTuple>*> (matchings_),
    usePersistenceMetric,
    alpha);

  {
    stringstream msg;
    msg << "[BottleneckDistance] Data-set processed in "
        << t.getElapsedTime() << " s. (" << threadNumber_
        << " thread(s))."
        << endl;
    msg << "[BottleneckDistance] distance = " << *(dataType*) distance_ << endl;
    // Using std:: with cout for CLion parser.
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
double BottleneckDistance::computeGeometricalRange(
  const vector<diagramTuple> *CTDiagram1,
  const vector<diagramTuple> *CTDiagram2,
  const int d1Size,
  const int d2Size) const
{
  float minX1, maxX1, minY1, maxY1, minZ1, maxZ1;
  float minX2, maxX2, minY2, maxY2, minZ2, maxZ2;
  float minX, minY, minZ, maxX, maxY, maxZ;
  minX1 = minY1 = minZ1 = minX2 = minY2 = minZ2 = numeric_limits<float>::max();
  maxX1 = maxY1 = maxZ1 = maxX2 = maxY2 = maxZ2 = numeric_limits<float>::min();

  for (int i = 0; i < d1Size; ++i) {
    diagramTuple t = CTDiagram1->at(i);
    float xa = get<7>(t), ya = get<8>(t), za = get<9>(t);
    float xb = get<11>(t), yb = get<12>(t), zb = get<13>(t);
    minX1 = min(min(minX1, xa), xb);
    minY1 = min(min(minY1, ya), yb);
    minZ1 = min(min(minZ1, za), zb);
    maxX1 = max(max(maxX1, xa), xb);
    maxY1 = max(max(maxY1, ya), yb);
    maxZ1 = max(max(maxZ1, za), zb);
  }

  for (int i = 0; i < d2Size; ++i) {
    diagramTuple t = CTDiagram2->at(i);
    float xa = get<7>(t), ya = get<8>(t), za = get<9>(t);
    float xb = get<11>(t), yb = get<12>(t), zb = get<13>(t);
    minX2 = min(min(minX2, xa), xb);
    minY2 = min(min(minY2, ya), yb);
    minZ2 = min(min(minZ2, za), zb);
    maxX2 = max(max(maxX2, xa), xb);
    maxY2 = max(max(maxY2, ya), yb);
    maxZ2 = max(max(maxZ2, za), zb);
  }

  minX = min(minX1, minX2); maxX = max(maxX1, maxX2);
  minY = min(minY1, minY2); maxY = max(maxY1, maxY2);
  minZ = min(minZ1, minZ2); maxZ = max(maxZ1, maxZ2);

  return sqrt(
      pow(maxX - minX, 2) + pow(maxY - minY, 2) + pow(maxZ - minZ, 2));
}

template <typename dataType>
double BottleneckDistance::computeMinimumRelevantPersistence(
  const vector<diagramTuple> *CTDiagram1,
  const vector<diagramTuple> *CTDiagram2,
  const int d1Size,
  const int d2Size) const
{
  vector<dataType> toSort;
  for (int i = 0; i < d1Size; ++i) {
    diagramTuple t = CTDiagram1->at(i);
    dataType persistence = abs<dataType>(get<4>(t));
    toSort.push_back(persistence);
  }
  for (int i = 0; i < d2Size; ++i) {
    diagramTuple t = CTDiagram2->at(i);
    dataType persistence = abs<dataType>(get<4>(t));
    toSort.push_back(persistence);
  }
  sort(toSort.begin(), toSort.end());
  double epsilon = 0.0000001;
  int largeSize = 2000;
  dataType zeroThresh = (dataType) epsilon;
  if (d1Size + d2Size > largeSize + 1) {
    zeroThresh = toSort.at(d1Size + d2Size - largeSize);
    if (toSort.at(d1Size + d2Size - (largeSize+1)) == zeroThresh)
      zeroThresh += (dataType) epsilon;
  }
  if (zeroThresh < epsilon) zeroThresh = epsilon;
  return zeroThresh;
}

template <typename dataType>
void BottleneckDistance::computeMinMaxSaddleNumberAndMapping(
  const vector<diagramTuple > *CTDiagram,
  int dSize,
  int &nbMin,
  int &nbMax,
  int &nbSaddle,
  vector<int> *minMap,
  vector<int> *maxMap,
  vector<int> *sadMap,
  const dataType zeroThresh)
{
  for (int i = 0; i < dSize; ++i) {
    diagramTuple t = CTDiagram->at(i);
    BNodeType nt1 = get<1>(t);
    BNodeType nt2 = get<3>(t);
    dataType dt = get<4>(t);
    if (abs<dataType>(dt) < zeroThresh) continue;

    if (nt1 == BLocalMin && nt2 == BLocalMax) {
      ++nbMax;
      maxMap->push_back(i);
    }
    else {
      if (nt1 == BLocalMax || nt2 == BLocalMax) {
        ++nbMax;
        maxMap->push_back(i);
      }
      if (nt1 == BLocalMin || nt2 == BLocalMin) {
        ++nbMin;
        minMap->push_back(i);
      }
      if ((nt1 == BSaddle1 && nt2 == BSaddle2)
          || (nt1 == BSaddle2 && nt2 == BSaddle1)) {
        ++nbSaddle;
        sadMap->push_back(i);
      }
    }
  }
}

template <typename dataType>
void BottleneckDistance::buildCostMatrices(
  const vector<diagramTuple> *CTDiagram1,
  const vector<diagramTuple> *CTDiagram2,
  const int d1Size,
  const int d2Size,
  function<dataType (const diagramTuple, const diagramTuple)>& distanceFunction,
  function<dataType (const diagramTuple)>& diagonalDistanceFunction,
  const double zeroThresh,
  dataType **minMatrix,
  dataType **maxMatrix,
  dataType **sadMatrix)
{

  int maxI = 0, minI = 0;
  int maxJ = 0, minJ = 0;
  int sadI = 0, sadJ = 0;

  for (int i = 0; i < d1Size; ++i)
  {
    diagramTuple t1 = CTDiagram1->at(i);
    if (abs<dataType>(get<4>(t1)) < zeroThresh) continue;

    BNodeType t11 = get<1>(t1);
    BNodeType t13 = get<3>(t1);
    bool isMin1 = (t11 == BLocalMin || t13 == BLocalMin);
    bool isMax1 = (t11 == BLocalMax || t13 == BLocalMax);
    bool isSad1 = (t11 == BSaddle1 && t13 == BSaddle2) ||
      (t11 == BSaddle2 && t13 == BSaddle1);
    if (t11 == BLocalMin && t13 == BLocalMax) {
      isMin1 = false;
      isMax1 = true;
    }

    minJ = 0;
    maxJ = 0;
    sadJ = 0;

    for (int j = 0; j < d2Size; ++j)
    {
      diagramTuple t2 = CTDiagram2->at(j);
      if (abs<dataType>(get<4>(t2)) < zeroThresh) continue;

      BNodeType t21 = get<1>(t2);
      BNodeType t23 = get<3>(t2);
      bool isMin2 = (t21 == BLocalMin || t23 == BLocalMin);
      bool isMax2 = (t21 == BLocalMax || t23 == BLocalMax);
      bool isSad2 = (t21 == BSaddle1 && t23 == BSaddle2) ||
        (t21 == BSaddle2 && t23 == BSaddle1);
      if (t21 == BLocalMin && t23 == BLocalMax) {
        isMin2 = false;
        isMax2 = true;
      }
      if ((isMin1 && !isMin2) || (isMax1 && !isMax2) || (isSad1 && !isSad2))
        continue;

      dataType distance = distanceFunction(t1, t2);

      if (isMin1 && isMin2)
        minMatrix[minI][minJ++] = distance;
      else if (isMax1 && isMax2)
        maxMatrix[maxI][maxJ++] = distance;
      else if (isSad1 && isSad2)
        sadMatrix[sadI][sadJ++] = distance;
    }

    dataType distanceToDiagonal = diagonalDistanceFunction(t1);
    if (isMin1) minMatrix[minI][minJ++] = distanceToDiagonal;
    if (isMax1) maxMatrix[maxI][maxJ++] = distanceToDiagonal;
    if (isSad1) sadMatrix[sadI][sadJ++] = distanceToDiagonal;

    if (isMin1) ++minI;
    if (isMax1) ++maxI;
    if (isSad1) ++sadI;
  }

  minJ = 0;
  maxJ = 0;
  sadJ = 0;

  // Last row: match remaining J components with diagonal.
  for (int j = 0; j < d2Size; ++j) {
    diagramTuple t2 = CTDiagram2->at(j);
    if (abs<dataType>(get<4>(t2)) < zeroThresh) continue;

    BNodeType t21 = get<1>(t2);
    BNodeType t23 = get<3>(t2);
    bool isMin2 = (t21 == BLocalMin || t23 == BLocalMin);
    bool isMax2 = (t21 == BLocalMax || t23 == BLocalMax);
    bool isSad2 = (t21 == BSaddle1 && t23 == BSaddle2) ||
                  (t21 == BSaddle2 && t23 == BSaddle1);
    if (t21 == BLocalMin && t23 == BLocalMax) {
      isMin2 = false;
      isMax2 = true;
    }

    dataType distanceToDiagonal = diagonalDistanceFunction(t2);
    if (isMin2) minMatrix[minI][minJ++] = distanceToDiagonal;
    if (isMax2) maxMatrix[maxI][maxJ++] = distanceToDiagonal;
    if (isSad2) sadMatrix[sadI][sadJ++] = distanceToDiagonal;
  }

  // Last cell
  minMatrix[minI][minJ] = numeric_limits<dataType>::max();
  maxMatrix[maxI][maxJ] = numeric_limits<dataType>::max();
  sadMatrix[sadI][sadJ] = numeric_limits<dataType>::max();
}

template <typename dataType>
dataType BottleneckDistance::findInitialThreshold(
  const int nbRow,
  const int nbCol,
  const dataType **matrix) const
{
  dataType col, row;
  vector<dataType> rowCols(0);
  for (int i = 0; i < nbRow; ++i) {
    row = matrix[i][0];
    for (int j = 1; j < nbCol; ++j) {
      row = min(row, matrix[i][j]);
    }
    rowCols.push_back(row);
  }

  for (int j = 0; j < nbCol; ++j) {
    col = matrix[0][j];
    for (int i = 0; i < nbRow; ++i) {
      col = min(col, matrix[i][j]);
    }
    rowCols.push_back(col);
  }

  dataType threshold = rowCols.size() > 0 ? rowCols.at(0) : 0;
  for (int i = 1, s = (int) rowCols.size(); i < s; ++i)
    threshold = max(threshold, rowCols.at(i));

  return threshold;
}

template <typename dataType>
void BottleneckDistance::filterFromThreshold(
  const dataType threshold,
  const int nbRow,
  const int nbCol,
  const dataType **matrix,
  dataType **bottleneckMatrix)
{
  for (int i = 0; i < nbRow; ++i)
    for (int j = 0; j < nbCol; ++j)
      bottleneckMatrix[i][j] =
        matrix[i][j] > threshold ? numeric_limits<dataType>::max() : matrix[i][j];
}

template <typename dataType>
bool BottleneckDistance::isValidMatching(
  const vector<matchingTuple>* matchings,
  const dataType thresholdMin) const
{
  for (int i = 0, s = matchings->size(); i < s; ++i) {
    matchingTuple t = matchings->at(i);
    // Safe enough (matrix minima are substracted as a first
    // Munkres step, but numeric_limits should be high enough)
    if (get<2>(t) > thresholdMin)
      return false;
  }
  return true;
}

template <typename dataType>
void BottleneckDistance::iterateSolving(
  vector<matchingTuple>* matchings,
  const dataType threshold,
  const int nbRow,
  const int nbCol,
  const dataType **matrix,
  dataType **bottleneckMatrix,
  Munkres *solver)
{
  // Progressively solve easier affectation problems until
  // a non-infinite cost is found.

  dataType adaptiveThreshold = threshold;

  // Check assignment maximum (infty => restart Munkres with a higher - next - threshold).
  int minIter = 2;
  while (!this->isValidMatching(const_cast<const vector<matchingTuple>* >(matchings), adaptiveThreshold)
         && minIter <= nbRow * nbCol)
  {
    matchings->clear();

    // Find next threshold.
    dataType newThreshold = numeric_limits<dataType>::max();
    for (int i = 0; i < nbRow; ++i) {
      for (int j = 0; j < nbCol; ++j) {
        dataType e = matrix[i][j];
        if (e > adaptiveThreshold && e < newThreshold) newThreshold = e;
      }
    }
    adaptiveThreshold = newThreshold;

    // Apply Munkres again.
    for (int i = 0; i < nbRow; ++i) {
      for (int j = 0; j < nbCol; ++j) {
        dataType e = matrix[i][j];
        bottleneckMatrix[i][j] =
          (e > adaptiveThreshold) ? numeric_limits<dataType>::max() : e;
      }
    }

    solver->run<dataType>(matchings);

    minIter++;
  }
}

template <typename dataType>
void BottleneckDistance::solvePWasserstein(
  const int nbRow,
  const int nbCol,
  dataType **matrix,
  vector<matchingTuple> *matchings,
  Munkres *solver)
{
  solver->setInput(nbRow+1, nbCol+1, (void*) matrix);
  solver->run<dataType>(matchings);
  solver->clearMatrix<dataType>();
}

template <typename dataType>
void BottleneckDistance::solveInfinityWasserstein(
  const int nbRow,
  const int nbCol,
  dataType **matrix,
  vector<matchingTuple> *matchings,
  Munkres *solver)
{
  // Cf. "A New Algorithm for Solving Linear Bottleneck Assignment Problem"
  // P. S. Pundir, S. K. Porwal, B. P. Singh
  // Journal of Institute of Science and Technology.

  // Filter values.
  dataType minThreshold =
    this->findInitialThreshold(
      nbRow+1, nbCol+1,
      const_cast<const dataType**>(matrix));

  auto **bottleneckMatrix = new dataType*[nbRow+1];
  for (int i = 0; i < nbRow+1; ++i)
    bottleneckMatrix[i] = new dataType[nbCol+1];

  this->filterFromThreshold(
    minThreshold, nbRow+1, nbCol+1,
    const_cast<const dataType**>(matrix), bottleneckMatrix);

  // Solve for the first time.
  solver->setInput(nbRow+1, nbCol+1, (void*) bottleneckMatrix);
  solver->run<dataType>(matchings);

  // Loop solve.
  this->iterateSolving(
    matchings, minThreshold, nbRow+1, nbCol+1,
    const_cast<const dataType**>(matrix), bottleneckMatrix, solver);

  solver->clearMatrix<dataType>();
}

template <typename dataType>
dataType BottleneckDistance::buildMappings(
  const vector<matchingTuple> inputMatchings,
  vector<matchingTuple> *outputMatchings,
  const vector<int> map1,
  const vector<int> map2,
  int wasserstein)
{
  dataType addedPersistence = 0;
  for (int i = 0, s = (int) inputMatchings.size(); i < s; ++i) {
    matchingTuple t = inputMatchings.at(i);
    dataType val = abs<dataType>(get<2>(t));

    if (get<0>(t) >= (int) map1.size()) {
      addedPersistence =
        (wasserstein > 0 ?
        addedPersistence + pow(val, wasserstein) :
        max(val, addedPersistence));
    }
    else if (get<1>(t) >= (int) map2.size()) {
      addedPersistence =
        (wasserstein > 0 ?
        addedPersistence + pow(val, wasserstein) :
        max(val, addedPersistence));
    } else {
      // Using std:: and casting for CLion parser.
      outputMatchings->push_back(
        std::make_tuple(map1.at(get<0>(t)), map2.at(get<1>(t)), val)
      );
    }
  }

  return addedPersistence;
}

#endif
