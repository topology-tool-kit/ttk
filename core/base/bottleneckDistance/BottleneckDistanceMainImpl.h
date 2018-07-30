#ifndef _BOTTLENECKDISTANCEIMPLMAIN_H
#define _BOTTLENECKDISTANCEIMPLMAIN_H

//  vector <   -- diagram
//    tuple <    -- pair of critical points
//      int , NodeType
//      int, NodeType
//      dataType  -- persistance of pair
//      int  -- type (0/min, 1/saddle, 2/max)
//      dataType            -- scalar value at vertex 1
//      float, float, float -- vertex 1 coordinates
//      dataType            -- scalar value at vertex 2
//      float, float, float -- vertex 2 coordinates
template <typename dataType>
int BottleneckDistance::computeBottleneck(
    const std::vector<diagramTuple> *d1,
    const std::vector<diagramTuple> *d2,
    std::vector<matchingTuple> *matchings,
    const bool usePersistenceMetric,
    const double alpha)
{
  int d1Size = (int) d1->size();
  int d2Size = (int) d2->size();

  bool transposeOriginal = d1Size > d2Size;
  const std::vector<diagramTuple> *CTDiagram1 = transposeOriginal ? d2 : d1;
  const std::vector<diagramTuple> *CTDiagram2 = transposeOriginal ? d1 : d2;
  if (transposeOriginal) {
    int temp = d1Size;
    d1Size = d2Size;
    d2Size = temp;
  }

  if (transposeOriginal) {
    std::stringstream msg;
    msg << "[BottleneckDistance] The first persistence diagram is larger than the second." << std::endl;
    msg << "[BottleneckDistance] Solving the transposed problem." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  // Check user parameters.
  const double geometricalFactor = alpha < 0.0 ? 0.0 : alpha > 1.0 ? 1.0 : alpha;
  const int wasserstein = (wasserstein_ == "inf") ? -1 : stoi(wasserstein_);
  if (wasserstein < 0 && wasserstein != -1) return -4;

  // Needed for distance computation.
  const double maxDistance =
      this->computeGeometricalRange<dataType>(CTDiagram1, CTDiagram2, d1Size, d2Size);

  // Needed to limit computation time.
  const dataType zeroThresh =
      this->computeMinimumRelevantPersistence<dataType>(CTDiagram1, CTDiagram2, d1Size, d2Size);

  // Initialize solvers.
  std::vector<matchingTuple> minMatchings;
  std::vector<matchingTuple> maxMatchings;
  std::vector<matchingTuple> sadMatchings;

  // Initialize cost matrices.
  int nbRowMin = 0, nbColMin = 0;
  int maxRowColMin = 0, minRowColMin = 0;
  int nbRowMax = 0, nbColMax = 0;
  int maxRowColMax = 0, minRowColMax = 0;
  int nbRowSad = 0, nbColSad = 0;
  int maxRowColSad = 0, minRowColSad = 0;

  // Remap for matchings.
  std::vector<int> minMap1; std::vector<int> minMap2;
  std::vector<int> maxMap1; std::vector<int> maxMap2;
  std::vector<int> sadMap1; std::vector<int> sadMap2;

  this->computeMinMaxSaddleNumberAndMapping(
      CTDiagram1, d1Size, nbRowMin, nbRowMax, nbRowSad,
      &minMap1, &maxMap1, &sadMap1, zeroThresh);
  this->computeMinMaxSaddleNumberAndMapping(
      CTDiagram2, d2Size, nbColMin, nbColMax, nbColSad,
      &minMap2, &maxMap2, &sadMap2, zeroThresh);

  // Automatically transpose if nb rows > nb cols
  maxRowColMin = std::max(nbRowMin+1, nbColMin+1);
  maxRowColMax = std::max(nbRowMax+1, nbColMax+1);
  maxRowColSad = std::max(nbRowSad+1, nbColSad+1);

  minRowColMin = std::min(nbRowMin+1, nbColMin+1);
  minRowColMax = std::min(nbRowMax+1, nbColMax+1);
  minRowColSad = std::min(nbRowSad+1, nbColSad+1);

  std::vector<std::vector<dataType>> minMatrix((unsigned long) minRowColMin, std::vector<dataType>(maxRowColMin));
  std::vector<std::vector<dataType>> maxMatrix((unsigned long) minRowColMax, std::vector<dataType>(maxRowColMax));
  std::vector<std::vector<dataType>> sadMatrix((unsigned long) minRowColSad, std::vector<dataType>(maxRowColSad));

  double px = px_;
  double py = py_;
  double pz = pz_;
  double pe = pe_;
  double ps = ps_;

  std::function<dataType (const diagramTuple, const diagramTuple)>
      distanceFunction =
      [maxDistance, wasserstein, geometricalFactor, px, py, pz, pe, ps]
          (const diagramTuple a, const diagramTuple b) -> dataType {
        dataType rX = std::get<6>(a);
        dataType rY = std::get<10>(a);
        dataType cX = std::get<6>(b);
        dataType cY = std::get<10>(b);
        dataType x = ps * abs_diff<dataType>(rX, cX);
        dataType y = pe * abs_diff<dataType>(rY, cY);
        double dist = sqrt(
            px * pow((/*std::get<7>(a)+*/std::get<11>(a))/*/2*/ - (/*std::get<7>(b)+*/std::get<11>(b))/*/2*/, 2) +
            py * pow((/*std::get<8>(a)+*/std::get<12>(a))/*/2*/ - (/*std::get<8>(b)+*/std::get<12>(b))/*/2*/, 2) +
            pz * pow((/*std::get<9>(a)+*/std::get<13>(a))/*/2*/ - (/*std::get<9>(b)+*/std::get<13>(b))/*/2*/, 2)
        );
        //dist /= maxDistance;
        double infDistance = std::max(x, y);
        double val = geometricalFactor * infDistance + (1.0 - geometricalFactor) * dist;
        val = wasserstein > 1 ? pow(val, wasserstein) : val;
        return val;
      };

  std::function<dataType (const diagramTuple)>
      diagonalDistanceFunction =
      [wasserstein, geometricalFactor, px, py, pz, ps, pe]
          (const diagramTuple a) -> dataType {
        dataType rX = std::get<6>(a);
        dataType rY = std::get<10>(a);
        double x1 = std::get<7>(a);
        double y1 = std::get<8>(a);
        double z1 = std::get<9>(a);
        double x2 = std::get<11>(a);
        double y2 = std::get<12>(a);
        double z2 = std::get<13>(a);
        double infDistance = ps * abs_diff<dataType>(rX, rY);
        double dist =
            sqrt(px * pow(x2 - x1, 2) + py * pow(y2 - y1, 2) + pz * pow(z2 - z1, 2));
        double val = geometricalFactor * infDistance + (1.0 - geometricalFactor) * dist;
        return wasserstein > 1 ? pow(val, wasserstein) : val;
      };

  const bool transposeMin = nbRowMin > nbColMin;
  const bool transposeMax = nbRowMax > nbColMax;
  const bool transposeSad = nbRowSad > nbColSad;

  Timer t;

  this->buildCostMatrices(
      CTDiagram1, CTDiagram2, d1Size, d2Size,
      distanceFunction, diagonalDistanceFunction,
      zeroThresh, minMatrix, maxMatrix, sadMatrix,
      transposeMin,
      transposeMax,
      transposeSad,
      wasserstein
  );

  if (wasserstein > 0) {

    if (nbRowMin > 0 && nbColMin > 0) {
      Munkres solverMin;
      dMsg(std::cout, "[BottleneckDistance] Affecting minima...\n", timeMsg);
      this->solvePWasserstein(minRowColMin, maxRowColMin, minMatrix, &minMatchings, &solverMin);
    }

    if (nbRowMax > 0 && nbColMax > 0) {
      Munkres solverMax;
      dMsg(std::cout, "[BottleneckDistance] Affecting maxima...\n", timeMsg);
      this->solvePWasserstein(minRowColMax, maxRowColMax, maxMatrix, &maxMatchings, &solverMax);
    }

    if (nbRowSad > 0 && nbColSad > 0) {
      Munkres solverSad;
      dMsg(std::cout, "[BottleneckDistance] Affecting saddles...\n", timeMsg);
      this->solvePWasserstein(minRowColSad, maxRowColSad, sadMatrix, &sadMatchings, &solverSad);
    }

  } else {

    // Launch solving for minima.
    if (nbRowMin > 0 && nbColMin > 0) {
      GabowTarjan solverMin;
      dMsg(std::cout, "[BottleneckDistance] Affecting minima...\n", timeMsg);
      this->solveInfinityWasserstein(
        minRowColMin, maxRowColMin, nbRowMin, nbColMin,
        minMatrix, &minMatchings, &solverMin);
    }

    // Launch solving for maxima.
    if (nbRowMax > 0 && nbColMax > 0) {
      GabowTarjan solverMax;
      dMsg(std::cout, "[BottleneckDistance] Affecting maxima...\n", timeMsg);
      this->solveInfinityWasserstein(
        minRowColMax, maxRowColMax, nbRowMax, nbColMax,
        maxMatrix, &maxMatchings, &solverMax);
    }

    // Launch solving for saddles.
    if (nbRowSad > 0 && nbColSad > 0) {
      GabowTarjan solverSad;
      dMsg(std::cout, "[BottleneckDistance] Affecting saddles...\n", timeMsg);
      this->solveInfinityWasserstein(
        minRowColSad, maxRowColSad, nbRowSad, nbColSad,
        sadMatrix, &sadMatchings, &solverSad);
    }
  }

  std::cout << "[BottleneckDistance] TTK CORE DONE IN " << t.getElapsedTime() << std::endl;

  // Rebuild mappings.
  // Begin cost computation for unpaired vertices.
  // std::cout << "Min" << std::endl;
  dataType addedMinPersistence = this->buildMappings<dataType>(
    minMatchings, transposeOriginal, transposeMin,
    matchings, minMap1, minMap2, wasserstein);

  // std::cout << "Max" << std::endl;
  dataType addedMaxPersistence = this->buildMappings<dataType>(
    maxMatchings, transposeOriginal, transposeMax,
    matchings, maxMap1, maxMap2, wasserstein);

  // std::cout << "Sad" << std::endl;
  dataType addedSadPersistence = this->buildMappings<dataType>(
    sadMatchings, transposeOriginal, transposeSad,
    matchings, sadMap1, sadMap2, wasserstein);

  // TODO [HIGH] do that for embeddings
  // Recompute matching weights for user-friendly distance.
  dataType d = 0;
  std::vector<bool> paired1(d1Size);
  std::vector<bool> paired2(d2Size);
  for (int b = 0; b < d1Size; ++b) paired1[b] = false;
  for (int b = 0; b < d2Size; ++b) paired2[b] = false;

  int numberOfMismatches = 0;
  for (int m = 0, ms = (int) matchings->size(); m < ms; ++m)
  {
    matchingTuple t = matchings->at(m);
    int i = transposeOriginal ? std::get<1>(t) : std::get<0>(t);
    int j = transposeOriginal ? std::get<0>(t) : std::get<1>(t);
    dataType val = std::get<2>(t);

    diagramTuple t1 = CTDiagram1->at(i);
    diagramTuple t2 = CTDiagram2->at(j);
    dataType rX = std::get<6>(t1); dataType rY = std::get<10>(t1);
    dataType cX = std::get<6>(t2); dataType cY = std::get<10>(t2);
    dataType x = rX - cX; dataType y = rY - cY;
    paired1[i] = true;
    paired2[j] = true;
    dataType lInf = std::max(abs<dataType>(x), abs<dataType>(y));

    if (geometricalFactor == 1.0 && ((wasserstein < 0 && lInf != val) || (wasserstein > 0 && pow(lInf, wasserstein) != val)))
      ++numberOfMismatches;

    dataType partialDistance = distanceFunction(t1, t2);
    //wasserstein > 0 ? pow(lInf, wasserstein) : std::max(d, lInf);

    if (wasserstein > 0) d += partialDistance;
    else d = partialDistance;
  }

  if (numberOfMismatches > 0) {
    std::stringstream msg;
    msg << "[BottleneckDistance] Distance mismatch when rebuilding " << numberOfMismatches << " matchings." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  dataType affectationD = d;
  d = wasserstein > 0 ?
      pow(d + addedMaxPersistence + addedMinPersistence + addedSadPersistence, (1.0 / (double) wasserstein)) :
      std::max(d, std::max(addedMaxPersistence, std::max(addedMinPersistence, addedSadPersistence)));

  {
    std::stringstream msg;
    msg << "[BottleneckDistance] Computed distance diagMax("
        << addedMaxPersistence << "), diagMin("
        << addedMinPersistence << "), diagSad("
        << addedSadPersistence << "), affAll(" << affectationD << "), res(" << d << ")" << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  distance_ = (double) d;
  return 0;
}

#endif
