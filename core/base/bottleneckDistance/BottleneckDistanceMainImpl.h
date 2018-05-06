#ifndef _BOTTLENECKDISTANCEIMPLMAIN_H
#define _BOTTLENECKDISTANCEIMPLMAIN_H

//  vector <   -- diagram
//    tuple <    -- pair of critical points
//      idVertex, NodeType
//      idVertex, NodeType
//      dataType  -- persistance of pair
//      idVertex  -- type (0/min, 1/saddle, 2/max)
//      dataType            -- scalar value at vertex 1
//      float, float, float -- vertex 1 coordinates
//      dataType            -- scalar value at vertex 2
//      float, float, float -- vertex 2 coordinates
template <typename dataType>
int ttk::BottleneckDistance::computeBottleneck(
  const std::vector<diagramTuple> *CTDiagram1,
  const std::vector<diagramTuple> *CTDiagram2,
  std::vector<matchingTuple> *matchings,
  const bool usePersistenceMetric,
  const double alpha)
{
  auto* distance = new dataType;

  const int d1Size = (int) CTDiagram1->size();
  const int d2Size = (int) CTDiagram2->size();

  if (d1Size > d2Size) {
    std::stringstream msg;
    msg << "[BottleneckDistance] The first persistence diagram is larger than the second." << std::endl;
    msg << "[BottleneckDistance] You should consider switching them for performance." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  // Check user parameters.
  const double geometricalFactor = (alpha < 0.0 || alpha > 1.0) ? 1.0 : alpha;
  const int wasserstein = (wasserstein_ == "inf") ? -1 : stoi(wasserstein_);
  if (wasserstein < 0 && wasserstein != -1) return -4;

  // Needed for distance computation.
  const double maxDistance =
    this->computeGeometricalRange(CTDiagram1, CTDiagram2, d1Size, d2Size);

  // Needed to limit computation time.
  const dataType zeroThresh =
    this->computeMinimumRelevantPersistence(CTDiagram1, CTDiagram2, d1Size, d2Size);

  // Initialize solvers.
  std::vector<matchingTuple> minMatchings;
  std::vector<matchingTuple> maxMatchings;
  std::vector<matchingTuple> sadMatchings;

  Munkres solverMin;
  Munkres solverMax;
  Munkres solverSad;

  // Initialize cost matrices.
  int nbRowMin = 0, nbColMin = 0;
  int nbRowMax = 0, nbColMax = 0;
  int nbRowSad = 0, nbColSad = 0;

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

  auto** minMatrix = new dataType*[nbRowMin+1];
  auto** maxMatrix = new dataType*[nbRowMax+1];
  auto** sadMatrix = new dataType*[nbRowSad+1];
  for (int i = 0; i < nbRowMin+1; ++i) minMatrix[i] = new dataType[nbColMin+1];
  for (int i = 0; i < nbRowMax+1; ++i) maxMatrix[i] = new dataType[nbColMax+1];
  for (int i = 0; i < nbRowSad+1; ++i) sadMatrix[i] = new dataType[nbColSad+1];

  std::function<dataType (const diagramTuple, const diagramTuple)>
    distanceFunction =
      [maxDistance, wasserstein, geometricalFactor] (const diagramTuple a, const diagramTuple b) -> dataType {
        dataType rX = std::get<6>(a);
        dataType rY = std::get<10>(a);
        dataType cX = std::get<6>(b);
        dataType cY = std::get<10>(b);
        dataType x = abs_diff<dataType>(rX, cX);
        dataType y = abs_diff<dataType>(rY, cY);
        double dist = sqrt(
          pow((std::get<7>(a)+std::get<11>(a))/2 - (std::get<7>(b)+std::get<11>(b))/2, 2) +
          pow((std::get<8>(a)+std::get<12>(a))/2 - (std::get<8>(b)+std::get<12>(b))/2, 2) +
          pow((std::get<9>(a)+std::get<13>(a))/2 - (std::get<9>(b)+std::get<13>(b))/2, 2)
        );
        dist /= maxDistance;
        double val =
          (wasserstein > 0
           ? pow(x, wasserstein) + pow(y, wasserstein)  // Wasserstein
           : std::max(x, y) // Bottleneck
          );
        return geometricalFactor * val + (1.0 - geometricalFactor) * dist;
      };

  std::function<dataType (const diagramTuple)>
    diagonalDistanceFunction =
      [wasserstein](const diagramTuple a) -> dataType {
        dataType rX = std::get<6>(a);
        dataType rY = std::get<10>(a);
        dataType dist = abs_diff<dataType>(rX, rY);
        return
          wasserstein > 0 ? pow(dist, wasserstein) : dist;
        //abs<dataType>(std::get<4>(a))
        // abs<dataType>(std::get<4>(a));
      };

  this->buildCostMatrices(
    CTDiagram1, CTDiagram2, d1Size, d2Size, distanceFunction, diagonalDistanceFunction,
    zeroThresh, minMatrix, maxMatrix, sadMatrix);


  if (wasserstein > 0) {

    if (nbRowMin > 0 && nbColMin > 0) {
      dMsg(std::cout, "[BottleneckDistance] Affecting minima...\n", timeMsg);
      this->solvePWasserstein(nbRowMin, nbColMin, minMatrix, &minMatchings, &solverMin);
    }

    if (nbRowMax > 0 && nbColMax > 0) {
      dMsg(std::cout, "[BottleneckDistance] Affecting maxima...\n", timeMsg);
      this->solvePWasserstein(nbRowMax, nbColMax, maxMatrix, &maxMatchings, &solverMax);
    }

    if (nbRowSad > 0 && nbColSad > 0) {
      dMsg(std::cout, "[BottleneckDistance] Affecting saddles...\n", timeMsg);
      this->solvePWasserstein(nbRowSad, nbColSad, sadMatrix, &sadMatchings, &solverSad);
    }

  } else {

    // Launch solving for minima.
    if (nbRowMin > 0 && nbColMin > 0) {
      dMsg(std::cout, "[BottleneckDistance] Affecting minima...\n", timeMsg);
      this->solveInfinityWasserstein(nbRowMin, nbColMin, minMatrix, &minMatchings, &solverMin);
    }

    // Launch solving for maxima.
    if (nbRowMax > 0 && nbColMax > 0) {
      dMsg(std::cout, "[BottleneckDistance] Affecting maxima...\n", timeMsg);
      this->solveInfinityWasserstein(nbRowMax, nbColMax, maxMatrix, &maxMatchings, &solverMax);
    }

    // Launch solving for saddles.
    if (nbRowSad > 0 && nbColSad > 0) {
      dMsg(std::cout, "[BottleneckDistance] Affecting saddles...\n", timeMsg);
      this->solveInfinityWasserstein(nbRowSad, nbColSad, sadMatrix, &sadMatchings, &solverSad);
    }
  }

  // Rebuild mappings.
  // Begin cost computation for unpaired vertices.
  dataType addedMinPersistence =
    this->buildMappings(minMatchings, matchings, minMap1, minMap2, wasserstein);

  dataType addedMaxPersistence =
    this->buildMappings(maxMatchings, matchings, maxMap1, maxMap2, wasserstein);

  dataType addedSadPersistence =
    this->buildMappings(sadMatchings, matchings, sadMap1, sadMap2, wasserstein);

  // Recompute matching weights for user-fristd::endly distance.
  dataType d = 0;
  std::vector<bool> paired1(d1Size);
  std::vector<bool> paired2(d2Size);
  for (int b = 0; b < d1Size; ++b) paired1[b] = false;
  for (int b = 0; b < d2Size; ++b) paired2[b] = false;

  for (int m = 0, ms = (int) matchings->size(); m < ms; ++m)
  {
    matchingTuple t = matchings->at(m);
    int i = std::get<0>(t);
    int j = std::get<1>(t);

    diagramTuple t1 = CTDiagram1->at(i);
    diagramTuple t2 = CTDiagram2->at(j);
    dataType rX = std::get<6>(t1);
    dataType rY = std::get<10>(t1);
    dataType cX = std::get<6>(t2);
    dataType cY = std::get<10>(t2);
    dataType x = rX - cX;
    dataType y = rY - cY;

    paired1[i] = true;
    paired2[j] = true;
    dataType linfty = std::max(abs<dataType>(x), abs<dataType>(y));
    if (wasserstein > 0) {
      d += pow(linfty, wasserstein);
    } else {
      d = std::max(d, linfty);
    }
  }

  d = wasserstein > 0
      ? pow(d + addedMaxPersistence + addedMinPersistence + addedSadPersistence, (1.0 / (double) wasserstein))
      : std::max(d, std::max(addedMaxPersistence, std::max(addedMinPersistence, addedSadPersistence)));

  {
    std::stringstream msg;
    msg << "[BottleneckDistance] Computed distance " << d << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  *distance = d;
  distance_ = (void*)(distance);
  return 0;
}

#endif
