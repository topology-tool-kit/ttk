#include <AssignmentMunkres.h>
#include <BottleneckDistance.h>
#include <GabowTarjan.h>
#include <Geometry.h>

ttk::BottleneckDistance::BottleneckDistance() {
  this->setDebugMsgPrefix("BottleneckDistance");
}

int ttk::BottleneckDistance::execute(const ttk::DiagramType &diag0,
                                     const ttk::DiagramType &diag1,
                                     std::vector<MatchingType> &matchings) {
  Timer t;

  if(diag0.empty() || diag1.empty()) {
    this->printErr("Empty input diagrams");
    return -1;
  }

  const bool fromParaView = this->PVAlgorithm >= 0;
  if(fromParaView) {
    switch(this->PVAlgorithm) {
      case 0:
        this->printMsg("Solving with the TTK approach");
        this->computeBottleneck(diag0, diag1, matchings);
        break;
      case 1: {
        this->printMsg("Solving with the legacy Dionysus exact approach.");
        this->printErr("Not supported");
      } break;
      case 2: {
        this->printMsg(
          "Solving with the approximate Dionysus geometric approach.");
        this->printErr("Not supported");
      } break;
      case 3: {
        this->printMsg("Solving with the parallel TTK approach");
        this->printErr("Not supported");
      } break;
      case 4: {
        this->printMsg("Benchmarking");
        this->printErr("Not supported");
      } break;
      default: {
        this->printErr("You must specify a valid assignment algorithm.");
      }
    }

  } else {
    switch(str2int(this->DistanceAlgorithm.c_str())) {
      case str2int("0"):
      case str2int("ttk"):
        this->printMsg("Solving with the TTK approach");
        this->computeBottleneck(diag0, diag1, matchings);
        break;
      case str2int("1"):
      case str2int("legacy"): {
        this->printMsg("Solving with the legacy Dionysus exact approach.");
        this->printErr("Not supported");
      } break;
      case str2int("2"):
      case str2int("geometric"): {
        this->printMsg(
          "Solving with the approximate Dionysus geometric approach.");
        this->printErr("Not supported");
      } break;
      case str2int("3"):
      case str2int("parallel"): {
        this->printMsg("Solving with the parallel TTK approach");
        this->printErr("Not supported");
      } break;
      case str2int("bench"): {
        this->printMsg("Benchmarking");
        this->printErr("Not supported");
      } break;
      default: {
        this->printErr("You must specify a valid assignment algorithm.");
      }
    }
  }

  this->printMsg("Complete", 1, t.getElapsedTime(), threadNumber_);

  return 0;
}

void ttk::BottleneckDistance::buildCostMatrices(
  const ttk::DiagramType &CTDiagram1,
  const ttk::DiagramType &CTDiagram2,
  const double zeroThresh,
  std::vector<std::vector<double>> &minMatrix,
  std::vector<std::vector<double>> &maxMatrix,
  std::vector<std::vector<double>> &sadMatrix,
  const bool reverseMin,
  const bool reverseMax,
  const bool reverseSad,
  const int wasserstein) const {

  int maxI = 0, minI = 0;
  int maxJ = 0, minJ = 0;
  int sadI = 0, sadJ = 0;

  for(const auto &p1 : CTDiagram1) {
    if(std::abs(p1.persistence()) < zeroThresh)
      continue;

    bool isMin1 = (p1.birth.type == CriticalType::Local_minimum
                   || p1.death.type == CriticalType::Local_minimum);
    bool isMax1 = (p1.birth.type == CriticalType::Local_maximum
                   || p1.death.type == CriticalType::Local_maximum);
    const bool isSad1 = (p1.birth.type == CriticalType::Saddle1
                         && p1.death.type == CriticalType::Saddle2)
                        || (p1.birth.type == CriticalType::Saddle2
                            && p1.death.type == CriticalType::Saddle1);
    if(p1.birth.type == CriticalType::Local_minimum
       && p1.death.type == CriticalType::Local_maximum) {
      isMin1 = false;
      isMax1 = true;
    }

    minJ = 0;
    maxJ = 0;
    sadJ = 0;

    for(const auto &p2 : CTDiagram2) {
      if(std::abs(p2.persistence()) < zeroThresh)
        continue;

      bool isMin2 = (p2.birth.type == CriticalType::Local_minimum
                     || p2.death.type == CriticalType::Local_minimum);
      bool isMax2 = (p2.birth.type == CriticalType::Local_maximum
                     || p2.death.type == CriticalType::Local_maximum);
      const bool isSad2 = (p2.birth.type == CriticalType::Saddle1
                           && p2.death.type == CriticalType::Saddle2)
                          || (p2.birth.type == CriticalType::Saddle2
                              && p2.death.type == CriticalType::Saddle1);
      if(p2.birth.type == CriticalType::Local_minimum
         && p2.death.type == CriticalType::Local_maximum) {
        isMin2 = false;
        isMax2 = true;
      }
      if((isMin1 && !isMin2) || (isMax1 && !isMax2) || (isSad1 && !isSad2))
        continue;

      double distance = this->distanceFunction(p1, p2, wasserstein);
      const double diag1 = this->diagonalDistanceFunction(p1, wasserstein);
      const double diag2 = this->diagonalDistanceFunction(p2, wasserstein);

      if(distance > diag1 + diag2)
        distance = std::numeric_limits<double>::max();

      if(isMin1 && isMin2) {
        if(reverseMin)
          minMatrix[minJ++][minI] = distance;
        else
          minMatrix[minI][minJ++] = distance;
      } else if(isMax1 && isMax2) {
        if(reverseMax)
          maxMatrix[maxJ++][maxI] = distance;
        else
          maxMatrix[maxI][maxJ++] = distance;
      } else if(isSad1 && isSad2) {
        if(reverseSad)
          sadMatrix[sadJ++][sadI] = distance;
        else
          sadMatrix[sadI][sadJ++] = distance;
      }
    }

    const double distanceToDiagonal
      = this->diagonalDistanceFunction(p1, wasserstein);
    if(isMin1) {
      if(reverseMin)
        minMatrix[minJ++][minI] = distanceToDiagonal;
      else
        minMatrix[minI][minJ++] = distanceToDiagonal;
    }
    if(isMax1) {
      if(reverseMax)
        maxMatrix[maxJ++][maxI] = distanceToDiagonal;
      else
        maxMatrix[maxI][maxJ++] = distanceToDiagonal;
    }
    if(isSad1) {
      if(reverseSad)
        sadMatrix[sadJ++][sadI] = distanceToDiagonal;
      else
        sadMatrix[sadI][sadJ++] = distanceToDiagonal;
    }

    if(isMin1)
      ++minI;
    if(isMax1)
      ++maxI;
    if(isSad1)
      ++sadI;
  }

  minJ = 0;
  maxJ = 0;
  sadJ = 0;

  // Last row: match remaining J components with diagonal.
  for(const auto &p3 : CTDiagram2) {
    if(std::abs(p3.persistence()) < zeroThresh)
      continue;

    bool isMin2 = (p3.birth.type == CriticalType::Local_minimum
                   || p3.death.type == CriticalType::Local_minimum);
    bool isMax2 = (p3.birth.type == CriticalType::Local_maximum
                   || p3.death.type == CriticalType::Local_maximum);
    const bool isSad2 = (p3.birth.type == CriticalType::Saddle1
                         && p3.death.type == CriticalType::Saddle2)
                        || (p3.birth.type == CriticalType::Saddle2
                            && p3.death.type == CriticalType::Saddle1);
    if(p3.birth.type == CriticalType::Local_minimum
       && p3.death.type == CriticalType::Local_maximum) {
      isMin2 = false;
      isMax2 = true;
    }

    const double distanceToDiagonal
      = this->diagonalDistanceFunction(p3, wasserstein);
    if(isMin2) {
      if(reverseMin)
        minMatrix[minJ++][minI] = distanceToDiagonal;
      else
        minMatrix[minI][minJ++] = distanceToDiagonal;
    }
    if(isMax2) {
      if(reverseMax)
        maxMatrix[maxJ++][maxI] = distanceToDiagonal;
      else
        maxMatrix[maxI][maxJ++] = distanceToDiagonal;
    }
    if(isSad2) {
      if(reverseSad)
        sadMatrix[sadJ++][sadI] = distanceToDiagonal;
      else
        sadMatrix[sadI][sadJ++] = distanceToDiagonal;
    }
  }

  // Last cell
  if(reverseMin)
    minMatrix[minJ][minI] = std::numeric_limits<double>::max();
  else
    minMatrix[minI][minJ] = std::numeric_limits<double>::max();

  if(reverseMax)
    maxMatrix[maxJ][maxI] = std::numeric_limits<double>::max();
  else
    maxMatrix[maxI][maxJ] = std::numeric_limits<double>::max();

  if(reverseSad)
    sadMatrix[sadJ][sadI] = std::numeric_limits<double>::max();
  else
    sadMatrix[sadI][sadJ] = std::numeric_limits<double>::max();
}

double ttk::BottleneckDistance::computeGeometricalRange(
  const ttk::DiagramType &CTDiagram1,
  const ttk::DiagramType &CTDiagram2) const {

  float minX, minY, minZ, maxX, maxY, maxZ;

  std::array<float, 3> min1{
    std::numeric_limits<float>::max(),
    std::numeric_limits<float>::max(),
    std::numeric_limits<float>::max(),
  },
    min2 = min1;
  std::array<float, 3> max1{
    std::numeric_limits<float>::lowest(),
    std::numeric_limits<float>::lowest(),
    std::numeric_limits<float>::lowest(),
  },
    max2 = max1;

  for(const auto &p : CTDiagram1) {
    min1[0] = std::min(std::min(min1[0], p.birth.coords[0]), p.death.coords[0]);
    min1[1] = std::min(std::min(min1[1], p.birth.coords[1]), p.death.coords[1]);
    min1[2] = std::min(std::min(min1[2], p.birth.coords[2]), p.death.coords[2]);
    max1[0] = std::max(std::max(max1[0], p.birth.coords[0]), p.birth.coords[0]);
    max1[1] = std::max(std::max(max1[1], p.birth.coords[1]), p.birth.coords[1]);
    max1[2] = std::max(std::max(max1[2], p.birth.coords[2]), p.birth.coords[2]);
  }

  for(const auto &p : CTDiagram2) {
    min2[0] = std::min(std::min(min2[0], p.birth.coords[0]), p.death.coords[0]);
    min2[1] = std::min(std::min(min2[1], p.birth.coords[1]), p.death.coords[1]);
    min2[2] = std::min(std::min(min2[2], p.birth.coords[2]), p.death.coords[2]);
    max2[0] = std::max(std::max(max2[0], p.birth.coords[0]), p.birth.coords[0]);
    max2[1] = std::max(std::max(max2[1], p.birth.coords[1]), p.birth.coords[1]);
    max2[2] = std::max(std::max(max2[2], p.birth.coords[2]), p.birth.coords[2]);
  }

  minX = std::min(min1[0], min2[0]);
  minY = std::min(min1[1], min2[1]);
  minZ = std::min(min1[2], min2[2]);
  maxX = std::max(max1[0], max2[0]);
  maxY = std::max(max1[1], max2[1]);
  maxZ = std::max(max1[2], max2[2]);

  const auto square = [](const double a) { return a * a; };

  return std::sqrt(square(maxX - minX) + square(maxY - minY)
                   + square(maxZ - minZ));
}

double ttk::BottleneckDistance::computeMinimumRelevantPersistence(
  const ttk::DiagramType &CTDiagram1,
  const ttk::DiagramType &CTDiagram2) const {

  const auto sp = this->Tolerance;
  const double s = sp > 0.0 && sp < 100.0 ? sp / 100.0 : 0;

  std::vector<double> toSort(CTDiagram1.size() + CTDiagram2.size());
  for(size_t i = 0; i < CTDiagram1.size(); ++i) {
    const auto &t = CTDiagram1[i];
    toSort[i] = std::abs(t.persistence());
  }
  for(size_t i = 0; i < CTDiagram2.size(); ++i) {
    const auto &t = CTDiagram2[i];
    toSort[CTDiagram1.size() + i] = std::abs(t.persistence());
  }

  const auto minVal = *std::min_element(toSort.begin(), toSort.end());
  const auto maxVal = *std::max_element(toSort.begin(), toSort.end());
  return s * (maxVal - minVal);
}

void ttk::BottleneckDistance::computeMinMaxSaddleNumberAndMapping(
  const ttk::DiagramType &CTDiagram,
  int &nbMin,
  int &nbMax,
  int &nbSaddle,
  std::vector<int> &minMap,
  std::vector<int> &maxMap,
  std::vector<int> &sadMap,
  const double zeroThresh) const {

  for(size_t i = 0; i < CTDiagram.size(); ++i) {
    const auto &t = CTDiagram[i];
    const auto nt1 = t.birth.type;
    const auto nt2 = t.death.type;
    const auto dt = t.persistence();
    if(std::abs(dt) < zeroThresh)
      continue;

    if(nt1 == CriticalType::Local_minimum
       && nt2 == CriticalType::Local_maximum) {
      nbMax++;
      maxMap.push_back(i);
    } else {
      if(nt1 == CriticalType::Local_maximum
         || nt2 == CriticalType::Local_maximum) {
        nbMax++;
        maxMap.push_back(i);
      }
      if(nt1 == CriticalType::Local_minimum
         || nt2 == CriticalType::Local_minimum) {
        nbMin++;
        minMap.push_back(i);
      }
      if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
         || (nt1 == CriticalType::Saddle2 && nt2 == CriticalType::Saddle1)) {
        nbSaddle++;
        sadMap.push_back(i);
      }
    }
  }
}

namespace {
  void solvePWasserstein(std::vector<std::vector<double>> &matrix,
                         std::vector<ttk::MatchingType> &matchings,
                         ttk::AssignmentMunkres<double> &solver) {

    solver.setInput(matrix);
    solver.run(matchings);
    solver.clearMatrix();
  }

  void solveInfinityWasserstein(const int nbRow,
                                const int nbCol,
                                std::vector<std::vector<double>> &matrix,
                                std::vector<ttk::MatchingType> &matchings,
                                ttk::GabowTarjan &solver) {

    // Copy input matrix.
    auto bottleneckMatrix = matrix;

    // Solve.
    solver.setInput(nbRow, nbCol, &bottleneckMatrix);
    solver.run(matchings);
    solver.clear();
  }
} // namespace

double ttk::BottleneckDistance::buildMappings(
  const std::vector<MatchingType> &inputMatchings,
  const bool transposeGlobal,
  const bool transposeLocal,
  std::vector<MatchingType> &outputMatchings,
  const std::vector<int> &m1,
  const std::vector<int> &m2,
  int wasserstein) const {

  // Input map permutation (so as to ignore transposition later on)
  const auto &map1 = transposeLocal ? m2 : m1;
  const auto &map2 = transposeLocal ? m1 : m2;

  double addedPersistence = 0;
  const auto doTranspose = transposeGlobal ^ transposeLocal;

  for(const auto &t : inputMatchings) {

    const auto val = std::abs(std::get<2>(t));
    const size_t p1 = std::get<0>(t);
    const size_t p2 = std::get<1>(t);

    if(p1 >= map1.size() || p2 >= map2.size()) {
      addedPersistence = (wasserstein > 0 ? addedPersistence + val
                                          : std::max(val, addedPersistence));
    } else {
      if(doTranspose) {
        outputMatchings.emplace_back(map2[p2], map1[p1], val);
      } else {
        outputMatchings.emplace_back(map1[p1], map2[p2], val);
      }
    }
  }

  return addedPersistence;
}

double ttk::BottleneckDistance::distanceFunction(const ttk::PersistencePair &a,
                                                 const ttk::PersistencePair &b,
                                                 const int wasserstein) const {

  const int w = std::max(wasserstein, 1); // L_inf not managed.

  // We don't match critical points of different index.
  // This must be ensured before calling the distance function.
  const bool isMin1 = a.birth.type == CriticalType::Local_minimum;
  const bool isMax1 = a.death.type == CriticalType::Local_maximum;

  const std::array<float, 3> coordsAbsDiff{
    std::abs(a.death.coords[0] - b.death.coords[0]),
    std::abs(a.death.coords[1] - b.death.coords[1]),
    std::abs(a.death.coords[2] - b.death.coords[2]),
  };

  const auto x
    = ((isMin1 && !isMax1) ? this->PE : this->PS)
      * Geometry::pow(std::abs(a.birth.sfValue - b.birth.sfValue), w);
  const auto y
    = (isMax1 ? this->PE : this->PS)
      * Geometry::pow(std::abs(a.death.sfValue - b.death.sfValue), w);
  const double geoDistance
    = isMax1 || isMin1
        ? (this->PX * Geometry::pow(coordsAbsDiff[0], w)
           + this->PY * Geometry::pow(coordsAbsDiff[1], w)
           + this->PZ * Geometry::pow(coordsAbsDiff[2], w))
        : (this->PX
             * Geometry::pow(
               std::abs(a.birth.coords[0] + a.death.coords[0]) / 2
                 - std::abs(b.birth.coords[0] + b.death.coords[0]) / 2,
               w)
           + this->PY
               * Geometry::pow(
                 std::abs(a.birth.coords[1] + a.death.coords[1]) / 2
                   - std::abs(b.birth.coords[1] + b.death.coords[1]) / 2,
                 w)
           + this->PZ
               * Geometry::pow(
                 std::abs(a.birth.coords[2] + a.death.coords[2]) / 2
                   - std::abs(b.birth.coords[2] + b.death.coords[2]) / 2,
                 w));

  const double persDistance = x + y;
  return Geometry::pow(persDistance + geoDistance, 1.0 / w);
}

double ttk::BottleneckDistance::diagonalDistanceFunction(
  const ttk::PersistencePair &a, const int wasserstein) const {

  const int w = std::max(wasserstein, 1);
  const bool isMin1 = a.birth.type == CriticalType::Local_minimum;
  const bool isMax1 = a.death.type == CriticalType::Local_maximum;
  const double infDistance = (isMin1 || isMax1 ? this->PE : this->PS)
                             * Geometry::pow(std::abs(a.persistence()), w);
  const double geoDistance
    = (this->PX
         * Geometry::pow(std::abs(a.death.coords[0] - a.birth.coords[0]), w)
       + this->PY
           * Geometry::pow(std::abs(a.death.coords[1] - a.birth.coords[1]), w)
       + this->PZ
           * Geometry::pow(std::abs(a.death.coords[2] - a.birth.coords[2]), w));

  return Geometry::pow(infDistance + geoDistance, 1.0 / w);
}

int ttk::BottleneckDistance::computeBottleneck(
  const ttk::DiagramType &d1,
  const ttk::DiagramType &d2,
  std::vector<MatchingType> &matchings) {

  const auto transposeOriginal = d1.size() > d2.size();
  if(transposeOriginal) {
    this->printMsg("The first persistence diagram is larger than the second.");
    this->printMsg("Solving the transposed problem.");
  }
  const auto &CTDiagram1 = transposeOriginal ? d2 : d1;
  const auto &CTDiagram2 = transposeOriginal ? d1 : d2;

  // Check user parameters.
  const auto isBottleneck = this->WassersteinMetric == "inf";
  const int wasserstein = isBottleneck ? -1 : stoi(this->WassersteinMetric);
  if(wasserstein < 0 && !isBottleneck) {
    this->printErr("Wrong value for the Wassertein power parameter");
    return -4;
  }

  // Needed to limit computation time.
  const auto zeroThresh
    = this->computeMinimumRelevantPersistence(CTDiagram1, CTDiagram2);

  // Initialize solvers.
  std::vector<MatchingType> minMatchings;
  std::vector<MatchingType> maxMatchings;
  std::vector<MatchingType> sadMatchings;

  // Initialize cost matrices.
  int nbRowMin = 0, nbColMin = 0;
  int nbRowMax = 0, nbColMax = 0;
  int nbRowSad = 0, nbColSad = 0;

  // Remap for matchings.
  std::array<std::vector<int>, 3> map1{};
  std::array<std::vector<int>, 3> map2{};

  this->computeMinMaxSaddleNumberAndMapping(CTDiagram1, nbRowMin, nbRowMax,
                                            nbRowSad, map1[0], map1[2], map1[1],
                                            zeroThresh);
  this->computeMinMaxSaddleNumberAndMapping(CTDiagram2, nbColMin, nbColMax,
                                            nbColSad, map2[0], map2[2], map2[1],
                                            zeroThresh);

  // Automatically transpose if nb rows > nb cols
  const auto maxRowColMin = std::max(nbRowMin + 1, nbColMin + 1);
  const auto maxRowColMax = std::max(nbRowMax + 1, nbColMax + 1);
  const auto maxRowColSad = std::max(nbRowSad + 1, nbColSad + 1);

  const auto minRowColMin = std::min(nbRowMin + 1, nbColMin + 1);
  const auto minRowColMax = std::min(nbRowMax + 1, nbColMax + 1);
  const auto minRowColSad = std::min(nbRowSad + 1, nbColSad + 1);

  std::vector<std::vector<double>> minMatrix(
    minRowColMin, std::vector<double>(maxRowColMin));
  std::vector<std::vector<double>> maxMatrix(
    minRowColMax, std::vector<double>(maxRowColMax));
  std::vector<std::vector<double>> sadMatrix(
    minRowColSad, std::vector<double>(maxRowColSad));

  const bool transposeMin = nbRowMin > nbColMin;
  const bool transposeMax = nbRowMax > nbColMax;
  const bool transposeSad = nbRowSad > nbColSad;

  Timer t;

  this->buildCostMatrices(CTDiagram1, CTDiagram2, zeroThresh, minMatrix,
                          maxMatrix, sadMatrix, transposeMin, transposeMax,
                          transposeSad, wasserstein);

  if(!isBottleneck) {

    if(nbRowMin > 0 && nbColMin > 0) {
      AssignmentMunkres<double> solverMin;
      this->printMsg("Affecting minima...");
      solvePWasserstein(minMatrix, minMatchings, solverMin);
    }

    if(nbRowMax > 0 && nbColMax > 0) {
      AssignmentMunkres<double> solverMax;
      this->printMsg("Affecting maxima...");
      solvePWasserstein(maxMatrix, maxMatchings, solverMax);
    }

    if(nbRowSad > 0 && nbColSad > 0) {
      AssignmentMunkres<double> solverSad;
      this->printMsg("Affecting saddles...");
      solvePWasserstein(sadMatrix, sadMatchings, solverSad);
    }

  } else {

    // Launch solving for minima.
    if(nbRowMin > 0 && nbColMin > 0) {
      GabowTarjan solverMin;
      this->printMsg("Affecting minima...");
      solveInfinityWasserstein(
        minRowColMin, maxRowColMin, minMatrix, minMatchings, solverMin);
    }

    // Launch solving for maxima.
    if(nbRowMax > 0 && nbColMax > 0) {
      GabowTarjan solverMax;
      this->printMsg("Affecting maxima...");
      solveInfinityWasserstein(
        minRowColMax, maxRowColMax, maxMatrix, maxMatchings, solverMax);
    }

    // Launch solving for saddles.
    if(nbRowSad > 0 && nbColSad > 0) {
      GabowTarjan solverSad;
      this->printMsg("Affecting saddles...");
      solveInfinityWasserstein(
        minRowColSad, maxRowColSad, sadMatrix, sadMatchings, solverSad);
    }
  }

  this->printMsg("TTK CORE DONE", 1, t.getElapsedTime());

  // Rebuild mappings.
  // Begin cost computation for unpaired vertices.
  const std::array<double, 3> addedPersistence{
    this->buildMappings(minMatchings, transposeOriginal, transposeMin,
                        matchings, map1[0], map2[0], wasserstein),
    this->buildMappings(sadMatchings, transposeOriginal, transposeSad,
                        matchings, map1[1], map2[1], wasserstein),
    this->buildMappings(maxMatchings, transposeOriginal, transposeMax,
                        matchings, map1[2], map2[2], wasserstein),
  };

  // TODO [HIGH] do that for embeddings
  // Recompute matching weights for user-friendly distance.
  std::array<double, 3> costs{};
  std::vector<bool> paired1(CTDiagram1.size(), false);
  std::vector<bool> paired2(CTDiagram2.size(), false);

  // store last matching distance for bottleneck
  double partialDistance{};

  for(const auto &mt : matchings) {
    const int i = transposeOriginal ? std::get<1>(mt) : std::get<0>(mt);
    const int j = transposeOriginal ? std::get<0>(mt) : std::get<1>(mt);

    const auto &t1 = CTDiagram1[i];
    const auto &t2 = CTDiagram2[j];
    paired1[i] = true;
    paired2[j] = true;

    partialDistance = this->distanceFunction(t1, t2, wasserstein);

    if(t1.death.type == CriticalType::Local_maximum) {
      if(!isBottleneck) {
        costs[2] += partialDistance;
      } else {
        costs[2] = std::max(costs[2], partialDistance);
      }
    } else if(t1.birth.type == CriticalType::Local_minimum) {
      if(!isBottleneck) {
        costs[0] += partialDistance;
      } else {
        costs[0] = std::max(costs[0], partialDistance);
      }
    } else if(t1.birth.type == CriticalType::Saddle1
              && t1.death.type == CriticalType::Saddle2) {
      if(!isBottleneck) {
        costs[1] += partialDistance;
      } else {
        costs[1] = std::max(costs[1], partialDistance);
      }
    }
  }

  const auto affectationD
    = !isBottleneck ? costs[0] + costs[1] + costs[2] : partialDistance;
  const auto addedPers
    = addedPersistence[0] + addedPersistence[1] + addedPersistence[2];
  this->distance_
    = !isBottleneck
        ? Geometry::pow(affectationD + addedPers, 1.0 / wasserstein)
        : std::max(affectationD, *std::max_element(addedPersistence.begin(),
                                                   addedPersistence.end()));

  this->printMsg("Computed distance:");
  this->printMsg("diagMax(" + std::to_string(addedPersistence[2])
                 + "), diagMin(" + std::to_string(addedPersistence[0])
                 + "), diagSad(" + std::to_string(addedPersistence[1]) + ")");
  this->printMsg("affAll(" + std::to_string(affectationD) + "), res("
                 + std::to_string(this->distance_) + ")");

  // aggregate costs per pair type
  if(!isBottleneck) {
    costs[0] = Geometry::pow(costs[0] + addedPersistence[0], 1.0 / wasserstein);
    costs[1] = Geometry::pow(costs[1] + addedPersistence[1], 1.0 / wasserstein);
    costs[2] = Geometry::pow(costs[2] + addedPersistence[2], 1.0 / wasserstein);
  } else {
    costs[0] = std::max(costs[0], addedPersistence[0]);
    costs[1] = std::max(costs[1], addedPersistence[1]);
    costs[2] = std::max(costs[2], addedPersistence[2]);
  }

  this->costs_ = costs;

  // display results
  const std::vector<std::vector<std::string>> rows{
    {" Min-saddle cost", std::to_string(this->costs_[0])},
    {" Saddle-saddle cost", std::to_string(this->costs_[1])},
    {" Saddle-max cost", std::to_string(this->costs_[2])},
    {isBottleneck ? "Bottleneck Distance" : "Wasserstein Distance",
     std::to_string(this->distance_)},
  };
  this->printMsg(rows);

  return 0;
}
