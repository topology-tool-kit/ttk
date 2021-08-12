#include <AssignmentMunkres.h>
#include <BottleneckDistance.h>
#include <GabowTarjan.h>
#include <Geometry.h>

using ttk::BottleneckDistance;
using ttk::MatchingType;

BottleneckDistance::BottleneckDistance() {
  this->setDebugMsgPrefix("BottleneckDistance");
}

constexpr unsigned long long str2int(const char *str, int h = 0) {
  return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
}

int BottleneckDistance::execute(const ttk::DiagramType &diag0,
                                const ttk::DiagramType &diag1,
                                std::vector<MatchingType> &matchings,
                                const bool usePersistenceMetric) {
  Timer t;

  bool fromParaView = pvAlgorithm_ >= 0;
  if(fromParaView) {
    switch(pvAlgorithm_) {
      case 0:
        this->printMsg("Solving with the TTK approach");
        this->computeBottleneck(diag0, diag1, matchings, usePersistenceMetric);
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
    switch(str2int(algorithm_.c_str())) {
      case str2int("0"):
      case str2int("ttk"):
        this->printMsg("Solving with the TTK approach");
        this->computeBottleneck(diag0, diag1, matchings, usePersistenceMetric);
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

template <typename distFuncType, typename diagFuncType>
void ttk::BottleneckDistance::buildCostMatrices(
  const ttk::DiagramType &CTDiagram1,
  const ttk::DiagramType &CTDiagram2,
  const int d1Size,
  const int d2Size,
  const distFuncType &distanceFunction,
  const diagFuncType &diagonalDistanceFunction,
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

  for(int i = 0; i < d1Size; ++i) {
    const auto &p1 = CTDiagram1[i];
    if(std::abs(p1.persistence) < zeroThresh)
      continue;

    bool isMin1 = (p1.birth.type == CriticalType::Local_minimum
                   || p1.death.type == CriticalType::Local_minimum);
    bool isMax1 = (p1.birth.type == CriticalType::Local_maximum
                   || p1.death.type == CriticalType::Local_maximum);
    bool isSad1 = (p1.birth.type == CriticalType::Saddle1
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

    for(int j = 0; j < d2Size; ++j) {
      const auto &p2 = CTDiagram2[j];
      if(std::abs(p2.persistence) < zeroThresh)
        continue;

      bool isMin2 = (p2.birth.type == CriticalType::Local_minimum
                     || p2.death.type == CriticalType::Local_minimum);
      bool isMax2 = (p2.birth.type == CriticalType::Local_maximum
                     || p2.death.type == CriticalType::Local_maximum);
      bool isSad2 = (p2.birth.type == CriticalType::Saddle1
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

      double distance = distanceFunction(p1, p2);
      double diag1 = diagonalDistanceFunction(p1);
      double diag2 = diagonalDistanceFunction(p2);

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

    double distanceToDiagonal = diagonalDistanceFunction(p1);
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
  for(int j = 0; j < d2Size; ++j) {
    const auto &p3 = CTDiagram2[j];
    if(std::abs(p3.persistence) < zeroThresh)
      continue;

    bool isMin2 = (p3.birth.type == CriticalType::Local_minimum
                   || p3.death.type == CriticalType::Local_minimum);
    bool isMax2 = (p3.birth.type == CriticalType::Local_maximum
                   || p3.death.type == CriticalType::Local_maximum);
    bool isSad2 = (p3.birth.type == CriticalType::Saddle1
                   && p3.death.type == CriticalType::Saddle2)
                  || (p3.birth.type == CriticalType::Saddle2
                      && p3.death.type == CriticalType::Saddle1);
    if(p3.birth.type == CriticalType::Local_minimum
       && p3.death.type == CriticalType::Local_maximum) {
      isMin2 = false;
      isMax2 = true;
    }

    double distanceToDiagonal = diagonalDistanceFunction(p3);
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
  {
    if(reverseMin)
      minMatrix[minJ][minI] = std::numeric_limits<double>::max();
    else
      minMatrix[minI][minJ] = std::numeric_limits<double>::max();
  }
  {
    if(reverseMax)
      maxMatrix[maxJ][maxI] = std::numeric_limits<double>::max();
    else
      maxMatrix[maxI][maxJ] = std::numeric_limits<double>::max();
  }
  {
    if(reverseSad)
      sadMatrix[sadJ][sadI] = std::numeric_limits<double>::max();
    else
      sadMatrix[sadI][sadJ] = std::numeric_limits<double>::max();
  }
}

double ttk::BottleneckDistance::computeGeometricalRange(
  const ttk::DiagramType &CTDiagram1,
  const ttk::DiagramType &CTDiagram2,
  const int d1Size,
  const int d2Size) const {

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

  return std::sqrt(Geometry::pow(maxX - minX, 2) + Geometry::pow(maxY - minY, 2)
                   + Geometry::pow(maxZ - minZ, 2));
}

double ttk::BottleneckDistance::computeMinimumRelevantPersistence(
  const ttk::DiagramType &CTDiagram1,
  const ttk::DiagramType &CTDiagram2,
  const int d1Size,
  const int d2Size) const {

  double sp = zeroThreshold_;
  double s = sp > 0.0 && sp < 100.0 ? sp / 100.0 : 0;

  std::vector<double> toSort;
  for(int i = 0; i < d1Size; ++i) {
    const auto &t = CTDiagram1[i];
    const auto persistence = std::abs(t.persistence);
    toSort.push_back(persistence);
  }
  for(int i = 0; i < d2Size; ++i) {
    const auto &t = CTDiagram2[i];
    const auto persistence = std::abs(t.persistence);
    toSort.push_back(persistence);
  }
  sort(toSort.begin(), toSort.end());

  double minVal = toSort.at(0);
  double maxVal = toSort.at(toSort.size() - 1);
  s *= (maxVal - minVal);

  // double epsilon = 0.0000001;
  // int largeSize = 2000;
  // dataType zeroThresh = (dataType) epsilon;
  // if (d1Size + d2Size > largeSize + 1) {
  //   zeroThresh = toSort.at(d1Size + d2Size - largeSize);
  //   if (toSort.at(d1Size + d2Size - (largeSize+1)) == zeroThresh)
  //     zeroThresh += (dataType) epsilon;
  // }
  // if (zeroThresh < epsilon) zeroThresh = epsilon;

  return s;
}

void ttk::BottleneckDistance::computeMinMaxSaddleNumberAndMapping(
  const ttk::DiagramType &CTDiagram,
  int dSize,
  int &nbMin,
  int &nbMax,
  int &nbSaddle,
  std::vector<int> &minMap,
  std::vector<int> &maxMap,
  std::vector<int> &sadMap,
  const double zeroThresh) const {

  for(int i = 0; i < dSize; ++i) {
    const auto &t = CTDiagram[i];
    const auto nt1 = t.birth.type;
    const auto nt2 = t.death.type;
    const auto dt = t.persistence;
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

void solvePWasserstein(const int nbRow,
                       const int nbCol,
                       std::vector<std::vector<double>> &matrix,
                       std::vector<MatchingType> &matchings,
                       ttk::AssignmentMunkres<double> &solver) {

  solver.setInput(matrix);
  solver.run(matchings);
  solver.clearMatrix();
}

void solveInfinityWasserstein(const int nbRow,
                              const int nbCol,
                              const int nbRowToCut,
                              const int nbColToCut,
                              std::vector<std::vector<double>> &matrix,
                              std::vector<MatchingType> &matchings,
                              ttk::GabowTarjan &solver) {

  // Copy input matrix.
  auto bottleneckMatrix = matrix;

  // Solve.
  solver.setInput(nbRow, nbCol, &bottleneckMatrix);
  solver.run(matchings);
  solver.clear();
}

double ttk::BottleneckDistance::buildMappings(
  const std::vector<MatchingType> &inputMatchings,
  const bool transposeGlobal,
  const bool transposeLocal,
  std::vector<MatchingType> &outputMatchings,
  const std::vector<int> &m1,
  const std::vector<int> &m2,
  int wasserstein) const {

  // Input map permutation (so as to ignore transposition later on)
  const std::vector<int> map1 = transposeLocal ? m2 : m1;
  const std::vector<int> map2 = transposeLocal ? m1 : m2;

  double addedPersistence = 0;
  for(int i = 0, s = (int)inputMatchings.size(); i < s; ++i) {
    const auto &t = inputMatchings.at(i);
    const auto val = std::abs(std::get<2>(t));

    int p1 = std::get<0>(t);
    int p2 = std::get<1>(t);

    if(p1 >= (int)map1.size() || p1 < 0 || p2 >= (int)map2.size() || p2 < 0) {
      addedPersistence = (wasserstein > 0 ? addedPersistence + val
                                          : std::max(val, addedPersistence));
    } else {
      int point1 = map1.at((unsigned long)p1);
      int point2 = map2.at((unsigned long)p2);
      bool doTranspose = transposeGlobal ^ transposeLocal;

      const auto &newT = doTranspose ? std::make_tuple(point2, point1, val)
                                     : std::make_tuple(point1, point2, val);

      outputMatchings.push_back(newT);
    }
  }

  return addedPersistence;
}

int ttk::BottleneckDistance::computeBottleneck(
  const ttk::DiagramType &d1,
  const ttk::DiagramType &d2,
  std::vector<MatchingType> &matchings,
  const bool usePersistenceMetric) {

  auto d1Size = (int)d1.size();
  auto d2Size = (int)d2.size();

  bool transposeOriginal = d1Size > d2Size;
  const ttk::DiagramType &CTDiagram1 = transposeOriginal ? d2 : d1;
  const ttk::DiagramType &CTDiagram2 = transposeOriginal ? d1 : d2;
  if(transposeOriginal) {
    int temp = d1Size;
    d1Size = d2Size;
    d2Size = temp;
  }

  if(transposeOriginal) {
    this->printMsg("The first persistence diagram is larger than the second.");
    this->printMsg("Solving the transposed problem.");
  }

  // Check user parameters.
  const int wasserstein = (wasserstein_ == "inf") ? -1 : stoi(wasserstein_);
  if(wasserstein < 0 && wasserstein != -1)
    return -4;

  // Needed to limit computation time.
  const auto zeroThresh = this->computeMinimumRelevantPersistence(
    CTDiagram1, CTDiagram2, d1Size, d2Size);

  // Initialize solvers.
  std::vector<MatchingType> minMatchings;
  std::vector<MatchingType> maxMatchings;
  std::vector<MatchingType> sadMatchings;

  // Initialize cost matrices.
  int nbRowMin = 0, nbColMin = 0;
  int maxRowColMin = 0, minRowColMin = 0;
  int nbRowMax = 0, nbColMax = 0;
  int maxRowColMax = 0, minRowColMax = 0;
  int nbRowSad = 0, nbColSad = 0;
  int maxRowColSad = 0, minRowColSad = 0;

  // Remap for matchings.
  std::vector<int> minMap1;
  std::vector<int> minMap2;
  std::vector<int> maxMap1;
  std::vector<int> maxMap2;
  std::vector<int> sadMap1;
  std::vector<int> sadMap2;

  this->computeMinMaxSaddleNumberAndMapping(CTDiagram1, d1Size, nbRowMin,
                                            nbRowMax, nbRowSad, minMap1,
                                            maxMap1, sadMap1, zeroThresh);
  this->computeMinMaxSaddleNumberAndMapping(CTDiagram2, d2Size, nbColMin,
                                            nbColMax, nbColSad, minMap2,
                                            maxMap2, sadMap2, zeroThresh);

  // Automatically transpose if nb rows > nb cols
  maxRowColMin = std::max(nbRowMin + 1, nbColMin + 1);
  maxRowColMax = std::max(nbRowMax + 1, nbColMax + 1);
  maxRowColSad = std::max(nbRowSad + 1, nbColSad + 1);

  minRowColMin = std::min(nbRowMin + 1, nbColMin + 1);
  minRowColMax = std::min(nbRowMax + 1, nbColMax + 1);
  minRowColSad = std::min(nbRowSad + 1, nbColSad + 1);

  std::vector<std::vector<double>> minMatrix(
    (unsigned long)minRowColMin, std::vector<double>(maxRowColMin));
  std::vector<std::vector<double>> maxMatrix(
    (unsigned long)minRowColMax, std::vector<double>(maxRowColMax));
  std::vector<std::vector<double>> sadMatrix(
    (unsigned long)minRowColSad, std::vector<double>(maxRowColSad));

  double px = px_;
  double py = py_;
  double pz = pz_;
  double pe = pe_;
  double ps = ps_;

  const auto distanceFunction = [wasserstein, px, py, pz, pe, ps](
                                  const ttk::PersistencePair &a,
                                  const ttk::PersistencePair &b) -> double {
    const auto ta1 = a.birth.type;
    const auto ta2 = a.death.type;
    const int w = wasserstein > 1 ? wasserstein : 1; // L_inf not managed.

    // We don't match critical points of different index.
    // This must be ensured before calling the distance function.
    // const auto tb1 = get<1>(b);
    // const auto tb2 = get<3>(b);
    bool isMin1 = ta1 == CriticalType::Local_minimum;
    bool isMax1 = ta2 == CriticalType::Local_maximum;
    // bool isBoth = isMin1 && isMax1;

    const auto rX = a.birth.sfValue;
    const auto rY = a.death.sfValue;
    const auto cX = b.birth.sfValue;
    const auto cY = b.death.sfValue;
    const auto x
      = ((isMin1 && !isMax1) ? pe : ps) * Geometry::pow(std::abs(rX - cX), w);
    const auto y = (isMax1 ? pe : ps) * Geometry::pow(std::abs(rY - cY), w);
    double geoDistance
      = isMax1 ? (
          px * Geometry::pow(abs(a.death.coords[0] - b.death.coords[0]), w)
          + py * Geometry::pow(abs(a.death.coords[1] - b.death.coords[1]), w)
          + pz * Geometry::pow(abs(a.death.coords[2] - b.death.coords[2]), w))
        : isMin1 ? (
            px * Geometry::pow(abs(a.birth.coords[0] - b.birth.coords[0]), w)
            + py * Geometry::pow(abs(a.birth.coords[1] - b.birth.coords[1]), w)
            + pz * Geometry::pow(abs(a.birth.coords[2] - b.birth.coords[2]), w))
                 : (px
                      * Geometry::pow(
                        abs(a.birth.coords[0] + a.death.coords[0]) / 2
                          - abs(b.birth.coords[0] + b.death.coords[0]) / 2,
                        w)
                    + py
                        * Geometry::pow(
                          abs(a.birth.coords[1] + a.death.coords[1]) / 2
                            - abs(b.birth.coords[1] + b.death.coords[1]) / 2,
                          w)
                    + pz
                        * Geometry::pow(
                          abs(a.birth.coords[2] + a.death.coords[2]) / 2
                            - abs(b.birth.coords[2] + b.death.coords[2]) / 2,
                          w));

    double persDistance = x + y;
    double val = persDistance + geoDistance;
    val = Geometry::pow(val, 1.0 / w);
    return val;
  };

  const auto diagonalDistanceFunction =
    [wasserstein, px, py, pz, ps, pe](const ttk::PersistencePair a) -> double {
    const auto ta1 = a.birth.type;
    const auto ta2 = a.death.type;
    const int w = wasserstein > 1 ? wasserstein : 1;
    bool isMin1 = ta1 == CriticalType::Local_minimum;
    bool isMax1 = ta2 == CriticalType::Local_maximum;

    const auto rX = a.birth.sfValue;
    const auto rY = a.death.sfValue;
    double x1 = a.birth.coords[0];
    double y1 = a.birth.coords[1];
    double z1 = a.birth.coords[2];
    double x2 = a.death.coords[0];
    double y2 = a.death.coords[1];
    double z2 = a.death.coords[2];

    double infDistance
      = (isMin1 || isMax1 ? pe : ps) * Geometry::pow(std::abs(rX - rY), w);
    double geoDistance = (px * Geometry::pow(abs(x2 - x1), w)
                          + py * Geometry::pow(abs(y2 - y1), w)
                          + pz * Geometry::pow(abs(z2 - z1), w));
    double val = infDistance + geoDistance;
    return Geometry::pow(val, 1.0 / w);
  };

  const bool transposeMin = nbRowMin > nbColMin;
  const bool transposeMax = nbRowMax > nbColMax;
  const bool transposeSad = nbRowSad > nbColSad;

  Timer t;

  this->buildCostMatrices(
    CTDiagram1, CTDiagram2, d1Size, d2Size, distanceFunction,
    diagonalDistanceFunction, zeroThresh, minMatrix, maxMatrix, sadMatrix,
    transposeMin, transposeMax, transposeSad, wasserstein);

  if(wasserstein > 0) {

    if(nbRowMin > 0 && nbColMin > 0) {
      AssignmentMunkres<double> solverMin;
      this->printMsg("Affecting minima...");
      solvePWasserstein(
        minRowColMin, maxRowColMin, minMatrix, minMatchings, solverMin);
    }

    if(nbRowMax > 0 && nbColMax > 0) {
      AssignmentMunkres<double> solverMax;
      this->printMsg("Affecting maxima...");
      solvePWasserstein(
        minRowColMax, maxRowColMax, maxMatrix, maxMatchings, solverMax);
    }

    if(nbRowSad > 0 && nbColSad > 0) {
      AssignmentMunkres<double> solverSad;
      this->printMsg("Affecting saddles...");
      solvePWasserstein(
        minRowColSad, maxRowColSad, sadMatrix, sadMatchings, solverSad);
    }

  } else {

    // Launch solving for minima.
    if(nbRowMin > 0 && nbColMin > 0) {
      GabowTarjan solverMin;
      this->printMsg("Affecting minima...");
      solveInfinityWasserstein(minRowColMin, maxRowColMin, nbRowMin, nbColMin,
                               minMatrix, minMatchings, solverMin);
    }

    // Launch solving for maxima.
    if(nbRowMax > 0 && nbColMax > 0) {
      GabowTarjan solverMax;
      this->printMsg("Affecting maxima...");
      solveInfinityWasserstein(minRowColMax, maxRowColMax, nbRowMax, nbColMax,
                               maxMatrix, maxMatchings, solverMax);
    }

    // Launch solving for saddles.
    if(nbRowSad > 0 && nbColSad > 0) {
      GabowTarjan solverSad;
      this->printMsg("Affecting saddles...");
      solveInfinityWasserstein(minRowColSad, maxRowColSad, nbRowSad, nbColSad,
                               sadMatrix, sadMatchings, solverSad);
    }
  }

  this->printMsg("TTK CORE DONE", 1, t.getElapsedTime());

  // Rebuild mappings.
  // Begin cost computation for unpaired vertices.
  // std::cout << "Min" << std::endl;
  const auto addedMinPersistence
    = this->buildMappings(minMatchings, transposeOriginal, transposeMin,
                          matchings, minMap1, minMap2, wasserstein);

  // std::cout << "Max" << std::endl;
  const auto addedMaxPersistence
    = this->buildMappings(maxMatchings, transposeOriginal, transposeMax,
                          matchings, maxMap1, maxMap2, wasserstein);

  // std::cout << "Sad" << std::endl;
  const auto addedSadPersistence
    = this->buildMappings(sadMatchings, transposeOriginal, transposeSad,
                          matchings, sadMap1, sadMap2, wasserstein);

  // TODO [HIGH] do that for embeddings
  // Recompute matching weights for user-friendly distance.
  double d = 0;
  std::vector<bool> paired1(d1Size);
  std::vector<bool> paired2(d2Size);
  for(int b = 0; b < d1Size; ++b)
    paired1[b] = false;
  for(int b = 0; b < d2Size; ++b)
    paired2[b] = false;

  int numberOfMismatches = 0;
  for(int m = 0, ms = (int)matchings.size(); m < ms; ++m) {
    MatchingType mt = matchings[m];
    int i = transposeOriginal ? std::get<1>(mt) : std::get<0>(mt);
    int j = transposeOriginal ? std::get<0>(mt) : std::get<1>(mt);
    // dataType val = std::get<2>(t);

    const auto &t1 = CTDiagram1[i];
    const auto &t2 = CTDiagram2[j];
    // dataType rX = std::get<6>(t1); dataType rY = std::get<10>(t1);
    // dataType cX = std::get<6>(t2); dataType cY = std::get<10>(t2);
    // dataType x = rX - cX; dataType y = rY - cY;
    paired1[i] = true;
    paired2[j] = true;
    // dataType lInf = std::max(abs<dataType>(x), abs<dataType>(y));

    // if (((wasserstein < 0 && lInf != val) || (wasserstein > 0 && pow(lInf,
    // wasserstein) != val)))
    //++numberOfMismatches;

    auto partialDistance = distanceFunction(t1, t2);
    // wasserstein > 0 ? pow(lInf, wasserstein) : std::max(d, lInf);

    if(wasserstein > 0)
      d += partialDistance;
    else
      d = partialDistance;
  }

  if(numberOfMismatches > 0) {
    this->printWrn("Distance mismatch when rebuilding "
                   + std::to_string(numberOfMismatches) + " matchings");
  }

  auto affectationD = d;
  d = wasserstein > 0
        ? Geometry::pow(
          d + addedMaxPersistence + addedMinPersistence + addedSadPersistence,
          (1.0 / (double)wasserstein))
        : std::max(
          d, std::max(addedMaxPersistence,
                      std::max(addedMinPersistence, addedSadPersistence)));

  {
    std::stringstream msg;
    this->printMsg("Computed distance:");
    this->printMsg("diagMax(" + std::to_string(addedMaxPersistence)
                   + "), diagMin(" + std::to_string(addedMinPersistence)
                   + "), diagSad(" + std::to_string(addedSadPersistence) + ")");
    this->printMsg("affAll(" + std::to_string(affectationD) + "), res("
                   + std::to_string(d) + ")");
  }

  distance_ = d;
  return 0;
}
