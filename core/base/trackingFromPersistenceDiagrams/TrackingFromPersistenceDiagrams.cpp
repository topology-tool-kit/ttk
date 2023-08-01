#include <TrackingFromPersistenceDiagrams.h>

ttk::TrackingFromPersistenceDiagrams::TrackingFromPersistenceDiagrams() {
  this->setDebugMsgPrefix("TrackingFromPersistenceDiagrams");
}

int ttk::TrackingFromPersistenceDiagrams::execute() {
  ttk::Timer t;

  // Check the consistency of the variables
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfInputs_)
    return -1;
  if(!inputData_)
    return -3;

  for(int i = 0; i < numberOfInputs_; i++) {
    if(!inputData_[i])
      return -4;
  }
#endif

  this->printMsg("Complete", 1, t.getElapsedTime(), threadNumber_);

  return 0;
}

int ttk::TrackingFromPersistenceDiagrams::performSingleMatching(
  int i,
  std::vector<ttk::DiagramType> &inputPersistenceDiagrams,
  std::vector<std::vector<MatchingType>> &outputMatchings,
  const std::string &algorithm,
  const std::string &wasserstein,
  double tolerance,
  double px,
  double py,
  double pz,
  double ps,
  double pe) {
  ttk::BottleneckDistance bottleneckDistance_;
  // bottleneckDistance_.setWrapper(this);
  bottleneckDistance_.setPersistencePercentThreshold(tolerance);
  bottleneckDistance_.setPX(px);
  bottleneckDistance_.setPY(py);
  bottleneckDistance_.setPZ(pz);
  bottleneckDistance_.setPS(ps);
  bottleneckDistance_.setPE(pe);
  bottleneckDistance_.setAlgorithm(algorithm);
  bottleneckDistance_.setWasserstein(wasserstein);

  bottleneckDistance_.execute(inputPersistenceDiagrams[i],
                              inputPersistenceDiagrams[i + 1],
                              outputMatchings[i]);

  return 0;
}

int ttk::TrackingFromPersistenceDiagrams::performMatchings(
  int numInputs,
  std::vector<ttk::DiagramType> &inputPersistenceDiagrams,
  std::vector<std::vector<MatchingType>> &outputMatchings,
  const std::string &algorithm,
  const std::string &wasserstein,
  double tolerance,
  double px,
  double py,
  double pz,
  double ps,
  double pe) {

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < numInputs - 1; ++i) {
    performSingleMatching(
      i, inputPersistenceDiagrams, outputMatchings,
      algorithm, // Not from paraview, from enclosing tracking plugin
      wasserstein, tolerance, px, py, pz, ps, pe // Coefficients
    );
  }

  // Never from PV,
  // Alaways activate output matchings.
  return 0;
}

int ttk::TrackingFromPersistenceDiagrams::performTracking(
  std::vector<ttk::DiagramType> &allDiagrams,
  std::vector<std::vector<MatchingType>> &allMatchings,
  std::vector<trackingTuple> &trackings) {
  auto numPersistenceDiagramsInput = (int)allDiagrams.size();

  for(int in = 1; in < numPersistenceDiagramsInput - 1; ++in) {
    std::vector<MatchingType> matchings1 = allMatchings[in - 1];
    std::vector<MatchingType> matchings2 = allMatchings[in];

    auto matchingsSize1 = (int)matchings1.size();
    auto matchingsSize2 = (int)matchings2.size();
    int const endIndex = numPersistenceDiagramsInput - 2;

    for(int i = 0; i < matchingsSize1; ++i) {
      const auto m1ai0 = std::get<0>(matchings1[i]);
      const auto m1ai1 = std::get<1>(matchings1[i]);

      for(int j = 0; j < matchingsSize2; ++j) {
        const auto m2aj0 = std::get<0>(matchings2[j]);
        const auto m2aj1 = std::get<1>(matchings2[j]);

        if(m1ai1 != m2aj0)
          continue;

        // Detect in trackings and push.
        bool found = false;
        for(trackingTuple &tt : trackings) {
          int const chainStart = std::get<0>(tt);
          int const chainEnd = std::get<1>(tt);
          std::vector<SimplexId> &chain = std::get<2>(tt);

          if(chainEnd == -1) {
            auto chainSize = (int)chain.size();
            if(chainSize == 0) {
              // Should not happen
              this->printErr("Brain error.");
            } else if(chainStart + chainSize == in
                      && chain.at((unsigned long)chainSize - 1) == m1ai0) {
              found = true;
              chain.push_back(m1ai1);
              int const numEnd = in == endIndex ? endIndex : -1;
              if(in == endIndex) {
                chain.push_back(m2aj1);
                std::get<1>(tt) = numEnd;
              }
              std::get<2>(tt) = chain;
            }
          }
          tt = std::make_tuple(chainStart, chainEnd, chain);
        }

        if(!found) {
          std::vector<SimplexId> chain;
          chain.push_back(m1ai0);
          chain.push_back(m1ai1);
          if(in == endIndex) {
            chain.push_back(m2aj1);
          }
          int const numEnd = in == endIndex ? endIndex : -1;
          trackingTuple const tt = std::make_tuple(in - 1, numEnd, chain);
          trackings.push_back(tt);
        }
        // Create new.
      }
    }

    // End non-matched chains.
    for(trackingTuple &tt : trackings) {
      int const chainStart = std::get<0>(tt);
      int const chainEnd = std::get<1>(tt);
      if(chainEnd == -1) {
        std::vector<SimplexId> const &chain = std::get<2>(tt);
        auto chainSize = (int)chain.size();
        if(chainStart + chainSize - 1 < in)
          std::get<1>(tt) = in - 1;
      }
    }
  }

  // Post-processing
  std::sort(trackings.begin(), trackings.end(),
            [](const trackingTuple &a, const trackingTuple &b) -> bool {
              return std::get<0>(a) < std::get<0>(b);
            });

  return 0;
}

int ttk::TrackingFromPersistenceDiagrams::performPostProcess(
  const std::vector<ttk::DiagramType> &allDiagrams,
  const std::vector<trackingTuple> &trackings,
  std::vector<std::set<int>> &trackingTupleToMerged,
  const double postProcThresh) {

  const int numPersistenceDiagramsInput = allDiagrams.size();

  // Merge close connected components with threshold.
  for(size_t k = 0; k < trackings.size(); ++k) {
    const auto &tk = trackings[k];
    int const startK = std::get<0>(tk);
    int endK = std::get<1>(tk);
    if(endK < 0)
      endK = numPersistenceDiagramsInput - 1;
    const std::vector<SimplexId> &chainK = std::get<2>(tk);
    const ttk::DiagramType &diagramStartK = allDiagrams[startK];
    const ttk::DiagramType &diagramEndK = allDiagrams[endK];

    const auto n1 = chainK.front();
    const auto n2 = chainK.back();
    const auto &tuple1 = diagramStartK[n1];
    const auto &tuple2 = diagramEndK[n2];

    double x1, y1, z1, x2, y2, z2;

    const auto point1Type1 = tuple1.birth.type;
    const auto point1Type2 = tuple1.death.type;
    bool const t11Min = point1Type1 == CriticalType::Local_minimum;
    bool const t11Max = point1Type1 == CriticalType::Local_maximum;
    bool const t12Min = point1Type2 == CriticalType::Local_minimum;
    bool const t12Max = point1Type2 == CriticalType::Local_maximum;
    // bool bothEx1 = t11Ex && t12Ex;
    bool const t1Max = t11Max || t12Max;
    bool const t1Min = !t1Max && (t11Min || t12Min);

    x1 = t1Max ? tuple1.death.coords[0] : t1Min ? tuple1.birth.coords[0] : 0;
    y1 = t1Max ? tuple1.death.coords[1] : t1Min ? tuple1.birth.coords[1] : 0;
    z1 = t1Max ? tuple1.death.coords[2] : t1Min ? tuple1.birth.coords[2] : 0;

    const auto point2Type1 = tuple2.birth.type;
    const auto point2Type2 = tuple2.death.type;
    bool const t21Min = point2Type1 == CriticalType::Local_minimum;
    bool const t21Max = point2Type1 == CriticalType::Local_maximum;
    bool const t22Min = point2Type2 == CriticalType::Local_minimum;
    bool const t22Max = point2Type2 == CriticalType::Local_maximum;
    // bool bothEx2 = t21Ex && t22Ex;
    bool const t2Max = t21Max || t22Max;
    bool const t2Min = !t2Max && (t21Min || t22Min);

    // if (bothEx2) {
    x2 = t2Max ? tuple2.death.coords[0] : t2Min ? tuple2.birth.coords[0] : 0;
    y2 = t2Max ? tuple2.death.coords[1] : t2Min ? tuple2.birth.coords[1] : 0;
    z2 = t2Max ? tuple2.death.coords[2] : t2Min ? tuple2.birth.coords[2] : 0;
    // }
    // if (!bothEx1 && !bothEx2)
    //  continue;

    // Saddle-saddle matching not supported.
    if(!t1Min && !t2Min && !t1Max && !t2Max)
      continue;

    // Check every other tracking trajectory.
    for(size_t m = k + 1; m < trackings.size(); ++m) {
      const auto &tm = trackings[m];
      int const startM = std::get<0>(tm);
      int const endM = std::get<1>(tm);
      const std::vector<SimplexId> &chainM = std::get<2>(tm);
      if((endK > 0 && startM > endK) || (endM > 0 && startK > endM))
        continue;

      for(int c = 0; c < (int)chainM.size(); ++c) {
        bool const doMatch1 = startM + c == startK;
        bool const doMatch2 = startM + c == endK;

        // if (startM + c != startK && startM + c != endK) continue;
        if(!doMatch1 && !doMatch2)
          continue;

        /// Check proximity.
        const auto n3 = chainM[c];
        const ttk::DiagramType &diagramM = allDiagrams[startM + c];
        const auto &tuple3 = diagramM[n3];
        double x3, y3, z3;
        const auto point3Type1 = tuple3.birth.type;
        const auto point3Type2 = tuple3.death.type;
        bool const t31Min = point3Type1 == CriticalType::Local_minimum;
        bool const t31Max = point3Type1 == CriticalType::Local_maximum;
        bool const t32Min = point3Type2 == CriticalType::Local_minimum;
        bool const t32Max = point3Type2 == CriticalType::Local_maximum;
        // bool bothEx3 = t31Ex && t32Ex;
        // if (!bothEx3)
        //  continue;
        bool const t3Max = t31Max || t32Max;
        bool const t3Min = !t3Max && (t31Min || t32Min);

        x3 = t3Max   ? tuple3.death.coords[0]
             : t3Min ? tuple3.birth.coords[0]
                     : 0;
        y3 = t3Max   ? tuple3.death.coords[1]
             : t3Min ? tuple3.birth.coords[1]
                     : 0;
        z3 = t3Max   ? tuple3.death.coords[2]
             : t3Min ? tuple3.birth.coords[2]
                     : 0;

        double dist = 0;
        bool hasMatched = false;
        const auto square = [](const double a) { return a * a; };
        if(doMatch1 && ((t3Max && t1Max) || (t3Min && t1Min))) {
          double const dist13
            = std::sqrt(square(x1 - x3) + square(y1 - y3) + square(z1 - z3));
          dist = dist13;
          if(dist13 >= postProcThresh)
            continue;
          hasMatched = true;
        }

        if(doMatch2 && ((t3Max && t2Max) || (t3Min && t2Min))) {
          double const dist23
            = std::sqrt(square(x2 - x3) + square(y2 - y3) + square(z2 - z3));
          dist = dist23;
          if(dist23 >= postProcThresh)
            continue;
          hasMatched = true;
        }

        if(!hasMatched)
          continue;

        /// Merge!
        std::stringstream msg;
        msg << "Merged " << m << " with " << k << ": d = " << dist << ".";
        printMsg(msg.str());

        // Get every other tracking trajectory.
        std::set<int> &mergedM = trackingTupleToMerged[m];
        // std::set<int> mergedK = trackingTupleToMerged[k];

        // Push for others to merge.
        // for (auto& i : mergedM) mergedK.insert(i);
        // for (auto& i : mergedK) mergedM.insert(i);
        // mergedK.insert(m);
        mergedM.insert(k);
        break;
      }
    }
  }

  return 0;
}
