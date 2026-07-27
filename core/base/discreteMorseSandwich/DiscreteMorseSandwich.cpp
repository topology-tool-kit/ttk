#include <DiscreteMorseSandwich.h>

ttk::DiscreteMorseSandwich::DiscreteMorseSandwich() {
  this->setDebugMsgPrefix("DiscreteMorseSandwich");
}

void ttk::DiscreteMorseSandwich::tripletsToPersistencePairs(
  std::vector<PersistencePair> &pairs,
  std::vector<bool> &pairedExtrema,
  std::vector<bool> &pairedSaddles,
  std::vector<SimplexId> &reps,
  std::vector<tripletType> &triplets,
  const SimplexId *const saddlesOrder,
  const SimplexId *const extremaOrder,
  const SimplexId pairDim) const {

  // comparison functions
  const auto cmpSadMax
    = [=](const tripletType &t0, const tripletType &t1) -> bool {
    const auto s0 = t0[0];
    const auto s1 = t1[0];
    const auto m0 = t0[2];
    const auto m1 = t1[2];

#ifdef _LIBCPP_VERSION
    // libc++'s std::sort compares an entry to itself
    if(&t0 == &t1) {
      return true;
    }
#endif // _LIBCPP_VERSION

    if(s0 != s1)
      return saddlesOrder[s0] > saddlesOrder[s1];
    else
      return extremaOrder[m0] < extremaOrder[m1];
  };

  const auto cmpSadMin
    = [=](const tripletType &t0, const tripletType &t1) -> bool {
    const auto s0 = t0[0];
    const auto s1 = t1[0];
    const auto m0 = t0[2];
    const auto m1 = t1[2];
    if(s0 != s1)
      return saddlesOrder[s0] < saddlesOrder[s1];
    else
      return extremaOrder[m0] > extremaOrder[m1];
  };

  // sort triplets
  if(pairDim == 0) {
    TTK_PSORT(this->threadNumber_, triplets.begin(), triplets.end(), cmpSadMin);
  } else {
    // saddle-saddle pairs from 1-saddles to 2-saddles
    TTK_PSORT(this->threadNumber_, triplets.begin(), triplets.end(), cmpSadMax);
  }

  // get representative of current extremum
  const auto getRep = [&reps](SimplexId v) -> SimplexId {
    auto r = reps[v];
    while(r != v) {
      v = r;
      r = reps[v];
    }
    return r;
  };

  const bool increasing = (pairDim > 0);

  const auto addPair = [&pairs, &pairedExtrema, &pairedSaddles, increasing,
                        pairDim](const SimplexId sad, const SimplexId extr) {
    if(increasing) {
      pairs.emplace_back(sad, extr, pairDim);
    } else {
      pairs.emplace_back(extr, sad, pairDim);
    }
    pairedSaddles[sad] = true;
    pairedExtrema[extr] = true;
  };

  for(const auto &t : triplets) {
    const auto sv = t[0];

    auto r1 = getRep(t[1]);

    if(t[2] < 0) {
      // deal with "shadow" triplets (a 2-saddle with only one
      // ascending 1-separatrix leading to an unique maximum)
      if(!pairedExtrema[r1] && !pairedSaddles[sv]) {
        // when considering the boundary, the "-1" of the triplets
        // indicate a virtual maximum of infinite persistence on the
        // boundary component. a pair is created with the other
        // maximum
        addPair(sv, r1);
      }

      continue;
    }

    auto r2 = getRep(t[2]);

    if(r1 != r2) {
      if(((extremaOrder[r1] > extremaOrder[r2]) == increasing
          || pairedExtrema[r1])
         && !pairedExtrema[r2]) {
        std::swap(r1, r2);
      }
      if(!pairedExtrema[r1]) {
        addPair(sv, r1);
        reps[t[1]] = r2;
        reps[r1] = r2;
      }
    }
  }
}

void ttk::DiscreteMorseSandwich::displayStats(
  const std::vector<PersistencePair> &pairs,
  const std::array<std::vector<SimplexId>, 4> &criticalCellsByDim,
  const std::vector<bool> &pairedMinima,
  const std::vector<bool> &paired1Saddles,
  const std::vector<bool> &paired2Saddles,
  const std::vector<bool> &pairedMaxima) const {

  const auto dim = this->dg_.getDimensionality();

  // display number of pairs per pair type
  std::vector<std::vector<std::string>> rows{
    {" #Min-saddle pairs",
     std::to_string(
       std::count_if(pairs.begin(), pairs.end(),
                     [](const PersistencePair &a) { return a.type == 0; }))},
    {" #Saddle-saddle pairs",
     std::to_string(dim == 3 ? std::count_if(
                      pairs.begin(), pairs.end(),
                      [](const PersistencePair &a) { return a.type == 1; })
                             : 0)},
    {" #Saddle-max pairs",
     std::to_string(std::count_if(
       pairs.begin(), pairs.end(),
       [dim](const PersistencePair &a) { return a.type == dim - 1; }))},
  };

  // display number of critical cells (paired and unpaired)
  std::vector<size_t> nCritCells(dim + 1);
  std::vector<size_t> nNonPairedCritCells(dim + 1);

  for(int i = 0; i < dim + 1; ++i) {
    nCritCells[i] = criticalCellsByDim[i].size();
    size_t nNonPaired{};
    for(size_t j = 0; j < criticalCellsByDim[i].size(); ++j) {
      const auto cell = criticalCellsByDim[i][j];
      if((i == 0 && !pairedMinima[cell]) || (i == 1 && !paired1Saddles[cell])
         || (i == 2 && dim == 3 && !paired2Saddles[cell])
         || (i == dim && !pairedMaxima[cell])) {
        nNonPaired++;
      }
    }
    nNonPairedCritCells[i] = nNonPaired;
  }

  std::vector<std::string> critCellsLabels{"Minima"};
  if(dim >= 2) {
    critCellsLabels.emplace_back("1-saddles");
  }
  if(dim >= 3) {
    critCellsLabels.emplace_back("2-saddles");
  }
  critCellsLabels.emplace_back("Maxima");

  for(int i = 0; i < dim + 1; ++i) {
    const std::string unpaired{nNonPairedCritCells[i] == 0
                                 ? " (all paired)"
                                 : " (" + std::to_string(nNonPairedCritCells[i])
                                     + " unpaired)"};

    rows.emplace_back(std::vector<std::string>{
      " #" + critCellsLabels[i], std::to_string(nCritCells[i]) + unpaired});
  }
  this->printMsg(rows, debug::Priority::DETAIL);
}
