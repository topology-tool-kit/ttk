#include <algorithm>

#include <PersistenceDiagramDistanceMatrix.h>

using namespace ttk;

std::vector<std::vector<double>> PersistenceDiagramDistanceMatrix::execute(
  std::vector<std::vector<DiagramTuple>> &intermediateDiagrams) const {

  Timer tm{};

  const auto nInputs = intermediateDiagrams.size();

  if(do_min_ && do_sad_ && do_max_) {
    this->printMsg("Processing all critical pairs types");
  } else if(do_min_) {
    this->printMsg("Processing only MIN-SAD pairs");
  } else if(do_sad_) {
    this->printMsg("Processing only SAD-SAD pairs");
  } else if(do_max_) {
    this->printMsg("Processing only SAD-MAX pairs");
  }

  std::vector<std::vector<DiagramTuple>> inputDiagramsMin(nInputs);
  std::vector<std::vector<DiagramTuple>> inputDiagramsSad(nInputs);
  std::vector<std::vector<DiagramTuple>> inputDiagramsMax(nInputs);

  std::vector<BidderDiagram<double>> bidder_diagrams_min{};
  std::vector<BidderDiagram<double>> bidder_diagrams_sad{};
  std::vector<BidderDiagram<double>> bidder_diagrams_max{};
  std::vector<BidderDiagram<double>> current_bidder_diagrams_min{};
  std::vector<BidderDiagram<double>> current_bidder_diagrams_sad{};
  std::vector<BidderDiagram<double>> current_bidder_diagrams_max{};

  // Store the persistence of the global min-max pair
  std::vector<double> maxDiagPersistence(nInputs);

  // Create diagrams for min, saddle and max persistence pairs
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs; i++) {
    std::vector<DiagramTuple> &CTDiagram = intermediateDiagrams[i];

    for(size_t j = 0; j < CTDiagram.size(); ++j) {
      const DiagramTuple &t = CTDiagram[j];
      const ttk::CriticalType nt1 = std::get<1>(t);
      const ttk::CriticalType nt2 = std::get<3>(t);
      const double pers = std::get<4>(t);
      maxDiagPersistence[i] = std::max(pers, maxDiagPersistence[i]);

      if(pers > 0) {
        if(nt1 == CriticalType::Local_minimum
           && nt2 == CriticalType::Local_maximum) {
          inputDiagramsMax[i].emplace_back(t);
        } else {
          if(nt1 == CriticalType::Local_maximum
             || nt2 == CriticalType::Local_maximum) {
            inputDiagramsMax[i].emplace_back(t);
          }
          if(nt1 == CriticalType::Local_minimum
             || nt2 == CriticalType::Local_minimum) {
            inputDiagramsMin[i].emplace_back(t);
          }
          if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
             || (nt1 == CriticalType::Saddle2
                 && nt2 == CriticalType::Saddle1)) {
            inputDiagramsSad[i].emplace_back(t);
          }
        }
      }
    }
  }

  if(this->do_min_) {
    setBidderDiagrams(nInputs, inputDiagramsMin, bidder_diagrams_min);
  }
  if(this->do_sad_) {
    setBidderDiagrams(nInputs, inputDiagramsSad, bidder_diagrams_sad);
  }
  if(this->do_max_) {
    setBidderDiagrams(nInputs, inputDiagramsMax, bidder_diagrams_max);
  }

  switch(this->Constraint) {
    case ConstraintType::FULL_DIAGRAMS:
      this->printMsg("Using all diagram pairs");
      break;
    case ConstraintType::NUMBER_PAIRS:
      this->printMsg("Using the " + std::to_string(this->MaxNumberOfPairs)
                     + " most persistent pairs");
      break;
    case ConstraintType::ABSOLUTE_PERSISTENCE: {
      std::stringstream pers{};
      pers << std::fixed << std::setprecision(2) << this->MinPersistence;
      this->printMsg("Using diagram pairs above a persistence threshold of "
                     + pers.str());
    } break;
    case ConstraintType::RELATIVE_PERSISTENCE:
      this->printMsg(
        "Use the "
        + std::to_string(static_cast<int>(100 * (1 - this->MinPersistence)))
        + "% most persistent pairs of each diagram");
      break;
  }

  std::vector<std::vector<double>> distMat(nInputs);
  if(this->Constraint == ConstraintType::FULL_DIAGRAMS) {
    getDiagramsDistMat(nInputs, distMat, bidder_diagrams_min,
                       bidder_diagrams_sad, bidder_diagrams_max);
  } else {
    if(this->do_min_) {
      enrichCurrentBidderDiagrams(
        bidder_diagrams_min, current_bidder_diagrams_min, maxDiagPersistence);
    }
    if(this->do_sad_) {
      enrichCurrentBidderDiagrams(
        bidder_diagrams_sad, current_bidder_diagrams_sad, maxDiagPersistence);
    }
    if(this->do_max_) {
      enrichCurrentBidderDiagrams(
        bidder_diagrams_max, current_bidder_diagrams_max, maxDiagPersistence);
    }
    getDiagramsDistMat(nInputs, distMat, current_bidder_diagrams_min,
                       current_bidder_diagrams_sad,
                       current_bidder_diagrams_max);
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);

  return distMat;
}

double PersistenceDiagramDistanceMatrix::getMostPersistent(
  const std::vector<BidderDiagram<double>> &bidder_diags) const {

  double max_persistence = 0;

  for(unsigned int i = 0; i < bidder_diags.size(); ++i) {
    for(int j = 0; j < bidder_diags[i].size(); ++j) {
      const double persistence = bidder_diags[i].get(j).getPersistence();
      if(persistence > max_persistence) {
        max_persistence = persistence;
      }
    }
  }

  return max_persistence;
}

double PersistenceDiagramDistanceMatrix::computeDistance(
  const BidderDiagram<double> &D1, const BidderDiagram<double> &D2) const {

  GoodDiagram<double> D2_bis{};
  for(int i = 0; i < D2.size(); i++) {
    const Bidder<double> &b = D2.get(i);
    Good<double> g(b.x_, b.y_, b.isDiagonal(), D2_bis.size());
    g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
    g.setPrice(0);
    D2_bis.addGood(g);
  }

  Auction<double> auction(
    this->Wasserstein, this->Alpha, this->Lambda, this->DeltaLim, true);
  auction.BuildAuctionDiagrams(&D1, &D2_bis);
  return auction.run();
}

void PersistenceDiagramDistanceMatrix::getDiagramsDistMat(
  const size_t nInputs,
  std::vector<std::vector<double>> &distanceMatrix,
  const std::vector<BidderDiagram<double>> &diags_min,
  const std::vector<BidderDiagram<double>> &diags_sad,
  const std::vector<BidderDiagram<double>> &diags_max) const {

  distanceMatrix.resize(nInputs);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs; ++i) {
    distanceMatrix[i].resize(nInputs);

    // matrix diagonal
    distanceMatrix[i][i] = 0.0;

    for(size_t j = i + 1; j < nInputs; ++j) {
      double distance{};

      if(this->do_min_) {
        auto &dimin = diags_min[i];
        auto &djmin = diags_min[j];
        distance += computeDistance(dimin, djmin);
      }
      if(this->do_sad_) {
        auto &disad = diags_sad[i];
        auto &djsad = diags_sad[j];
        distance += computeDistance(disad, djsad);
      }
      if(this->do_max_) {
        auto &dimax = diags_max[i];
        auto &djmax = diags_max[j];
        distance += computeDistance(dimax, djmax);
      }

      distanceMatrix[i][j] = distance;
    }
  }

  // distance matrix is symmetric
  for(size_t i = 0; i < nInputs; ++i) {
    for(size_t j = i + 1; j < nInputs; ++j) {
      distanceMatrix[j][i] = distanceMatrix[i][j];
    }
  }
}

void PersistenceDiagramDistanceMatrix::setBidderDiagrams(
  const size_t nInputs,
  std::vector<std::vector<DiagramTuple>> &inputDiagrams,
  std::vector<BidderDiagram<double>> &bidder_diags) const {

  bidder_diags.resize(nInputs);

  for(size_t i = 0; i < nInputs; i++) {
    auto &diag = inputDiagrams[i];
    auto &bidders = bidder_diags[i];

    for(size_t j = 0; j < diag.size(); j++) {
      // Add bidder to bidders
      Bidder<double> b(diag[j], j, this->Lambda);
      b.setPositionInAuction(bidders.size());
      bidders.addBidder(b);
      if(b.isDiagonal() || b.x_ == b.y_) {
        this->printMsg("Diagonal point in diagram " + std::to_string(i) + "!",
                       ttk::debug::Priority::DETAIL);
      }
    }
  }
}

void PersistenceDiagramDistanceMatrix::enrichCurrentBidderDiagrams(
  const std::vector<BidderDiagram<double>> &bidder_diags,
  std::vector<BidderDiagram<double>> &current_bidder_diags,
  const std::vector<double> &maxDiagPersistence) const {

  current_bidder_diags.resize(bidder_diags.size());
  const auto nInputs = current_bidder_diags.size();

  if(this->Constraint == ConstraintType::ABSOLUTE_PERSISTENCE
     || this->Constraint == ConstraintType::RELATIVE_PERSISTENCE) {
    for(size_t i = 0; i < nInputs; ++i) {
      for(int j = 0; j < bidder_diags[i].size(); ++j) {
        auto b = bidder_diags[i].get(j);

        if( // filter out pairs below absolute persistence threshold
          (this->Constraint == ConstraintType::ABSOLUTE_PERSISTENCE
           && b.getPersistence() > this->MinPersistence)
          || // filter out pairs below persistence threshold relative to
          // current diagram global min-max pair
          (this->Constraint == ConstraintType::RELATIVE_PERSISTENCE
           && b.getPersistence()
                > this->MinPersistence * maxDiagPersistence[i])) {
          b.id_ = current_bidder_diags[i].size();
          b.setPositionInAuction(current_bidder_diags[i].size());
          current_bidder_diags[i].addBidder(b);
        }
      }
    }

  } else if(this->Constraint == ConstraintType::NUMBER_PAIRS) {
    for(size_t i = 0; i < nInputs; ++i) {
      // index of pair in bidder_diags + corresponding persistence value
      using PairIV = std::pair<size_t, double>;
      std::vector<PairIV> pairs(bidder_diags[i].size());
      for(int j = 0; j < bidder_diags[i].size(); ++j) {
        pairs.emplace_back(j, bidder_diags[i].get(j).getPersistence());
      }
      // sort diagram pairs by persistence value
      std::sort(
        pairs.begin(), pairs.end(),
        [](const PairIV &a, const PairIV &b) { return a.second < b.second; });
      const auto firstId = pairs.size() > this->MaxNumberOfPairs
                             ? pairs.size() - this->MaxNumberOfPairs
                             : 0;
      // take the last (most persistent) pairs
      for(size_t j = firstId; j < pairs.size(); ++j) {
        auto b = bidder_diags[i].get(pairs[j].first);
        b.id_ = current_bidder_diags[i].size();
        b.setPositionInAuction(current_bidder_diags[i].size());
        current_bidder_diags[i].addBidder(b);
      }
    }
  }
}
