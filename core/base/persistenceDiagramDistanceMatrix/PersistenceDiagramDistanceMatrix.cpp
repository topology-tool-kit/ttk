#include <algorithm>
#include <limits>

#include <PersistenceDiagramAuction.h>
#include <PersistenceDiagramDistanceMatrix.h>

using namespace ttk;

std::vector<std::vector<double>> PersistenceDiagramDistanceMatrix::execute(
  const std::vector<Diagram> &intermediateDiagrams,
  const std::array<size_t, 2> &nInputs) const {

  Timer tm{};

  const auto nDiags = intermediateDiagrams.size();

  if(do_min_ && do_sad_ && do_max_) {
    this->printMsg("Processing all critical pairs types");
  } else if(do_min_) {
    this->printMsg("Processing only MIN-SAD pairs");
  } else if(do_sad_) {
    this->printMsg("Processing only SAD-SAD pairs");
  } else if(do_max_) {
    this->printMsg("Processing only SAD-MAX pairs");
  }

  std::vector<Diagram> inputDiagramsMin(nDiags);
  std::vector<Diagram> inputDiagramsSad(nDiags);
  std::vector<Diagram> inputDiagramsMax(nDiags);

  std::vector<BidderDiagram<double>> bidder_diagrams_min{};
  std::vector<BidderDiagram<double>> bidder_diagrams_sad{};
  std::vector<BidderDiagram<double>> bidder_diagrams_max{};
  std::vector<BidderDiagram<double>> current_bidder_diagrams_min{};
  std::vector<BidderDiagram<double>> current_bidder_diagrams_sad{};
  std::vector<BidderDiagram<double>> current_bidder_diagrams_max{};

  // Store the persistence of the global min-max pair
  std::vector<double> maxDiagPersistence(nDiags);

  // Create diagrams for min, saddle and max persistence pairs
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDiags; i++) {
    const Diagram &CTDiagram = intermediateDiagrams[i];

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
    setBidderDiagrams(nDiags, inputDiagramsMin, bidder_diagrams_min);
  }
  if(this->do_sad_) {
    setBidderDiagrams(nDiags, inputDiagramsSad, bidder_diagrams_sad);
  }
  if(this->do_max_) {
    setBidderDiagrams(nDiags, inputDiagramsMax, bidder_diagrams_max);
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
    case ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG:
      this->printMsg(
        "Using the "
        + std::to_string(static_cast<int>(100 * (1 - this->MinPersistence)))
        + "% most persistent pairs of every diagram");
      break;
    case ConstraintType::RELATIVE_PERSISTENCE_GLOBAL:
      this->printMsg(
        "Using the "
        + std::to_string(static_cast<int>(100 * (1 - this->MinPersistence)))
        + "% most persistent pairs of all diagrams");
      break;
  }

  std::vector<std::vector<double>> distMat{};
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

template <typename T>
double PersistenceDiagramDistanceMatrix::getMostPersistent(
  const std::vector<T> &bidder_diags) const {

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

template <typename T>
double
  PersistenceDiagramDistanceMatrix::computePowerDistance(const T &D1,
                                                         const T &D2) const {

  GoodDiagram<double> D2_bis{};
  for(int i = 0; i < D2.size(); i++) {
    const auto &b = D2.get(i);
    Good<double> g(b.x_, b.y_, b.isDiagonal(), D2_bis.size());
    g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
    g.setPrice(0);
    D2_bis.addGood(g);
  }

  PersistenceDiagramAuction<double> auction(
    this->Wasserstein, this->Alpha, this->Lambda, this->DeltaLim, true);
  auction.BuildAuctionDiagrams(&D1, &D2_bis);
  return auction.run();
}

template <typename T>
void PersistenceDiagramDistanceMatrix::getDiagramsDistMat(
  const std::array<size_t, 2> &nInputs,
  std::vector<std::vector<double>> &distanceMatrix,
  const std::vector<T> &diags_min,
  const std::vector<T> &diags_sad,
  const std::vector<T> &diags_max) const {

  distanceMatrix.resize(nInputs[0]);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs[0]; ++i) {

    if(nInputs[1] == 0) {
      distanceMatrix[i].resize(nInputs[0]);
      // set the matrix diagonal
      distanceMatrix[i][i] = 0.0;
    } else {
      distanceMatrix[i].resize(nInputs[1]);
    }

    const auto getDist = [&](const size_t a, const size_t b) -> double {
      double distance{};
      if(this->do_min_) {
        auto &dimin = diags_min[a];
        auto &djmin = diags_min[b];
        distance += computePowerDistance(dimin, djmin);
      }
      if(this->do_sad_) {
        auto &disad = diags_sad[a];
        auto &djsad = diags_sad[b];
        distance += computePowerDistance(disad, djsad);
      }
      if(this->do_max_) {
        auto &dimax = diags_max[a];
        auto &djmax = diags_max[b];
        distance += computePowerDistance(dimax, djmax);
      }
      return Geometry::pow(distance, 1.0 / this->Wasserstein);
    };

    if(nInputs[1] == 0) {
      // square matrix: only compute the upper triangle (i < j < nInputs[0])
      for(size_t j = i + 1; j < nInputs[0]; ++j) {
        distanceMatrix[i][j] = getDist(i, j);
      }
    } else {
      // rectangular matrix: compute the whole line/column (0 <= j < nInputs[1])
      for(size_t j = 0; j < nInputs[1]; ++j) {
        distanceMatrix[i][j] = getDist(i, j + nInputs[0]);
      }
    }
  }

  if(nInputs[1] == 0) {
    // square distance matrix is symmetric: complete the lower triangle
    for(size_t i = 0; i < nInputs[0]; ++i) {
      for(size_t j = i + 1; j < nInputs[0]; ++j) {
        distanceMatrix[j][i] = distanceMatrix[i][j];
      }
    }
  }
}

template <typename T>
void PersistenceDiagramDistanceMatrix::setBidderDiagrams(
  const size_t nInputs,
  std::vector<Diagram> &inputDiagrams,
  std::vector<T> &bidder_diags) const {

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

template <typename T>
void PersistenceDiagramDistanceMatrix::enrichCurrentBidderDiagrams(
  const std::vector<T> &bidder_diags,
  std::vector<T> &current_bidder_diags,
  const std::vector<double> &maxDiagPersistence) const {

  current_bidder_diags.resize(bidder_diags.size());
  const auto nInputs = current_bidder_diags.size();
  const auto maxPersistence
    = *std::max_element(maxDiagPersistence.begin(), maxDiagPersistence.end());

  if(this->Constraint == ConstraintType::ABSOLUTE_PERSISTENCE
     || this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG
     || this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_GLOBAL) {
    for(size_t i = 0; i < nInputs; ++i) {
      for(int j = 0; j < bidder_diags[i].size(); ++j) {
        auto b = bidder_diags[i].get(j);

        if( // filter out pairs below absolute persistence threshold
          (this->Constraint == ConstraintType::ABSOLUTE_PERSISTENCE
           && b.getPersistence() > this->MinPersistence)
          || // filter out pairs below persistence threshold relative to
          // the most persistent pair *of each diagrams*
          (this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG
           && b.getPersistence() > this->MinPersistence * maxDiagPersistence[i])
          || // filter out pairs below persistence threshold relative to the
             // most persistence pair *in all diagrams*
          (this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_GLOBAL
           && b.getPersistence() > this->MinPersistence * maxPersistence)) {
          b.id_ = current_bidder_diags[i].size();
          b.setPositionInAuction(current_bidder_diags[i].size());
          current_bidder_diags[i].addBidder(b);
        }
      }
    }
    return;
  }

  const double prev_min_persistence = 2.0 * getMostPersistent(bidder_diags);
  double new_min_persistence = 0.0;

  // 1. Get size of the largest current diagram, deduce the maximal number
  // of points to append
  size_t max_diagram_size = 0;
  for(const auto &diag : current_bidder_diags) {
    max_diagram_size
      = std::max(static_cast<size_t>(diag.size()), max_diagram_size);
  }
  size_t max_points_to_add = std::max(
    this->MaxNumberOfPairs, this->MaxNumberOfPairs + max_diagram_size / 10);
  // 2. Get which points can be added, deduce the new minimal persistence
  std::vector<std::vector<int>> candidates_to_be_added(nInputs);
  std::vector<std::vector<size_t>> idx(nInputs);

  for(size_t i = 0; i < nInputs; i++) {
    double local_min_persistence = std::numeric_limits<double>::min();
    std::vector<double> persistences;
    for(int j = 0; j < bidder_diags[i].size(); j++) {
      const auto b = bidder_diags[i].get(j);
      double persistence = b.getPersistence();
      if(persistence >= 0.0 && persistence <= prev_min_persistence) {
        candidates_to_be_added[i].emplace_back(j);
        idx[i].emplace_back(idx[i].size());
        persistences.emplace_back(persistence);
      }
    }
    const auto cmp = [&persistences](const size_t a, const size_t b) {
      return ((persistences[a] > persistences[b])
              || ((persistences[a] == persistences[b]) && (a > b)));
    };
    std::sort(idx[i].begin(), idx[i].end(), cmp);
    const auto size = candidates_to_be_added[i].size();
    if(size >= max_points_to_add) {
      double last_persistence_added
        = persistences[idx[i][max_points_to_add - 1]];
      if(last_persistence_added > local_min_persistence) {
        local_min_persistence = last_persistence_added;
      }
    }
    if(i == 0) {
      new_min_persistence = local_min_persistence;
    } else {
      if(local_min_persistence < new_min_persistence) {
        new_min_persistence = local_min_persistence;
      }
    }
    // 3. Add the points to the current diagrams
    const auto s = candidates_to_be_added[i].size();
    for(size_t j = 0; j < std::min(max_points_to_add, s); j++) {
      auto b = bidder_diags[i].get(candidates_to_be_added[i][idx[i][j]]);
      const double persistence = b.getPersistence();
      if(persistence >= new_min_persistence) {
        b.id_ = current_bidder_diags[i].size();
        b.setPositionInAuction(current_bidder_diags[i].size());
        current_bidder_diags[i].addBidder(b);
      }
    }
  }
}
