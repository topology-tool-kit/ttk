/// \ingroup base
/// \class ttk::PDClustering
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date September 2019
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa PersistenceDiagramClustering

#include <PDClustering.h>

#include <algorithm>
#include <cmath>
#include <cstdlib> /* srand, rand */
#include <iostream>
#include <iterator>
#include <random>

std::vector<int> ttk::PDClustering::execute(
  std::vector<DiagramType> &final_centroids,
  std::vector<std::vector<std::vector<std::vector<MatchingType>>>>
    &all_matchings_per_type_and_cluster) {

  all_matchings_per_type_and_cluster.resize(k_);
  for(int c = 0; c < k_; c++) {
    all_matchings_per_type_and_cluster[c].resize(3);
    for(int i = 0; i < 3; i++) {
      all_matchings_per_type_and_cluster[c][i].resize(numberOfInputs_);
    }
  }
  int matchings_only = false;
  Timer const tm;
  {
    // PARTICULARITIES FOR THE CASE OF ONE UNIQUE CLUSTER
    if(k_ <= 1) {
      use_accelerated_ = false;
      use_kmeanspp_ = false;
      if(numberOfInputs_ == 2 and !forceUseOfAlgorithm_) {
        use_progressive_ = false;
        deterministic_ = true;
        matchings_only = true;
        time_limit_ = 99999999999;
      }
    }

    std::vector<bool *> current_prec;
    current_prec.emplace_back(&precision_min_);
    current_prec.emplace_back(&precision_sad_);
    current_prec.emplace_back(&precision_max_);

    std::vector<bool *> current_dos;
    current_dos.emplace_back(&do_min_);
    current_dos.emplace_back(&do_sad_);
    current_dos.emplace_back(&do_max_);
    bool converged = false;
    std::vector<bool> diagrams_complete(3);
    for(int c = 0; c < 3; c++) {
      diagrams_complete[c] = (!use_progressive_) || (!original_dos[c]);
    }
    n_iterations_ = 0;
    double total_time = 0;

    // double cost = std::numeric_limits<double>::max();
    setBidderDiagrams();
    cost_ = std::numeric_limits<double>::max();
    double min_cost_min = std::numeric_limits<double>::max();
    double min_cost_max = std::numeric_limits<double>::max();
    double min_cost_sad = std::numeric_limits<double>::max();
    // double last_min_cost_obtained = -1;
    double last_min_cost_obtained_min = -1;
    double last_min_cost_obtained_sad = -1;
    double last_min_cost_obtained_max = -1;
    std::vector<double> epsilon0(3);
    std::vector<double> epsilon_candidate(3);
    std::vector<double> rho(3);

    // Getting current diagrams (with only at most min_points_to_add points)
    std::vector<double> max_persistence(3);
    std::vector<double> lowest_persistence(3);
    std::vector<double> min_persistence(3);

    const auto square = [](const double a) { return a * a; };

    for(int i_crit = 0; i_crit < 3; i_crit++) {
      max_persistence[i_crit] = 2 * getMostPersistent(i_crit);
      lowest_persistence[i_crit] = getLessPersistent(i_crit);
      min_persistence[i_crit] = 0;
      epsilon_[i_crit] = square(0.5 * max_persistence[i_crit])
                         / 8.; // max_persistence actually holds 2 times the
                               // highest persistence
      epsilon0[i_crit] = epsilon_[i_crit];
    }

    std::vector<int> min_points_to_add(3);
    min_points_to_add[0] = 10;
    min_points_to_add[1] = 10;
    min_points_to_add[2] = 10;

    if(use_progressive_) {
    } else {
      min_points_to_add[0] = std::numeric_limits<int>::max();
      min_points_to_add[1] = std::numeric_limits<int>::max();
      min_points_to_add[2] = std::numeric_limits<int>::max();
    }
    std::vector<std::vector<double>> min_diag_price(3);
    std::vector<std::vector<double>> min_off_diag_price(3);
    for(int c = 0; c < 3; ++c) {
      for(int i = 0; i < numberOfInputs_; i++) {
        min_diag_price[c].emplace_back(0);
        min_off_diag_price[c].emplace_back(0);
      }
    }
    min_persistence = enrichCurrentBidderDiagrams(
      max_persistence, min_persistence, min_diag_price, min_off_diag_price,
      min_points_to_add, false, true);

    for(int c = 0; c < 3; c++) {
      if(min_persistence[c] <= lowest_persistence[c]) {
        diagrams_complete[c] = true;
      } // max_persistence actually holds 2 times the highest persistence
    }
    bool all_diagrams_complete
      = diagrams_complete[0] && diagrams_complete[1] && diagrams_complete[2];
    if(all_diagrams_complete) {
      use_progressive_ = false;
    }

    // Initializing centroids and clusters
    if(use_kmeanspp_) {
      initializeCentroidsKMeanspp();
    } else {
      initializeCentroids();
    }
    initializeEmptyClusters();
    if(use_accelerated_) {
      initializeAcceleratedKMeans();
      getCentroidDistanceMatrix();
      acceleratedUpdateClusters();
    } else {
      updateClusters();
      old_clustering_ = clustering_;
    }
    if(debugLevel_ > 3 && k_ > 1) {
      printMsg("Initial Clustering: ");
      printClustering();
    }
    initializeBarycenterComputers(min_persistence);
    while(!converged || (!all_diagrams_complete && use_progressive_)) {
      Timer t_inside;
      {
        n_iterations_++;

        for(int i_crit = 0; i_crit < 3; i_crit++) {
          if(*(current_dos[i_crit])) {
            rho[i_crit] = min_persistence[i_crit] > 0
                            ? std::sqrt(8.0 * epsilon_[i_crit])
                            : -1;
          }
        }

        if(use_progressive_ && n_iterations_ > 1) {

          do_min_ = do_min_ && (min_persistence[0] > rho[0]);
          do_sad_ = do_sad_ && (min_persistence[1] > rho[1]);
          do_max_ = do_max_ && (min_persistence[2] > rho[2]);

          for(int i_crit = 0; i_crit < 3; i_crit++) {
            if(*(current_dos[i_crit])) {
              epsilon_candidate[i_crit] = square(min_persistence[i_crit]) / 8.;
              if(epsilon_candidate[i_crit] > epsilon_[i_crit]) {
                // Should always be the case except if min_persistence is
                // equal to zero
                epsilon_[i_crit] = epsilon_candidate[i_crit];
              }
            }
          }

          if(epsilon_[0] < 5e-5) {
            // Add all remaining points for final convergence.
            // rho[0] = 0;
            min_persistence[0] = 0;
            min_points_to_add[0] = std::numeric_limits<int>::max();
          }
          if(epsilon_[1] < 5e-5) {
            // Add all remaining points for final convergence.
            // rho[1] = 0;
            min_persistence[1] = 0;
            min_points_to_add[1] = std::numeric_limits<int>::max();
          }
          if(epsilon_[2] < 5e-5) {
            // Add all remaining points for final convergence.
            // rho[2] = 0;
            min_persistence[2] = 0;
            min_points_to_add[2] = std::numeric_limits<int>::max();
          }

          if(do_min_ || do_sad_ || do_max_) {
            min_persistence = enrichCurrentBidderDiagrams(
              min_persistence, rho, min_diag_price, min_off_diag_price,
              min_points_to_add, true, false);
          }
          barycenter_inputs_reset_flag = true;

          for(int i_crit = 0; i_crit < 3; i_crit++) {
            if(*(current_dos[i_crit])) {
              if(min_persistence[i_crit] <= lowest_persistence[i_crit]) {
                diagrams_complete[i_crit] = true;
              }
            }
          }

          if(diagrams_complete[0] && diagrams_complete[1]
             && diagrams_complete[2]) {
            use_progressive_ = false;
            all_diagrams_complete = true;
          }

          resetDosToOriginalValues();
        }
        std::vector<double> max_shift_vec = updateCentroidsPosition(
          &min_off_diag_price, &min_diag_price,
          all_matchings_per_type_and_cluster, matchings_only);
        if(do_min_ && !UseDeltaLim_) {
          precision_min_ = (epsilon_[0] < epsilon0[0] / 500.);
        }
        if(do_sad_ && !UseDeltaLim_) {
          precision_sad_ = (epsilon_[1] < epsilon0[1] / 500.);
        }
        if(do_max_ && !UseDeltaLim_) {
          precision_max_ = (epsilon_[2] < epsilon0[2] / 500.);
        }

        for(int i_crit = 0; i_crit < 3; i_crit++) {
          if(*(current_dos[i_crit]) /*&& (!*(current_prec[i_crit]) || !diagrams_complete[i_crit] ) */) {
            epsilon_candidate[i_crit] = std::min(
              std::max(max_shift_vec[i_crit] / 8., epsilon_[i_crit] / 5.),
              epsilon0[i_crit] / square(n_iterations_));

            if((epsilon_candidate[i_crit] < epsilon_[i_crit]
                && !diagrams_complete[i_crit])
               || diagrams_complete[i_crit]) {
              epsilon_[i_crit] = epsilon_candidate[i_crit];
            } else {
              epsilon_[i_crit] *= 0.95;
            }
          }
        }

        if(epsilon_[0] < epsilon_min_ /*&& diagrams_complete[0]*/) {
          this->printMsg("[min barycenter] epsilon under minimal value ",
                         debug::Priority::VERBOSE);
          do_min_ = false;
          epsilon_[0] = epsilon_min_;
          diagrams_complete[0] = true;
        }
        if(epsilon_[1] < epsilon_min_ /*&& diagrams_complete[1]*/) {
          this->printMsg("[sad barycenter] epsilon under minimal value ",
                         debug::Priority::VERBOSE);
          do_sad_ = false;
          epsilon_[1] = epsilon_min_;
          diagrams_complete[1] = true;
        }
        if(epsilon_[2] < epsilon_min_ /*&& diagrams_complete[2]*/) {
          this->printWrn("[max barycenter] epsilon under minimal value ");
          do_max_ = false;
          epsilon_[2] = epsilon_min_;
          diagrams_complete[2] = true;
        }

        if(diagrams_complete[0] && diagrams_complete[1]
           && diagrams_complete[2]) {
          use_progressive_ = false;
          all_diagrams_complete = true;
        }
        if(use_accelerated_) {
          acceleratedUpdateClusters();
        } else {
          // updateClusters();
        }

        precision_criterion_
          = precision_min_ && precision_sad_ && precision_max_;
        bool const precision_criterion_reached = precision_criterion_;

        this->printMsg("Iteration " + std::to_string(n_iterations_)
                         + " epsilon " + std::to_string(epsilon_[0]) + " "
                         + std::to_string(epsilon_[1]) + " "
                         + std::to_string(epsilon_[2]),
                       debug::Priority::VERBOSE);
        this->printMsg(" complete " + std::to_string(diagrams_complete[0]) + " "
                         + std::to_string(diagrams_complete[1]) + " "
                         + std::to_string(diagrams_complete[2]),
                       debug::Priority::VERBOSE);
        this->printMsg(" precision " + std::to_string(precision_min_) + " "
                         + std::to_string(precision_sad_) + " "
                         + std::to_string(precision_max_),
                       debug::Priority::VERBOSE);
        this->printMsg(" cost " + std::to_string(cost_min_) + " "
                         + std::to_string(cost_sad_) + " "
                         + std::to_string(cost_max_),
                       debug::Priority::VERBOSE);

        if(cost_min_ < min_cost_min && n_iterations_ > 2
           && diagrams_complete[0] /*&& precision_min_*/) {
          min_cost_min = cost_min_;
          last_min_cost_obtained_min = 0;
        } else if(n_iterations_ > 2 && precision_min_ && diagrams_complete[0]) {
          last_min_cost_obtained_min += 1;
          if(last_min_cost_obtained_min > 1) {
            do_min_ = false;
          }
        }

        if(cost_sad_ < min_cost_sad && n_iterations_ > 2
           && diagrams_complete[1] /*&& precision_sad_*/) {
          min_cost_sad = cost_sad_;
          last_min_cost_obtained_sad = 0;
        } else if(n_iterations_ > 2 && precision_sad_ && diagrams_complete[1]) {
          last_min_cost_obtained_sad += 1;
          if(last_min_cost_obtained_sad > 1 && diagrams_complete[1]) {
            do_sad_ = false;
          }
        }

        if(cost_max_ < min_cost_max && n_iterations_ > 2
           && diagrams_complete[2] /*&& precision_max_*/) {
          min_cost_max = cost_max_;
          last_min_cost_obtained_max = 0;
        } else if(n_iterations_ > 2 && precision_max_ && diagrams_complete[2]) {
          last_min_cost_obtained_max += 1;
          if(last_min_cost_obtained_max > 1 && diagrams_complete[2]) {
            do_max_ = false;
          }
        }

        if(debugLevel_ > 5) {
          this->printMsg("Clustering result:", debug::Priority::DETAIL);
          printClustering();
        }
        converged = converged
                    || (all_diagrams_complete && !do_min_ && !do_sad_
                        && !do_max_ && (precision_criterion_reached))
                    || (precision_criterion_reached && UseDeltaLim_
                        && numberOfInputs_ == 2);
      }

      total_time
        += t_inside.getElapsedTime(); // - t_real_cost.getElapsedTime();

      if(total_time + t_inside.getElapsedTime() > 0.9 * time_limit_) {
        min_cost_min = cost_min_;
        min_cost_sad = cost_sad_;
        min_cost_max = cost_max_;
        converged = true;
      }
      if(total_time > 0.1 * time_limit_) {
        all_diagrams_complete = true;
        diagrams_complete[0] = true;
        diagrams_complete[1] = true;
        diagrams_complete[2] = true;
        use_progressive_ = false;
      }
      if(debugLevel_ > 4) {
        std::cout << "== Iteration " << n_iterations_
                  << " == complete : " << all_diagrams_complete
                  << " , progressive : " << use_progressive_
                  << " , converged : " << converged << std::endl;
      }
    }
    resetDosToOriginalValues();

    // display results
    std::vector<std::vector<std::string>> const rows{
      {" Min-saddle cost", std::to_string(cost_min_)},
      {" Saddle-saddle cost", std::to_string(cost_sad_)},
      {" Saddle-max cost", std::to_string(cost_max_)},
      {matchings_only ? "Wasserstein Distance" : "Final Cost",
       std::to_string(cost_min_ + cost_sad_ + cost_max_)},
    };
    this->printMsg(rows);

    if(!use_progressive_ && k_ > 1) {
      clustering_ = old_clustering_; // reverting to last clustering
    }
    invertClusters(); // this is to pass the old inverse clustering to the VTK
                      // wrapper
    if(k_ > 1) {
      this->printMsg("Clustering result:");
      printClustering();
    }
  } // End of timer

  // CORRECT MATCHINGS :
  // correctMatchings(all_matchings);
  if(matchings_only) {
    computeBarycenterForTwoGlobal(all_matchings_per_type_and_cluster);
  }
  correctMatchings(all_matchings_per_type_and_cluster);
  // Filling the final centroids for output

  final_centroids.resize(k_);
  centroids_sizes_.resize(k_);
  for(int c = 0; c < k_; c++) {
    centroids_sizes_[c].resize(3);
    if(do_min_)
      centroids_sizes_[c][0] = centroids_min_[c].size();
    if(do_sad_)
      centroids_sizes_[c][1] = centroids_saddle_[c].size();
    if(do_max_)
      centroids_sizes_[c][2] = centroids_max_[c].size();
  }

  for(int c = 0; c < k_; ++c) {
    // if NumberOfClusters > 1, the global pair was duplicated
    // and needs to be removed from the min-saddle problem
    // It is the first pair.
    int const removeFirstPairMin
      = (k_ > 1 and original_dos[0] and original_dos[2]) ? 1 : 0;
    int addedFirstPairMax = 0;
    int addedFirstPairMin = removeFirstPairMin;

    // min-max Pair
    if(removeFirstPairMin or (!do_min_ and do_max_)) {
      Good const &g = centroids_max_[c].at(0);
      const auto &critCoords = g.GetCriticalCoordinates();
      final_centroids[c].emplace_back(PersistencePair{
        CriticalVertex{0, CriticalType::Local_minimum, g.x_, critCoords},
        CriticalVertex{0, CriticalType::Local_maximum, g.y_, critCoords}, 0,
        false});
      addedFirstPairMax = 1;
    } else if(do_min_) {
      Good const &g = centroids_min_[c].at(0);
      const auto &critCoords = g.GetCriticalCoordinates();
      final_centroids[c].emplace_back(PersistencePair{
        CriticalVertex{0, CriticalType::Local_minimum, g.x_, critCoords},
        CriticalVertex{0, CriticalType::Local_maximum, g.y_, critCoords}, 0,
        false});
      addedFirstPairMin = 1;
    }

    if(do_min_) {
      for(size_t i = addedFirstPairMin; i < centroids_min_[c].size(); ++i) {
        Good const &g = centroids_min_[c].at(i);
        const auto &critCoords = g.GetCriticalCoordinates();
        final_centroids[c].emplace_back(PersistencePair{
          CriticalVertex{0, CriticalType::Local_minimum, g.x_, critCoords},
          CriticalVertex{0, CriticalType::Saddle1, g.y_, critCoords}, 0, true});
        if(g.getPersistence() > 1000) {
          this->printMsg("Found a abnormally high persistence in min diagram",
                         debug::Priority::WARNING);
        }
      }
    }

    if(do_sad_) {
      for(size_t i = 0; i < centroids_saddle_[c].size(); ++i) {
        Good const &g = centroids_saddle_[c].at(i);
        const auto &critCoords = g.GetCriticalCoordinates();
        final_centroids[c].emplace_back(PersistencePair{
          CriticalVertex{0, CriticalType::Saddle1, g.x_, critCoords},
          CriticalVertex{0, CriticalType::Saddle2, g.y_, critCoords}, 1, true});
        if(g.getPersistence() > 1000) {
          this->printMsg("Found a abnormally high persistence in sad diagram",
                         debug::Priority::WARNING);
        }
      }
    }

    if(do_max_) {
      for(size_t i = addedFirstPairMax; i < centroids_max_[c].size(); ++i) {
        Good const &g = centroids_max_[c].at(i);
        const auto &critCoords = g.GetCriticalCoordinates();
        ttk::CriticalType saddle_type;
        if(do_sad_)
          saddle_type = ttk::CriticalType::Saddle2;
        else
          saddle_type = ttk::CriticalType::Saddle1;

        final_centroids[c].emplace_back(PersistencePair{
          CriticalVertex{0, saddle_type, g.x_, critCoords},
          CriticalVertex{0, CriticalType::Local_maximum, g.y_, critCoords}, 2,
          true});

        if(g.getPersistence() > 1000) {
          this->printMsg("Found a abnormally high persistence in min diagram",
                         debug::Priority::WARNING);
        }
      }
    }
  }

  if(distanceWritingOptions_ == 1) {
    printDistancesToFile();
  } else if(distanceWritingOptions_ == 2) {
    printRealDistancesToFile();
  }

  return inv_clustering_;
}

void ttk::PDClustering::correctMatchings(
  std::vector<std::vector<std::vector<std::vector<MatchingType>>>>
    &previous_matchings) {
  for(int c = 0; c < k_; c++) {
    for(size_t i = 0; i < clustering_[c].size(); i++) {
      int const diagram_id = clustering_[c][i];
      if(original_dos[0]) {
        // 1. Invert the current_bidder_ids_ vector
        std::vector<int> new_to_old_id(
          current_bidder_diagrams_min_[diagram_id].size(), -1);
        for(size_t j = 0; j < current_bidder_ids_min_[diagram_id].size(); j++) {
          int const new_id = current_bidder_ids_min_[diagram_id][j];
          if(new_id >= 0) {
            new_to_old_id[new_id] = j;
          }
        }
        // 2. Reconstruct the matchings
        std::vector<MatchingType> matchings_diagram_i;
        for(size_t j = 0; j < previous_matchings[c][0][i].size(); j++) {
          MatchingType m = previous_matchings[c][0][i][j];
          int const new_id = std::get<0>(m);
          if(new_id >= 0 && std::get<1>(m) >= 0) {
            std::get<0>(m) = new_to_old_id[new_id];
            matchings_diagram_i.emplace_back(m);
          } else if(std::get<1>(m)
                    >= 0) { // new_id < 0 corresponds to a diagonal matching
            std::get<0>(m) = -1;
            matchings_diagram_i.emplace_back(m);
          }
        }
        previous_matchings[c][0][i].resize(matchings_diagram_i.size());
        previous_matchings[c][0][i] = matchings_diagram_i;
      }
      if(original_dos[1]) {
        // 1. Invert the current_bidder_ids_ vector
        std::vector<int> new_to_old_id(
          current_bidder_diagrams_saddle_[diagram_id].size());
        for(size_t j = 0; j < current_bidder_ids_sad_[diagram_id].size(); j++) {
          int const new_id = current_bidder_ids_sad_[diagram_id][j];
          if(new_id >= 0) {
            new_to_old_id[new_id] = j;
          }
        }
        // 2. Reconstruct the matchings
        int zero_done = 0;
        std::vector<MatchingType> matchings_diagram_i;
        for(size_t j = 0; j < previous_matchings[c][1][i].size(); j++) {
          MatchingType m = previous_matchings[c][1][i][j];
          int const new_id = std::get<0>(m);
          if(new_id >= 0 && std::get<1>(m) >= 0) {
            int const old_id = new_to_old_id[new_id];
            if(old_id > 0) {
              std::get<0>(m) = old_id;
              matchings_diagram_i.emplace_back(m);
            } else {
              if(!zero_done) {
                zero_done = 1;
                std::get<0>(m) = old_id;
                matchings_diagram_i.emplace_back(m);
              }
            }
          } else if(std::get<1>(m)
                    >= 0) { // new_id < 0 corresponds to a diagonal matching
            std::get<0>(m) = -1;
            matchings_diagram_i.emplace_back(m);
          }
        }
        previous_matchings[c][1][i].resize(matchings_diagram_i.size());
        previous_matchings[c][1][i] = matchings_diagram_i;
      }
      if(original_dos[2]) {
        // 1. Invert the current_bidder_ids_ vector
        std::vector<int> new_to_old_id(
          current_bidder_diagrams_max_[diagram_id].size());
        for(size_t j = 0; j < current_bidder_ids_max_[diagram_id].size(); j++) {
          int const new_id = current_bidder_ids_max_[diagram_id][j];
          if(new_id >= 0) {
            new_to_old_id[new_id] = j;
          }
        }
        // 2. Reconstruct the matchings
        int zero_done = 0;
        std::vector<MatchingType> matchings_diagram_i;
        for(size_t j = 0; j < previous_matchings[c][2][i].size(); j++) {
          MatchingType m = previous_matchings[c][2][i][j];
          int const new_id = std::get<0>(m);
          if(new_id >= 0 && std::get<1>(m) >= 0) {
            int const old_id = new_to_old_id[new_id];
            if(old_id > 0) {
              std::get<0>(m) = old_id;
              matchings_diagram_i.emplace_back(m);
            } else {
              if(!zero_done) {
                zero_done = 1;
                std::get<0>(m) = old_id;
                matchings_diagram_i.emplace_back(m);
              }
            }
          } else if(std::get<1>(m)
                    >= 0) { // new_id < 0 corresponds to a diagonal matching
            std::get<0>(m) = -1;
            matchings_diagram_i.emplace_back(m);
          }
        }
        previous_matchings[c][2][i].resize(matchings_diagram_i.size());
        previous_matchings[c][2][i] = matchings_diagram_i;
      }
    }
  }
}

void ttk::PDClustering::printMatchings(
  std::vector<std::vector<std::vector<MatchingType>>> matchings) {
  std::cout << "\n MATCHINGS : " << std::endl;
  for(int d = 0; d < 3; d++) {
    if(original_dos[d]) {
      std::cout << "\n Diagram type : " << d << std::endl;
      for(size_t i = 0; i < matchings[d].size(); i++) {
        std::cout << " diagram " << i << " : ";
        for(size_t j = 0; j < matchings[d][i].size(); j++) {
          std::cout << std::get<0>(matchings[d][i][j]) << " ";
          std::cout << std::get<1>(matchings[d][i][j]) << " ";
          std::cout << std::get<2>(matchings[d][i][j]) << "  |   ";
        }
        std::cout << "\n";
      }
    }
  }
}

std::vector<std::vector<int>> ttk::PDClustering::get_centroids_sizes() {
  return centroids_sizes_;
}

double ttk::PDClustering::getMostPersistent(int type) {
  double max_persistence = 0;

  if(do_min_ && (type == -1 || type == 0)) {
    for(size_t i = 0; i < bidder_diagrams_min_.size(); ++i) {
      for(size_t j = 0; j < bidder_diagrams_min_[i].size(); ++j) {
        Bidder const b = bidder_diagrams_min_[i].at(j);
        double const persistence = b.getPersistence();
        if(persistence > max_persistence) {
          max_persistence = persistence;
        }
      }
    }
  }

  if(do_sad_ && (type == -1 || type == 1)) {
    for(size_t i = 0; i < bidder_diagrams_saddle_.size(); ++i) {
      for(size_t j = 0; j < bidder_diagrams_saddle_[i].size(); ++j) {
        Bidder const b = bidder_diagrams_saddle_[i].at(j);
        double const persistence = b.getPersistence();
        if(persistence > max_persistence) {
          max_persistence = persistence;
        }
      }
    }
  }

  if(do_max_ && (type == -1 || type == 2)) {
    for(size_t i = 0; i < bidder_diagrams_max_.size(); ++i) {
      for(size_t j = 0; j < bidder_diagrams_max_[i].size(); ++j) {
        Bidder const b = bidder_diagrams_max_[i].at(j);
        double const persistence = b.getPersistence();
        if(persistence > max_persistence) {
          max_persistence = persistence;
        }
      }
    }
  }
  return max_persistence;
}

double ttk::PDClustering::getLessPersistent(int type) {
  // type == -1 : query the min of all the types of diagrams.
  // type = 0 : min,  1 : sad,   2 : max

  double min_persistence = std::numeric_limits<double>::max();
  if(do_min_ && (type == -1 || type == 0)) {
    for(size_t i = 0; i < bidder_diagrams_min_.size(); ++i) {
      for(size_t j = 0; j < bidder_diagrams_min_[i].size(); ++j) {
        Bidder const b = bidder_diagrams_min_[i].at(j);
        double const persistence = b.getPersistence();
        if(persistence < min_persistence) {
          min_persistence = persistence;
        }
      }
    }
  }

  if(do_sad_ && (type == -1 || type == 1)) {
    for(size_t i = 0; i < bidder_diagrams_saddle_.size(); ++i) {
      for(size_t j = 0; j < bidder_diagrams_saddle_[i].size(); ++j) {
        Bidder const b = bidder_diagrams_saddle_[i].at(j);
        double const persistence = b.getPersistence();
        if(persistence < min_persistence) {
          min_persistence = persistence;
        }
      }
    }
  }

  if(do_max_ && (type == -1 || type == 2)) {
    for(size_t i = 0; i < bidder_diagrams_max_.size(); ++i) {
      for(size_t j = 0; j < bidder_diagrams_max_[i].size(); ++j) {
        Bidder const b = bidder_diagrams_max_[i].at(j);
        double const persistence = b.getPersistence();
        if(persistence < min_persistence) {
          min_persistence = persistence;
        }
      }
    }
  }
  return min_persistence;
}

std::vector<std::vector<double>> ttk::PDClustering::getMinPrices() {
  std::vector<std::vector<double>> min_prices(3);

  if(original_dos[0]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[0].emplace_back(std::numeric_limits<double>::max());
      for(size_t j = 0; j < centroids_with_price_min_[i].size(); ++j) {
        Good const g = centroids_with_price_min_[i].at(j);
        double const price = g.getPrice();
        if(price < min_prices[0][i]) {
          min_prices[0][i] = price;
        }
      }
    }
  }

  if(original_dos[1]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[1].emplace_back(std::numeric_limits<double>::max());
      for(size_t j = 0; j < centroids_with_price_saddle_[i].size(); ++j) {
        Good const g = centroids_with_price_saddle_[i].at(j);
        double const price = g.getPrice();
        if(price < min_prices[1][i]) {
          min_prices[1][i] = price;
        }
      }
    }
  }

  if(original_dos[2]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[2].emplace_back(std::numeric_limits<double>::max());
      for(size_t j = 0; j < centroids_with_price_max_[i].size(); ++j) {
        Good const g = centroids_with_price_max_[i].at(j);
        double const price = g.getPrice();
        if(price < min_prices[2][i]) {
          min_prices[2][i] = price;
        }
      }
    }
  }

  return min_prices;
}

std::vector<std::vector<double>> ttk::PDClustering::getMinDiagonalPrices() {
  std::vector<std::vector<double>> min_prices(3);
  if(original_dos[0]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[0].emplace_back(std::numeric_limits<double>::max());
      for(size_t j = 0; j < current_bidder_diagrams_min_[i].size(); ++j) {
        Bidder const b = current_bidder_diagrams_min_[i].at(j);
        double const price = b.diagonal_price_;
        if(price < min_prices[0][i]) {
          min_prices[0][i] = price;
        }
      }
      if(min_prices[0][i] >= std::numeric_limits<double>::max() / 2.) {
        min_prices[0][i] = 0;
      }
    }
  }

  if(original_dos[1]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[1].emplace_back(std::numeric_limits<double>::max());
      for(size_t j = 0; j < current_bidder_diagrams_saddle_[i].size(); ++j) {
        Bidder const b = current_bidder_diagrams_saddle_[i].at(j);
        double const price = b.diagonal_price_;
        if(price < min_prices[1][i]) {
          min_prices[1][i] = price;
        }
      }
      if(min_prices[1][i] >= std::numeric_limits<double>::max() / 2.) {
        min_prices[1][i] = 0;
      }
    }
  }

  if(original_dos[2]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[2].emplace_back(std::numeric_limits<double>::max());
      for(size_t j = 0; j < current_bidder_diagrams_max_[i].size(); ++j) {
        Bidder const b = current_bidder_diagrams_max_[i].at(j);
        double const price = b.diagonal_price_;
        if(price < min_prices[2][i]) {
          min_prices[2][i] = price;
        }
      }
      if(min_prices[2][i] >= std::numeric_limits<double>::max() / 2.) {
        min_prices[2][i] = 0;
      }
    }
  }
  return min_prices;
}

double ttk::PDClustering::computeDistance(const BidderDiagram &D1,
                                          const BidderDiagram &D2,
                                          const double delta_lim) {
  GoodDiagram const D2_bis = diagramToCentroid(D2);
  return computeDistance(D1, D2_bis, delta_lim);
}

double ttk::PDClustering::computeDistance(const BidderDiagram &D1,
                                          const GoodDiagram &D2,
                                          const double delta_lim) {
  std::vector<MatchingType> matchings;
  const auto D2_bis = centroidWithZeroPrices(D2);
  PersistenceDiagramAuction auction(
    wasserstein_, geometrical_factor_, lambda_, delta_lim, use_kdtree_);
  auction.BuildAuctionDiagrams(D1, D2_bis);
  double const cost = auction.run(matchings);
  return cost;
}

double ttk::PDClustering::computeDistance(BidderDiagram *const D1,
                                          const GoodDiagram *const D2,
                                          const double delta_lim) {
  std::vector<MatchingType> matchings;
  PersistenceDiagramAuction auction(
    wasserstein_, geometrical_factor_, lambda_, delta_lim, use_kdtree_);
  int const size1 = D1->size();
  auction.BuildAuctionDiagrams(*D1, *D2);
  double const cost = auction.run(matchings);
  // Diagonal Points were added in the original diagram. The following line
  // removes them.
  D1->resize(size1);
  return cost;
}

double ttk::PDClustering::computeDistance(const GoodDiagram &D1,
                                          const GoodDiagram &D2,
                                          const double delta_lim) {
  BidderDiagram const D1_bis = centroidToDiagram(D1);
  return computeDistance(D1_bis, D2, delta_lim);
}

ttk::GoodDiagram
  ttk::PDClustering::centroidWithZeroPrices(const GoodDiagram &centroid) {
  GoodDiagram GD = GoodDiagram();
  for(size_t i = 0; i < centroid.size(); i++) {
    Good g = centroid.at(i);
    g.setPrice(0);
    GD.emplace_back(g);
  }
  return GD;
}

ttk::BidderDiagram
  ttk::PDClustering::diagramWithZeroPrices(const BidderDiagram &diagram) {
  BidderDiagram BD = BidderDiagram();
  for(size_t i = 0; i < diagram.size(); i++) {
    Bidder b = diagram.at(i);
    b.setDiagonalPrice(0);
    BD.emplace_back(b);
  }
  return BD;
}

ttk::BidderDiagram
  ttk::PDClustering::centroidToDiagram(const GoodDiagram &centroid) {
  BidderDiagram BD = BidderDiagram();
  for(size_t i = 0; i < centroid.size(); i++) {
    Good g = centroid.at(i);

    Bidder b = Bidder(g.x_, g.y_, g.isDiagonal(), BD.size());
    b.SetCriticalCoordinates(g.coords_[0], g.coords_[1], g.coords_[2]);
    b.setPositionInAuction(BD.size());
    BD.emplace_back(b);
  }
  return BD;
}

ttk::GoodDiagram
  ttk::PDClustering::diagramToCentroid(const BidderDiagram &diagram) {
  GoodDiagram GD = GoodDiagram();
  for(size_t i = 0; i < diagram.size(); i++) {
    Bidder b = diagram.at(i);

    Good g = Good(b.x_, b.y_, b.isDiagonal(), GD.size());
    g.SetCriticalCoordinates(b.coords_[0], b.coords_[1], b.coords_[2]);
    GD.emplace_back(g);
  }
  return GD;
}

void ttk::PDClustering::initializeEmptyClusters() {
  clustering_ = std::vector<std::vector<int>>(k_);
}

void ttk::PDClustering::initializeCentroids() {
  std::vector<int> idx(numberOfInputs_);
  // To perform a random draw with replacement, the vector {1, 2, ...,
  // numberOfInputs_} is shuffled, and we consider its k_ first elements to be
  // the initial centroids.
  for(int i = 0; i < numberOfInputs_; i++) {
    idx[i] = i;
  }
  if(!deterministic_)
    std::shuffle(idx.begin(), idx.end(), std::random_device());

  for(int c = 0; c < k_; c++) {
    if(do_min_) {
      GoodDiagram const centroid_min
        = diagramToCentroid(current_bidder_diagrams_min_[idx[c]]);
      centroids_min_.emplace_back(centroid_min);
    }
    if(do_sad_) {
      GoodDiagram const centroid_sad
        = diagramToCentroid(current_bidder_diagrams_saddle_[idx[c]]);
      centroids_saddle_.emplace_back(centroid_sad);
    }
    if(do_max_) {
      GoodDiagram const centroid_max
        = diagramToCentroid(current_bidder_diagrams_max_[idx[c]]);
      centroids_max_.emplace_back(centroid_max);
    }
  }
}

void ttk::PDClustering::initializeCentroidsKMeanspp() {
  std::vector<int> indexes_clusters;
  int const random_idx = deterministic_ ? 0 : rand() % numberOfInputs_;
  indexes_clusters.emplace_back(random_idx);

  if(do_min_) {
    GoodDiagram const centroid_min
      = diagramToCentroid(current_bidder_diagrams_min_[random_idx]);
    centroids_min_.emplace_back(centroid_min);
  }
  if(do_sad_) {
    GoodDiagram const centroid_sad
      = diagramToCentroid(current_bidder_diagrams_saddle_[random_idx]);
    centroids_saddle_.emplace_back(centroid_sad);
  }
  if(do_max_) {
    GoodDiagram const centroid_max
      = diagramToCentroid(current_bidder_diagrams_max_[random_idx]);
    centroids_max_.emplace_back(centroid_max);
  }

  const auto square = [](const double a) { return a * a; };

  while((int)indexes_clusters.size() < k_) {
    std::vector<double> min_distance_to_centroid(numberOfInputs_);
    std::vector<double> probabilities(numberOfInputs_);

    // Uncomment for a deterministic algorithm
    double maximal_distance = 0;
    int candidate_centroid = 0;

    for(int i = 0; i < numberOfInputs_; i++) {

      min_distance_to_centroid[i] = std::numeric_limits<double>::max();
      if(std::find(indexes_clusters.begin(), indexes_clusters.end(), i)
         != indexes_clusters.end()) {

        min_distance_to_centroid[i] = 0;
      } else {

        for(size_t j = 0; j < indexes_clusters.size(); ++j) {

          double distance = 0;
          if(do_min_) {
            GoodDiagram const centroid_min
              = centroidWithZeroPrices(centroids_min_[j]);
            distance += computeDistance(
              current_bidder_diagrams_min_[i], centroid_min, 0.01);
          }
          if(do_sad_) {
            GoodDiagram const centroid_saddle
              = centroidWithZeroPrices(centroids_saddle_[j]);
            distance += computeDistance(
              current_bidder_diagrams_saddle_[i], centroid_saddle, 0.01);
          }
          if(do_max_) {
            GoodDiagram const centroid_max
              = centroidWithZeroPrices(centroids_max_[j]);
            distance += computeDistance(
              current_bidder_diagrams_max_[i], centroid_max, 0.01);
          }
          if(distance < min_distance_to_centroid[i]) {
            min_distance_to_centroid[i] = distance;
          }
        }
      }
      probabilities[i] = square(min_distance_to_centroid[i]);

      // The following block is useful in case of need for a deterministic
      // algorithm
      if(deterministic_ && min_distance_to_centroid[i] > maximal_distance) {
        maximal_distance = min_distance_to_centroid[i];
        candidate_centroid = i;
      }
    }
    // Comment the following four lines to make it deterministic
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<int> distribution(
      probabilities.begin(), probabilities.end());

    if(!deterministic_) {
      candidate_centroid = distribution(gen);
    }

    indexes_clusters.emplace_back(candidate_centroid);
    if(do_min_) {
      GoodDiagram const centroid_min
        = diagramToCentroid(current_bidder_diagrams_min_[candidate_centroid]);
      centroids_min_.emplace_back(centroid_min);
    }
    if(do_sad_) {
      GoodDiagram const centroid_sad = diagramToCentroid(
        current_bidder_diagrams_saddle_[candidate_centroid]);
      centroids_saddle_.emplace_back(centroid_sad);
    }
    if(do_max_) {
      GoodDiagram const centroid_max
        = diagramToCentroid(current_bidder_diagrams_max_[candidate_centroid]);
      centroids_max_.emplace_back(centroid_max);
    }
  }
}

void ttk::PDClustering::initializeAcceleratedKMeans() {
  // r_ is a vector stating for each diagram if its distance to its centroid is
  // up to date (false) or needs to be recomputed (true)
  r_ = std::vector<bool>(numberOfInputs_);
  // u_ is a vector of upper bounds of the distance of each diagram to its
  // closest centroid
  u_ = std::vector<double>(numberOfInputs_);
  inv_clustering_ = std::vector<int>(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; i++) {
    r_[i] = true;
    u_[i] = std::numeric_limits<double>::max();
    inv_clustering_[i] = -1;
  }
  // l_ is the matrix of lower bounds for the distance from each diagram
  // to each centroid
  l_ = std::vector<std::vector<double>>(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; ++i) {
    l_[i] = std::vector<double>(k_);
    for(int c = 0; c < k_; ++c) {
      l_[i][c] = 0;
    }
  }

  // And d_ is a (K x K) matrix storing the distances between each pair of
  // centroids
  centroidsDistanceMatrix_.resize(k_);
  for(int i = 0; i < k_; ++i) {
    centroidsDistanceMatrix_[i].resize(k_, 0.0);
  }
}

std::vector<std::vector<double>> ttk::PDClustering::getDistanceMatrix() {
  std::vector<std::vector<double>> D(numberOfInputs_);

  for(int i = 0; i < numberOfInputs_; ++i) {
    BidderDiagram D1_min, D1_sad, D1_max;
    if(do_min_) {
      D1_min = diagramWithZeroPrices(current_bidder_diagrams_min_[i]);
    }
    if(do_sad_) {
      D1_sad = diagramWithZeroPrices(current_bidder_diagrams_saddle_[i]);
    }
    if(do_max_) {
      D1_max = diagramWithZeroPrices(current_bidder_diagrams_max_[i]);
    }
    for(int c = 0; c < k_; ++c) {
      GoodDiagram D2_min, D2_sad, D2_max;
      double distance = 0;
      if(do_min_) {
        D2_min = centroids_min_[c];
        distance += computeDistance(D1_min, D2_min, 0.01);
      }
      if(do_sad_) {
        D2_sad = centroids_saddle_[c];
        distance += computeDistance(D1_sad, D2_sad, 0.01);
      }
      if(do_max_) {
        D2_max = centroids_max_[c];
        distance += computeDistance(D1_max, D2_max, 0.01);
      }
      D[i].emplace_back(distance);
    }
  }
  return D;
}

void ttk::PDClustering::getCentroidDistanceMatrix() {
  for(int i = 0; i < k_; ++i) {
    GoodDiagram D1_min, D1_sad, D1_max;
    if(do_min_) {
      D1_min = centroidWithZeroPrices(centroids_min_[i]);
    }
    if(do_sad_) {
      D1_sad = centroidWithZeroPrices(centroids_saddle_[i]);
    }
    if(do_max_) {
      D1_max = centroidWithZeroPrices(centroids_max_[i]);
    }
    for(int j = i + 1; j < k_; ++j) {
      double distance{};
      GoodDiagram D2_min, D2_sad, D2_max;
      if(do_min_) {
        D2_min = centroidWithZeroPrices(centroids_min_[j]);
        distance += computeDistance(D1_min, D2_min, 0.01);
      }
      if(do_sad_) {
        D2_sad = centroidWithZeroPrices(centroids_saddle_[j]);
        distance += computeDistance(D1_sad, D2_sad, 0.01);
      }
      if(do_max_) {
        D2_max = centroidWithZeroPrices(centroids_max_[j]);
        distance += computeDistance(D1_max, D2_max, 0.01);
      }

      centroidsDistanceMatrix_[i][j] = distance;
      centroidsDistanceMatrix_[j][i] = distance;
    }
  }
}

void ttk::PDClustering::computeDistanceToCentroid() {
  distanceToCentroid_.resize(numberOfInputs_);

  for(int i = 0; i < numberOfInputs_; ++i) {
    double const delta_lim{0.01};
    double distance{};
    auto c = inv_clustering_[i];
    if(original_dos[0]) {
      GoodDiagram const centroid_min
        = centroidWithZeroPrices(centroids_min_[c]);
      BidderDiagram const bidder_diag
        = diagramWithZeroPrices(current_bidder_diagrams_min_[i]);
      distance += computeDistance(bidder_diag, centroid_min, delta_lim);
    }
    if(original_dos[1]) {
      GoodDiagram const centroid_saddle
        = centroidWithZeroPrices(centroids_saddle_[c]);
      BidderDiagram const bidder_diag
        = diagramWithZeroPrices(current_bidder_diagrams_saddle_[i]);
      distance += computeDistance(bidder_diag, centroid_saddle, delta_lim);
    }
    if(original_dos[2]) {
      GoodDiagram const centroid_max
        = centroidWithZeroPrices(centroids_max_[c]);
      BidderDiagram const bidder_diag
        = diagramWithZeroPrices(current_bidder_diagrams_max_[i]);
      distance += computeDistance(bidder_diag, centroid_max, delta_lim);
    }
    distanceToCentroid_[i] = distance;
  }
}

void ttk::PDClustering::updateClusters() {
  if(k_ > 1) {
    std::vector<std::vector<double>> distance_matrix = getDistanceMatrix();
    old_clustering_ = clustering_;
    invertClusters();
    initializeEmptyClusters();

    for(int i = 0; i < numberOfInputs_; ++i) {
      double min_distance_to_centroid = std::numeric_limits<double>::max();
      int cluster = -1;
      for(int c = 0; c < k_; ++c) {
        if(distance_matrix[i][c] < min_distance_to_centroid) {
          min_distance_to_centroid = distance_matrix[i][c];
          cluster = c;
        }
      }

      clustering_[cluster].emplace_back(i);
      if(cluster != inv_clustering_[i]) {
        // New centroid attributed to this diagram
        resetDosToOriginalValues();
        barycenter_inputs_reset_flag = true;
        if(do_min_) {
          centroids_with_price_min_[i]
            = centroidWithZeroPrices(centroids_min_[cluster]);
        }
        if(do_sad_) {
          centroids_with_price_saddle_[i]
            = centroidWithZeroPrices(centroids_saddle_[cluster]);
        }
        if(do_max_) {
          centroids_with_price_max_[i]
            = centroidWithZeroPrices(centroids_max_[cluster]);
        }
        inv_clustering_[i] = cluster;
      }
    }

  } else {
    old_clustering_ = clustering_;
    invertClusters();
    initializeEmptyClusters();

    for(int i = 0; i < numberOfInputs_; i++) {
      clustering_[0].emplace_back(i);
      if(n_iterations_ < 1) {
        if(do_min_) {
          centroids_with_price_min_[i]
            = centroidWithZeroPrices(centroids_min_[0]);
        }
        if(do_sad_) {
          centroids_with_price_saddle_[i]
            = centroidWithZeroPrices(centroids_saddle_[0]);
        }
        if(do_max_) {
          centroids_with_price_max_[i]
            = centroidWithZeroPrices(centroids_max_[0]);
        }
      }
      inv_clustering_[i] = 0;
    }
  }
}

void ttk::PDClustering::invertClusters() {
  /// Converts the clustering (vector of vector of diagram's id) into
  /// a vector of size numberOfInputs_ containing the cluster of each input
  /// diagram.

  // Initializes clusters with -1
  inv_clustering_ = std::vector<int>(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; ++i) {
    inv_clustering_[i] = -1;
  }

  // Fill in the clusters
  for(int c = 0; c < k_; ++c) {
    for(size_t j = 0; j < clustering_[c].size(); ++j) {
      int const idx = clustering_[c][j];
      inv_clustering_[idx] = c;
    }
  }
}

void ttk::PDClustering::invertInverseClusters() {
  clustering_ = std::vector<std::vector<int>>(k_);
  for(int i = 0; i < numberOfInputs_; ++i) {
    clustering_[inv_clustering_[i]].emplace_back(i);
  }

  // Check if a cluster was left without diagram
  for(int c = 0; c < k_; ++c) {
    if(clustering_[c].empty()) {
      std::cout << "Problem in invertInverseClusters()... \nCluster " << c
                << " was left with no diagram attached to it... " << std::endl;
    }
  }
}

void ttk::PDClustering::acceleratedUpdateClusters() {
  // Step 1
  getCentroidDistanceMatrix();
  old_clustering_ = clustering_;
  // self.old_clusters = copy.copy(self.clusters)
  invertClusters();
  initializeEmptyClusters();
  bool const do_min = original_dos[0];
  bool const do_sad = original_dos[1];
  bool const do_max = original_dos[2];

  for(int i = 0; i < numberOfInputs_; ++i) {
    // Step 3 find potential changes of clusters
    BidderDiagram D1_min, D1_sad, D1_max;
    if(do_min) {
      D1_min = diagramWithZeroPrices(current_bidder_diagrams_min_[i]);
    }
    if(do_sad) {
      D1_sad = diagramWithZeroPrices(current_bidder_diagrams_saddle_[i]);
    }
    if(do_max) {
      D1_max = diagramWithZeroPrices(current_bidder_diagrams_max_[i]);
    }

    for(int c = 0; c < k_; ++c) {
      if(inv_clustering_[i] == -1) {
        // If not yet assigned, assign it first to a random cluster

        if(deterministic_) {
          inv_clustering_[i] = i % k_;
        } else {
          std::cout << " - ASSIGNED TO A RANDOM CLUSTER " << '\n';
          inv_clustering_[i] = rand() % (k_);
        }

        r_[i] = true;
        if(do_min) {
          centroids_with_price_min_[i]
            = centroidWithZeroPrices(centroids_min_[inv_clustering_[i]]);
        }
        if(do_sad) {
          centroids_with_price_saddle_[i]
            = centroidWithZeroPrices(centroids_saddle_[inv_clustering_[i]]);
        }
        if(do_max) {
          centroids_with_price_max_[i]
            = centroidWithZeroPrices(centroids_max_[inv_clustering_[i]]);
        }
      }

      if(c != inv_clustering_[i] && u_[i] > l_[i][c]
         && u_[i] > 0.5 * centroidsDistanceMatrix_[inv_clustering_[i]][c]) {
        // Step 3a, If necessary, recompute the distance to centroid
        if(r_[i]) {
          double distance = 0;
          GoodDiagram centroid_min, centroid_sad, centroid_max;
          if(do_min) {
            centroid_min
              = centroidWithZeroPrices(centroids_min_[inv_clustering_[i]]);
            distance += computeDistance(D1_min, centroid_min, 0.01);
          }
          if(do_sad) {
            centroid_sad
              = centroidWithZeroPrices(centroids_saddle_[inv_clustering_[i]]);
            distance += computeDistance(D1_sad, centroid_sad, 0.01);
          }
          if(do_max) {
            centroid_max
              = centroidWithZeroPrices(centroids_max_[inv_clustering_[i]]);
            distance += computeDistance(D1_max, centroid_max, 0.01);
          }
          r_[i] = false;
          u_[i] = distance;
          l_[i][inv_clustering_[i]] = distance;
        }
        // Step 3b, check if still potential change of clusters
        if((n_iterations_ > 2 || n_iterations_ < 1)
           && (u_[i] > l_[i][c]
               || u_[i]
                    > 0.5 * centroidsDistanceMatrix_[inv_clustering_[i]][c])) {
          BidderDiagram diagram_min, diagram_sad, diagram_max;
          GoodDiagram centroid_min, centroid_sad, centroid_max;
          double distance = 0;

          if(do_min) {
            centroid_min = centroidWithZeroPrices(centroids_min_[c]);
            diagram_min
              = diagramWithZeroPrices(current_bidder_diagrams_min_[i]);
            distance += computeDistance(diagram_min, centroid_min, 0.01);
          }
          if(do_sad) {
            centroid_sad = centroidWithZeroPrices(centroids_saddle_[c]);
            diagram_sad
              = diagramWithZeroPrices(current_bidder_diagrams_saddle_[i]);
            distance += computeDistance(diagram_sad, centroid_sad, 0.01);
          }
          if(do_max) {
            centroid_max = centroidWithZeroPrices(centroids_max_[c]);
            diagram_max
              = diagramWithZeroPrices(current_bidder_diagrams_max_[i]);
            distance += computeDistance(diagram_max, centroid_max, 0.01);
          }
          l_[i][c] = distance;
          // TODO Prices are lost here... If distance<self.u[i], we should keep
          // the prices
          if(distance < u_[i]) {
            // Changing cluster
            resetDosToOriginalValues();
            barycenter_inputs_reset_flag = true;
            u_[i] = distance;
            inv_clustering_[i] = c;

            if(do_min) {
              centroids_with_price_min_[i]
                = centroidWithZeroPrices(centroids_min_[c]);
            }
            if(do_sad) {
              centroids_with_price_saddle_[i]
                = centroidWithZeroPrices(centroids_saddle_[c]);
            }
            if(do_max) {
              centroids_with_price_max_[i]
                = centroidWithZeroPrices(centroids_max_[c]);
            }
          }
        }
      }
    }
  }
  invertInverseClusters();
  for(int c = 0; c < k_; ++c) {
    if(clustering_[c].empty()) {
      std::cout << "Adding artificial centroid because a cluster was empty"
                << std::endl;
      bool idx_acceptable = false;
      int idx = 0;

      std::vector<double> copy_of_u(u_.size());
      copy_of_u = u_;
      while(!idx_acceptable) {
        auto argMax = std::max_element(copy_of_u.begin(), copy_of_u.end());
        idx = std::distance(copy_of_u.begin(), argMax);
        if(inv_clustering_[idx] < k_ && inv_clustering_[idx] >= 0
           && clustering_[inv_clustering_[idx]].size() > 1) {
          idx_acceptable = true;
          int const cluster_removal = inv_clustering_[idx];
          // Removing the index to remove
          clustering_[cluster_removal].erase(
            std::remove(clustering_[cluster_removal].begin(),
                        clustering_[cluster_removal].end(), idx),
            clustering_[cluster_removal].end());
        } else {
          if(copy_of_u.size() > (size_t)idx) {
            copy_of_u.erase(argMax);
          } else {
            idx_acceptable = true;
            int cluster_max = 0;
            if(!clustering_[cluster_max].empty()) {
              idx = clustering_[cluster_max][0];
            }
            for(int i_test = 1; i_test < k_; i_test++) {
              if(clustering_[i_test].size() > clustering_[cluster_max].size()) {
                cluster_max = i_test;
                idx = clustering_[cluster_max][0];
              }
            }
            int const cluster_removal = inv_clustering_[idx];
            clustering_[cluster_removal].erase(
              std::remove(clustering_[cluster_removal].begin(),
                          clustering_[cluster_removal].end(), idx),
              clustering_[cluster_removal].end());
          }
        }
      }

      clustering_[c].emplace_back(idx);
      inv_clustering_[idx] = c;

      if(do_min) {
        centroids_min_[c]
          = diagramToCentroid(current_bidder_diagrams_min_[idx]);
        centroids_with_price_min_[idx]
          = centroidWithZeroPrices(centroids_min_[c]);
      }
      if(do_sad) {
        centroids_saddle_[c]
          = diagramToCentroid(current_bidder_diagrams_saddle_[idx]);
        centroids_with_price_saddle_[idx]
          = centroidWithZeroPrices(centroids_saddle_[c]);
      }
      if(do_max) {
        centroids_max_[c]
          = diagramToCentroid(current_bidder_diagrams_max_[idx]);
        centroids_with_price_max_[idx]
          = centroidWithZeroPrices(centroids_max_[c]);
      }
      resetDosToOriginalValues();
      barycenter_inputs_reset_flag = true;
    }
  }
}

std::vector<double> ttk::PDClustering::updateCentroidsPosition(
  std::vector<std::vector<double>> *min_price,
  std::vector<std::vector<double>> *min_diag_price,
  std::vector<std::vector<std::vector<std::vector<MatchingType>>>>
    &all_matchings_per_type_and_cluster,
  int only_matchings) {
  barycenter_inputs_reset_flag = true;
  std::vector<double> max_shift_vector(3);
  max_shift_vector[0] = 0;
  max_shift_vector[1] = 0;
  max_shift_vector[2] = 0;
  double max_shift_c_min = 0;
  double max_shift_c_sad = 0;
  double max_shift_c_max = 0;
  double max_wasserstein_shift = 0;
  bool precision_min = true;
  bool precision_sad = true;
  bool precision_max = true;
  cost_ = 0;
  double const sq_dist_min = cost_min_;
  double const sq_dist_sad = cost_sad_;
  double const sq_dist_max = cost_max_;
  if(do_min_) {
    cost_min_ = 0;
  }
  if(do_sad_) {
    cost_sad_ = 0;
  }
  if(do_max_) {
    cost_max_ = 0;
  }

  for(int c = 0; c < k_; ++c) {
    if(!clustering_[c].empty()) {
      std::vector<GoodDiagram> centroids_with_price_min,
        centroids_with_price_sad, centroids_with_price_max;
      int count = 0;
      for(int const idx : clustering_[c]) {
        // Timer time_first_thing;
        // Find the position of diagrams[idx] in old cluster c
        std::vector<int>::iterator const i = std::find(
          old_clustering_[c].begin(), old_clustering_[c].end(), idx);
        int const pos = (i == old_clustering_[c].end())
                          ? -1
                          : std::distance(old_clustering_[c].begin(), i);
        if(pos >= 0) {
          // Diagram was already linked to this centroid before
          if(do_min_) {
            centroids_with_price_min.emplace_back(
              centroids_with_price_min_[idx]);
          }
          if(do_sad_) {
            centroids_with_price_sad.emplace_back(
              centroids_with_price_saddle_[idx]);
          }
          if(do_max_) {
            centroids_with_price_max.emplace_back(
              centroids_with_price_max_[idx]);
          }
        } else {
          int number_of_points_min = 0;
          int number_of_points_max = 0;
          int number_of_points_sad = 0;

          // Otherwise, centroid is given 0 prices and the diagram is given 0
          // diagonal-prices
          if(do_min_) {
            centroids_with_price_min.emplace_back(
              centroidWithZeroPrices(centroids_min_[c]));
            current_bidder_diagrams_min_[idx]
              = diagramWithZeroPrices(current_bidder_diagrams_min_[idx]);
            number_of_points_min += centroids_with_price_min_[idx].size()
                                    + current_bidder_diagrams_min_[idx].size();
          }
          if(do_sad_) {
            centroids_with_price_sad.emplace_back(
              centroidWithZeroPrices(centroids_saddle_[c]));
            current_bidder_diagrams_saddle_[idx]
              = diagramWithZeroPrices(current_bidder_diagrams_saddle_[idx]);
            number_of_points_sad
              += centroids_with_price_saddle_[idx].size()
                 + current_bidder_diagrams_saddle_[idx].size();
          }
          if(do_max_) {
            centroids_with_price_max.emplace_back(
              centroidWithZeroPrices(centroids_max_[c]));
            current_bidder_diagrams_max_[idx]
              = diagramWithZeroPrices(current_bidder_diagrams_max_[idx]);
            number_of_points_max += centroids_with_price_max_[idx].size()
                                    + current_bidder_diagrams_max_[idx].size();
          }

          if(n_iterations_ > 1) {
            // If diagram new to cluster and we're not at first iteration,
            // precompute prices for the objects via compute_distance()
            // number_of_points /= (int)do_min_ + (int)do_sad_ + (int)do_max_;
            // double d_estimated = pow(cost_ / numberOfInputs_, 1. /
            // wasserstein_) + 1e-7; We use pointer in the auction in order to
            // keep the prices at the end

            if(do_min_) {
              // double estimated_delta_lim = number_of_points_min *
              // epsilon_[0] / (2*sq_dist_min) ;
              double estimated_delta_lim
                = 1.
                    / sqrt(1 - number_of_points_min * epsilon_[0] / sq_dist_min)
                  - 1;
              if(estimated_delta_lim > 1) {
                estimated_delta_lim = 1;
              }
              computeDistance(&(current_bidder_diagrams_min_[idx]),
                              &(centroids_with_price_min[count]),
                              estimated_delta_lim);
            }
            if(do_sad_) {
              // double estimated_delta_lim = number_of_points_sad *
              // epsilon_[1] / (2*sq_dist_sad);
              double estimated_delta_lim
                = 1.
                    / sqrt(1 - number_of_points_sad * epsilon_[1] / sq_dist_sad)
                  - 1;
              if(estimated_delta_lim > 1) {
                estimated_delta_lim = 1;
              }
              computeDistance(&(current_bidder_diagrams_saddle_[idx]),
                              &(centroids_with_price_sad[count]),
                              estimated_delta_lim);
            }
            if(do_max_) {
              // double estimated_delta_lim = number_of_points_max *
              // epsilon_[2] / (2*sq_dist_max);
              double estimated_delta_lim
                = 1.
                    / sqrt(1 - number_of_points_max * epsilon_[2] / sq_dist_max)
                  - 1;
              if(estimated_delta_lim > 1) {
                estimated_delta_lim = 1;
              }
              computeDistance(&(current_bidder_diagrams_max_[idx]),
                              &(centroids_with_price_max[count]),
                              estimated_delta_lim);
            }
          }
        }
        count++;
      }
      double total_cost = 0;
      double wasserstein_shift = 0;

      using KDTreePair = PDBarycenter::KDTreePair;

      if(do_min_) {
        std::vector<std::vector<MatchingType>> all_matchings;
        std::vector<int> sizes;
        Timer const time_preprocess_bary;
        std::vector<BidderDiagram> diagrams_c_min;
        if(barycenter_inputs_reset_flag) {
          for(int const idx : clustering_[c]) {
            diagrams_c_min.emplace_back(current_bidder_diagrams_min_[idx]);
          }
          sizes.resize(diagrams_c_min.size());
          for(size_t i = 0; i < diagrams_c_min.size(); i++) {
            sizes[i] = diagrams_c_min[i].size();
          }
          barycenter_computer_min_[c].setNumberOfInputs(diagrams_c_min.size());
          barycenter_computer_min_[c].setCurrentBidders(diagrams_c_min);

          std::vector<GoodDiagram> barycenter_goods(clustering_[c].size());
          for(size_t i_diagram = 0; i_diagram < clustering_[c].size();
              i_diagram++) {
            barycenter_goods[i_diagram]
              = centroids_with_price_min_[clustering_[c][i_diagram]];
          }
          barycenter_computer_min_[c].setCurrentBarycenter(barycenter_goods);
          all_matchings.resize(diagrams_c_min.size());
        } else {
          sizes.resize(barycenter_computer_min_[c].getCurrentBidders().size());
          for(size_t i = 0;
              i < barycenter_computer_min_[c].getCurrentBidders().size(); i++) {
            sizes[i]
              = barycenter_computer_min_[c].getCurrentBidders().at(i).size();
          }
          diagrams_c_min.resize(
            barycenter_computer_min_[c].getCurrentBidders().size());
          all_matchings.resize(
            barycenter_computer_min_[c].getCurrentBidders().size());
        }

        KDTreePair pair;
        bool use_kdt = false;
        if(!barycenter_computer_min_[c].getCurrentBarycenter()[0].empty()) {
          pair = barycenter_computer_min_[c].getKDTree();
          use_kdt = true;
        }

        barycenter_computer_min_[c].runMatching(
          &total_cost, epsilon_[0], sizes, *pair.first, pair.second,
          &(min_diag_price->at(0)), &(min_price->at(0)), &(all_matchings),
          use_kdt, only_matchings);
        for(size_t ii = 0; ii < all_matchings.size(); ii++) {
          all_matchings_per_type_and_cluster[c][0][ii].resize(
            all_matchings[ii].size());
          all_matchings_per_type_and_cluster[c][0][ii] = all_matchings[ii];
        }
        for(int ii = all_matchings.size(); ii < numberOfInputs_; ii++) {
          all_matchings_per_type_and_cluster[c][0][ii].resize(0);
        }

        precision_min
          = barycenter_computer_min_[c].isPrecisionObjectiveMet(deltaLim_, 0);
        cost_min_ += total_cost;
        Timer const time_update;
        if(!only_matchings) {
          max_shift_c_min
            = barycenter_computer_min_[c].updateBarycenter(all_matchings);
        }

        if(max_shift_c_min > max_shift_vector[0]) {
          max_shift_vector[0] = max_shift_c_min;
        }

        // Now that barycenters and diagrams are updated in
        // PDBarycenter class, we import the results here.
        diagrams_c_min = barycenter_computer_min_[c].getCurrentBidders();
        centroids_with_price_min
          = barycenter_computer_min_[c].getCurrentBarycenter();
        int i = 0;
        for(int const idx : clustering_[c]) {
          current_bidder_diagrams_min_[idx] = diagrams_c_min[i];
          centroids_with_price_min_[idx] = centroids_with_price_min[i];
          i++;
        }

        GoodDiagram const old_centroid = centroids_min_[c];
        centroids_min_[c] = centroidWithZeroPrices(
          centroids_with_price_min_[clustering_[c][0]]);

        if(use_accelerated_) {
          wasserstein_shift
            += computeDistance(old_centroid, centroids_min_[c], 0.01);
        }
      }

      if(do_sad_) {
        std::vector<std::vector<MatchingType>> all_matchings;
        total_cost = 0;
        std::vector<int> sizes;

        std::vector<BidderDiagram> diagrams_c_min;
        if(barycenter_inputs_reset_flag) {
          for(int const idx : clustering_[c]) {
            diagrams_c_min.emplace_back(current_bidder_diagrams_saddle_[idx]);
          }
          sizes.resize(diagrams_c_min.size());
          for(size_t i = 0; i < diagrams_c_min.size(); i++) {
            sizes[i] = diagrams_c_min[i].size();
          }
          barycenter_computer_sad_[c].setNumberOfInputs(diagrams_c_min.size());
          barycenter_computer_sad_[c].setCurrentBidders(diagrams_c_min);
          std::vector<GoodDiagram> barycenter_goods(clustering_[c].size());
          for(size_t i_diagram = 0; i_diagram < clustering_[c].size();
              i_diagram++) {
            barycenter_goods[i_diagram]
              = centroids_with_price_saddle_[clustering_[c][i_diagram]];
          }
          barycenter_computer_sad_[c].setCurrentBarycenter(barycenter_goods);
          all_matchings.resize(diagrams_c_min.size());
        } else {
          sizes.resize(barycenter_computer_sad_[c].getCurrentBidders().size());
          for(size_t i = 0;
              i < barycenter_computer_sad_[c].getCurrentBidders().size(); i++) {
            sizes[i]
              = barycenter_computer_sad_[c].getCurrentBidders().at(i).size();
          }
          all_matchings.resize(
            barycenter_computer_sad_[c].getCurrentBidders().size());
          diagrams_c_min.resize(
            barycenter_computer_sad_[c].getCurrentBidders().size());
        }

        KDTreePair pair;
        bool use_kdt = false;
        if(!barycenter_computer_sad_[c].getCurrentBarycenter()[0].empty()) {
          pair = barycenter_computer_sad_[c].getKDTree();
          use_kdt = true;
        }

        barycenter_computer_sad_[c].runMatching(
          &total_cost, epsilon_[1], sizes, *pair.first, pair.second,
          &(min_diag_price->at(1)), &(min_price->at(1)), &(all_matchings),
          use_kdt, only_matchings);
        for(size_t ii = 0; ii < all_matchings.size(); ii++) {
          all_matchings_per_type_and_cluster[c][1][ii].resize(
            all_matchings[ii].size());
          all_matchings_per_type_and_cluster[c][1][ii] = all_matchings[ii];
        }
        for(int ii = all_matchings.size(); ii < numberOfInputs_; ii++) {
          all_matchings_per_type_and_cluster[c][1][ii].resize(0);
        }

        precision_sad
          = barycenter_computer_sad_[c].isPrecisionObjectiveMet(deltaLim_, 0);
        if(!only_matchings) {
          max_shift_c_sad
            = barycenter_computer_sad_[c].updateBarycenter(all_matchings);
        }

        cost_sad_ += total_cost;

        if(max_shift_c_sad > max_shift_vector[1]) {
          max_shift_vector[1] = max_shift_c_sad;
        }

        // Now that barycenters and diagrams are updated in PDBarycenter class,
        // we import the results here.
        diagrams_c_min = barycenter_computer_sad_[c].getCurrentBidders();
        centroids_with_price_sad
          = barycenter_computer_sad_[c].getCurrentBarycenter();
        int i = 0;
        for(int const idx : clustering_[c]) {
          current_bidder_diagrams_saddle_[idx] = diagrams_c_min[i];
          centroids_with_price_saddle_[idx] = centroids_with_price_sad[i];
          i++;
        }
        GoodDiagram const old_centroid = centroids_saddle_[c];
        centroids_saddle_[c] = centroidWithZeroPrices(
          centroids_with_price_saddle_[clustering_[c][0]]);
        if(use_accelerated_)
          wasserstein_shift
            += computeDistance(old_centroid, centroids_saddle_[c], 0.01);
      }

      if(do_max_) {
        std::vector<std::vector<MatchingType>> all_matchings;
        Timer const time_preprocess_bary;
        total_cost = 0;
        std::vector<int> sizes;
        std::vector<BidderDiagram> diagrams_c_min;
        if(barycenter_inputs_reset_flag) {
          for(int const idx : clustering_[c]) {
            diagrams_c_min.emplace_back(current_bidder_diagrams_max_[idx]);
          }
          sizes.resize(diagrams_c_min.size());
          for(size_t i = 0; i < diagrams_c_min.size(); i++) {
            sizes[i] = diagrams_c_min[i].size();
          }
          barycenter_computer_max_[c].setNumberOfInputs(diagrams_c_min.size());
          barycenter_computer_max_[c].setCurrentBidders(diagrams_c_min);
          std::vector<GoodDiagram> barycenter_goods(clustering_[c].size());
          for(size_t i_diagram = 0; i_diagram < clustering_[c].size();
              i_diagram++) {
            barycenter_goods[i_diagram]
              = centroids_with_price_max_[clustering_[c][i_diagram]];
          }
          barycenter_computer_max_[c].setCurrentBarycenter(barycenter_goods);
          all_matchings.resize(diagrams_c_min.size());
        } else {
          sizes.resize(barycenter_computer_max_[c].getCurrentBidders().size());
          for(size_t i = 0;
              i < barycenter_computer_max_[c].getCurrentBidders().size(); i++) {
            sizes[i]
              = barycenter_computer_max_[c].getCurrentBidders().at(i).size();
          }

          diagrams_c_min.resize(
            barycenter_computer_max_[c].getCurrentBidders().size());
          all_matchings.resize(
            barycenter_computer_max_[c].getCurrentBidders().size());
        }

        KDTreePair pair;
        bool use_kdt = false;
        if(!barycenter_computer_max_[c].getCurrentBarycenter()[0].empty()) {
          pair = barycenter_computer_max_[c].getKDTree();
          use_kdt = true;
        }

        barycenter_computer_max_[c].runMatching(
          &total_cost, epsilon_[2], sizes, *pair.first, pair.second,
          &(min_diag_price->at(2)), &(min_price->at(2)), &(all_matchings),
          use_kdt, only_matchings);
        for(size_t ii = 0; ii < all_matchings.size(); ii++) {
          all_matchings_per_type_and_cluster[c][2][ii].resize(
            all_matchings[ii].size());
          all_matchings_per_type_and_cluster[c][2][ii] = all_matchings[ii];
        }
        for(int ii = all_matchings.size(); ii < numberOfInputs_; ii++) {
          all_matchings_per_type_and_cluster[c][2][ii].resize(0);
        }
        precision_max
          = barycenter_computer_max_[c].isPrecisionObjectiveMet(deltaLim_, 0);

        cost_max_ += total_cost;
        Timer const time_update;
        if(!only_matchings) {
          max_shift_c_max
            = barycenter_computer_max_[c].updateBarycenter(all_matchings);
        }
        if(max_shift_c_max > max_shift_vector[2]) {
          max_shift_vector[2] = max_shift_c_max;
        }

        // Now that barycenters and diagrams are updated in PDBarycenter class,
        // we import the results here.
        diagrams_c_min = barycenter_computer_max_[c].getCurrentBidders();
        centroids_with_price_max
          = barycenter_computer_max_[c].getCurrentBarycenter();
        int i = 0;
        for(int const idx : clustering_[c]) {
          current_bidder_diagrams_max_[idx] = diagrams_c_min[i];
          centroids_with_price_max_[idx] = centroids_with_price_max[i];
          i++;
        }
        GoodDiagram const old_centroid = centroids_max_[c];
        centroids_max_[c] = centroidWithZeroPrices(
          centroids_with_price_max_[clustering_[c][0]]);
        if(use_accelerated_) {
          wasserstein_shift
            += computeDistance(old_centroid, centroids_max_[c], 0.01);
        }
      }

      cost_ = cost_min_ + cost_sad_ + cost_max_;
      if(wasserstein_shift > max_wasserstein_shift) {
        max_wasserstein_shift = wasserstein_shift;
      }
      if(use_accelerated_) {
        for(int i = 0; i < numberOfInputs_; ++i) {
          // Step 5 of Accelerated KMeans: Update the lower bound on distance
          // thanks to the triangular inequality
          l_[i][c] = Geometry::pow(
            Geometry::pow(l_[i][c], 1. / wasserstein_)
              - Geometry::pow(wasserstein_shift, 1. / wasserstein_),
            wasserstein_);
          if(l_[i][c] < 0) {
            l_[i][c] = 0;
          }
        }
        for(int const idx : clustering_[c]) {
          // Step 6, update the upper bound on the distance to the centroid
          // thanks to the triangle inequality
          u_[idx] = Geometry::pow(
            Geometry::pow(u_[idx], 1. / wasserstein_)
              + Geometry::pow(wasserstein_shift, 1. / wasserstein_),
            wasserstein_);
          r_[idx] = true;
        }
      }
    }
  }
  // Normally return max_shift, but it seems there is a bug
  // yielding max_shift > 100 * max_wasserstein_shift
  // which should logically not really happen...
  // This is supposed to be only a temporary patch...
  precision_min_ = precision_min;
  precision_sad_ = precision_sad;
  precision_max_ = precision_max;
  precision_criterion_ = precision_min && precision_sad && precision_max;
  barycenter_inputs_reset_flag = false;
  return max_shift_vector; // std::min(max_shift, max_wasserstein_shift);
}

void ttk::PDClustering::setBidderDiagrams() {
  for(int i = 0; i < numberOfInputs_; i++) {
    if(do_min_) {
      DiagramType *CTDiagram = &((*inputDiagramsMin_)[i]);
      BidderDiagram bidders;
      for(size_t j = 0; j < CTDiagram->size(); j++) {
        // Add bidder to bidders
        Bidder b((*CTDiagram)[j], j, lambda_);

        b.setPositionInAuction(bidders.size());
        bidders.emplace_back(b);
        if(b.isDiagonal() || b.x_ == b.y_) {
          this->printMsg("Diagonal point in diagram", debug::Priority::DETAIL);
        }
      }
      bidder_diagrams_min_.emplace_back(bidders);
      current_bidder_diagrams_min_.emplace_back();
      centroids_with_price_min_.emplace_back();
      std::vector<int> ids(bidders.size());
      for(size_t j = 0; j < ids.size(); j++) {
        ids[j] = -1;
      }
      current_bidder_ids_min_.emplace_back(ids);
    }

    if(do_sad_) {
      DiagramType *CTDiagram = &((*inputDiagramsSaddle_)[i]);

      BidderDiagram bidders;
      for(size_t j = 0; j < CTDiagram->size(); j++) {
        // Add bidder to bidders
        Bidder b((*CTDiagram)[j], j, lambda_);

        b.setPositionInAuction(bidders.size());
        bidders.emplace_back(b);
        if(b.isDiagonal() || b.x_ == b.y_) {
          this->printMsg("Diagonal point in diagram", debug::Priority::DETAIL);
        }
      }
      bidder_diagrams_saddle_.emplace_back(bidders);
      current_bidder_diagrams_saddle_.emplace_back();
      centroids_with_price_saddle_.emplace_back();
      std::vector<int> ids(bidders.size());
      for(size_t j = 0; j < ids.size(); j++) {
        ids[j] = -1;
      }
      current_bidder_ids_sad_.emplace_back(ids);
    }

    if(do_max_) {
      DiagramType *CTDiagram = &((*inputDiagramsMax_)[i]);

      BidderDiagram bidders;
      for(size_t j = 0; j < CTDiagram->size(); j++) {
        // Add bidder to bidders
        Bidder b((*CTDiagram)[j], j, lambda_);

        b.setPositionInAuction(bidders.size());
        bidders.emplace_back(b);
        if(b.isDiagonal() || b.x_ == b.y_) {
          this->printMsg("Diagonal point in diagram", debug::Priority::DETAIL);
        }
      }
      bidder_diagrams_max_.emplace_back(bidders);
      current_bidder_diagrams_max_.emplace_back();
      centroids_with_price_max_.emplace_back();
      std::vector<int> ids(bidders.size());
      for(size_t j = 0; j < ids.size(); j++) {
        ids[j] = -1;
      }
      current_bidder_ids_max_.emplace_back(ids);
    }
  }
}

std::vector<double> ttk::PDClustering::enrichCurrentBidderDiagrams(
  std::vector<double> &previous_min_persistence,
  std::vector<double> &min_persistence,
  std::vector<std::vector<double>> &initial_diagonal_prices,
  std::vector<std::vector<double>> &initial_off_diagonal_prices,
  std::vector<int> &min_points_to_add,
  bool add_points_to_barycenter,
  bool first_enrichment) {

  std::vector<double> new_min_persistence = min_persistence;

  if(!do_min_) {
    new_min_persistence[0] = previous_min_persistence[0];
  }
  if(!do_sad_) {
    new_min_persistence[1] = previous_min_persistence[1];
  }
  if(!do_max_) {
    new_min_persistence[2] = previous_min_persistence[2];
  }

  // 1. Get size of the largest current diagram, deduce the maximal number of
  // points to append
  size_t max_diagram_size_min{};
  size_t max_diagram_size_sad{};
  size_t max_diagram_size_max{};
  if(do_min_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      if(current_bidder_diagrams_min_[i].size() > max_diagram_size_min) {
        max_diagram_size_min = current_bidder_diagrams_min_[i].size();
      }
    }
  }
  if(do_sad_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      if(current_bidder_diagrams_saddle_[i].size() > max_diagram_size_sad) {
        max_diagram_size_sad = current_bidder_diagrams_saddle_[i].size();
      }
    }
  }
  if(do_max_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      if(current_bidder_diagrams_max_[i].size() > max_diagram_size_max) {
        max_diagram_size_max = current_bidder_diagrams_max_[i].size();
      }
    }
  }
  int const max_points_to_add_min
    = std::max(min_points_to_add[0],
               min_points_to_add[0] + (int)(max_diagram_size_min / 10));
  int const max_points_to_add_sad
    = std::max(min_points_to_add[1],
               min_points_to_add[1] + (int)(max_diagram_size_sad / 10));
  int const max_points_to_add_max
    = std::max(min_points_to_add[2],
               min_points_to_add[2] + (int)(max_diagram_size_max / 10));

  // 2. Get which points can be added, deduce the new minimal persistence
  std::vector<std::vector<int>> candidates_to_be_added_min(numberOfInputs_);
  std::vector<std::vector<int>> candidates_to_be_added_sad(numberOfInputs_);
  std::vector<std::vector<int>> candidates_to_be_added_max(numberOfInputs_);
  std::vector<std::vector<int>> idx_min(numberOfInputs_);
  std::vector<std::vector<int>> idx_sad(numberOfInputs_);
  std::vector<std::vector<int>> idx_max(numberOfInputs_);

  if(do_min_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      std::vector<double> persistences;
      for(size_t j = 0; j < bidder_diagrams_min_[i].size(); j++) {
        Bidder const b = bidder_diagrams_min_[i].at(j);
        double const persistence = b.getPersistence();
        if(persistence >= min_persistence[0]
           && persistence <= previous_min_persistence[0]) {
          candidates_to_be_added_min[i].emplace_back(j);
          idx_min[i].emplace_back(idx_min[i].size());
          persistences.emplace_back(persistence);
        }
      }
      sort(
        idx_min[i].begin(), idx_min[i].end(), [&persistences](int &a, int &b) {
          return ((persistences[a] > persistences[b])
                  || ((persistences[a] == persistences[b]) && (a > b)));
        });
      int const size = candidates_to_be_added_min[i].size();
      if(size >= max_points_to_add_min) {
        double const last_persistence_added_min
          = persistences[idx_min[i][max_points_to_add_min - 1]];
        if(first_enrichment) { // a minima min_point_to_add (=max_point_to_add)
                               // added per diagram
          if(i == 0) {
            new_min_persistence[0] = last_persistence_added_min;
          } else {
            if(last_persistence_added_min < new_min_persistence[0])
              new_min_persistence[0] = last_persistence_added_min;
          }
        } else { // a maxima max_point_to_add added per diagram
          if(last_persistence_added_min > new_min_persistence[0]) {
            new_min_persistence[0] = last_persistence_added_min;
          }
        }
      }
    }
  }

  if(do_sad_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      std::vector<double> persistences;
      for(size_t j = 0; j < bidder_diagrams_saddle_[i].size(); j++) {
        Bidder const b = bidder_diagrams_saddle_[i].at(j);
        double const persistence = b.getPersistence();
        if(persistence >= min_persistence[1]
           && persistence <= previous_min_persistence[1]) {
          candidates_to_be_added_sad[i].emplace_back(j);
          idx_sad[i].emplace_back(idx_sad[i].size());
          persistences.emplace_back(persistence);
        }
      }
      sort(
        idx_sad[i].begin(), idx_sad[i].end(), [&persistences](int &a, int &b) {
          return ((persistences[a] > persistences[b])
                  || ((persistences[a] == persistences[b]) && (a > b)));
        });
      int const size = candidates_to_be_added_sad[i].size();
      if(size >= max_points_to_add_sad) {
        double const last_persistence_added_sad
          = persistences[idx_sad[i][max_points_to_add_sad - 1]];
        if(first_enrichment) { // a minima min_point_to_add (=max_point_to_add)
                               // added per diagram
          if(i == 0) {
            new_min_persistence[1] = last_persistence_added_sad;
          } else {
            if(last_persistence_added_sad < new_min_persistence[1])
              new_min_persistence[1] = last_persistence_added_sad;
          }
        } else { // a maxima max_point_to_add added per diagram
          if(last_persistence_added_sad > new_min_persistence[1]) {
            new_min_persistence[1] = last_persistence_added_sad;
          }
        }
      }
    }
  }
  if(do_max_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      std::vector<double> persistences;
      for(size_t j = 0; j < bidder_diagrams_max_[i].size(); j++) {
        Bidder const b = bidder_diagrams_max_[i].at(j);
        double const persistence = b.getPersistence();
        if(persistence >= min_persistence[2]
           && persistence <= previous_min_persistence[2]) {
          candidates_to_be_added_max[i].emplace_back(j);
          idx_max[i].emplace_back(idx_max[i].size());
          persistences.emplace_back(persistence);
        }
      }
      sort(
        idx_max[i].begin(), idx_max[i].end(), [&persistences](int &a, int &b) {
          return ((persistences[a] > persistences[b])
                  || ((persistences[a] == persistences[b]) && (a > b)));
        });
      int const size = candidates_to_be_added_max[i].size();
      if(size >= max_points_to_add_max) {
        double const last_persistence_added_max
          = persistences[idx_max[i][max_points_to_add_max - 1]];
        if(first_enrichment) { // a minima min_point_to_add (=max_point_to_add)
                               // added per diagram
          if(i == 0) {
            new_min_persistence[2] = last_persistence_added_max;
          } else {
            if(last_persistence_added_max < new_min_persistence[2])
              new_min_persistence[2] = last_persistence_added_max;
          }
        } else { // a maxima max_point_to_add added per diagram
          if(last_persistence_added_max > new_min_persistence[2]) {
            new_min_persistence[2] = last_persistence_added_max;
          }
        }
      }
    }
  }

  // 3. Add the points to the current diagrams
  if(do_min_) {
    int compteur_for_adding_points = 0;
    for(int i = 0; i < numberOfInputs_; i++) {
      int const size = candidates_to_be_added_min[i].size();
      for(int j = 0; j < std::min(max_points_to_add_min, size); j++) {
        Bidder b = bidder_diagrams_min_[i].at(
          candidates_to_be_added_min[i][idx_min[i][j]]);
        double const persistence = b.getPersistence();
        if(persistence >= new_min_persistence[0]) {
          b.id_ = current_bidder_diagrams_min_[i].size();
          b.setPositionInAuction(current_bidder_diagrams_min_[i].size());
          b.setDiagonalPrice(initial_diagonal_prices[0][i]);
          current_bidder_diagrams_min_[i].emplace_back(b);
          current_bidder_ids_min_[i]
                                 [candidates_to_be_added_min[i][idx_min[i][j]]]
            = current_bidder_diagrams_min_[i].size() - 1;

          if(use_accelerated_ && n_iterations_ > 0) {
            for(int c = 0; c < k_; ++c) {
              // Step 5 of Accelerated KMeans: Update the lower bound on
              // distance thanks to the triangular inequality
              l_[i][c]
                = Geometry::pow(Geometry::pow(l_[i][c], 1. / wasserstein_)
                                  - persistence / sqrt(2),
                                wasserstein_);
              if(l_[i][c] < 0) {
                l_[i][c] = 0;
              }
            }
            // Step 6, update the upper bound on the distance to the centroid
            // thanks to the triangle inequality
            u_[i] = Geometry::pow(
              Geometry::pow(u_[i], 1. / wasserstein_) + persistence / sqrt(2),
              wasserstein_);
            r_[i] = true;
          }
          int const to_be_added_to_barycenter
            = deterministic_ ? compteur_for_adding_points % numberOfInputs_
                             : rand() % numberOfInputs_;
          if(to_be_added_to_barycenter == 0 && add_points_to_barycenter) {
            for(int k = 0; k < numberOfInputs_; k++) {
              if(inv_clustering_[i] == inv_clustering_[k]) {
                Good g = Good(
                  b.x_, b.y_, false, centroids_with_price_min_[k].size());
                g.setPrice(initial_off_diagonal_prices[0][k]);
                g.SetCriticalCoordinates(
                  b.coords_[0], b.coords_[1], b.coords_[2]);
                centroids_with_price_min_[k].emplace_back(g);
              }
            }
            Good g = Good(
              b.x_, b.y_, false, centroids_min_[inv_clustering_[i]].size());
            g.SetCriticalCoordinates(b.coords_[0], b.coords_[1], b.coords_[2]);
            centroids_min_[inv_clustering_[i]].emplace_back(g);
          }
        }
        compteur_for_adding_points++;
      }
    }
  }
  if(do_sad_) {
    int compteur_for_adding_points = 0;
    for(int i = 0; i < numberOfInputs_; i++) {
      int const size = candidates_to_be_added_sad[i].size();
      for(int j = 0; j < std::min(max_points_to_add_sad, size); j++) {
        Bidder b = bidder_diagrams_saddle_[i].at(
          candidates_to_be_added_sad[i][idx_sad[i][j]]);
        double const persistence = b.getPersistence();
        if(persistence >= new_min_persistence[1]) {
          b.id_ = current_bidder_diagrams_saddle_[i].size();
          b.setPositionInAuction(current_bidder_diagrams_saddle_[i].size());
          b.setDiagonalPrice(initial_diagonal_prices[1][i]);
          current_bidder_diagrams_saddle_[i].emplace_back(b);
          current_bidder_ids_sad_[i]
                                 [candidates_to_be_added_sad[i][idx_sad[i][j]]]
            = current_bidder_diagrams_saddle_[i].size() - 1;

          if(use_accelerated_ && n_iterations_ > 0) {
            for(int c = 0; c < k_; ++c) {
              // Step 5 of Accelerated KMeans: Update the lower bound on
              // distance thanks to the triangular inequality
              l_[i][c]
                = Geometry::pow(Geometry::pow(l_[i][c], 1. / wasserstein_)
                                  - persistence / sqrt(2),
                                wasserstein_);
              if(l_[i][c] < 0) {
                l_[i][c] = 0;
              }
            }
            // Step 6, update the upper bound on the distance to the centroid
            // thanks to the triangle inequality
            u_[i] = Geometry::pow(
              Geometry::pow(u_[i], 1. / wasserstein_) + persistence / sqrt(2),
              wasserstein_);
            r_[i] = true;
          }
          int const to_be_added_to_barycenter
            = deterministic_ ? compteur_for_adding_points % numberOfInputs_
                             : rand() % numberOfInputs_;
          if(to_be_added_to_barycenter == 0 && add_points_to_barycenter) {
            for(int k = 0; k < numberOfInputs_; k++) {
              if(inv_clustering_[i] == inv_clustering_[k]) {
                Good g = Good(
                  b.x_, b.y_, false, centroids_with_price_saddle_[k].size());
                g.setPrice(initial_off_diagonal_prices[1][k]);
                g.SetCriticalCoordinates(
                  b.coords_[0], b.coords_[1], b.coords_[2]);
                centroids_with_price_saddle_[k].emplace_back(g);
              }
            }
          }
        }
        compteur_for_adding_points++;
      }
    }
  }
  if(do_max_) {
    int compteur_for_adding_points = 0;
    for(int i = 0; i < numberOfInputs_; i++) {
      int const size = candidates_to_be_added_max[i].size();
      for(int j = 0; j < std::min(max_points_to_add_max, size); j++) {
        Bidder b = bidder_diagrams_max_[i].at(
          candidates_to_be_added_max[i][idx_max[i][j]]);
        double const persistence = b.getPersistence();
        if(persistence >= new_min_persistence[2]) {
          b.id_ = current_bidder_diagrams_max_[i].size();
          b.setPositionInAuction(current_bidder_diagrams_max_[i].size());
          b.setDiagonalPrice(initial_diagonal_prices[2][i]);
          current_bidder_diagrams_max_[i].emplace_back(b);
          current_bidder_ids_max_[i]
                                 [candidates_to_be_added_max[i][idx_max[i][j]]]
            = current_bidder_diagrams_max_[i].size() - 1;

          if(use_accelerated_ && n_iterations_ > 0) {
            for(int c = 0; c < k_; ++c) {
              // Step 5 of Accelerated KMeans: Update the lower bound on
              // distance thanks to the triangular inequality
              l_[i][c]
                = Geometry::pow(Geometry::pow(l_[i][c], 1. / wasserstein_)
                                  - persistence / sqrt(2),
                                wasserstein_);
              if(l_[i][c] < 0) {
                l_[i][c] = 0;
              }
            }
            // Step 6, update the upper bound on the distance to the centroid
            // thanks to the triangle inequality
            u_[i] = Geometry::pow(
              Geometry::pow(u_[i], 1. / wasserstein_) + persistence / sqrt(2),
              wasserstein_);
            r_[i] = true;
          }
          int const to_be_added_to_barycenter
            = deterministic_ ? compteur_for_adding_points % numberOfInputs_
                             : rand() % numberOfInputs_;
          if(to_be_added_to_barycenter == 0 && add_points_to_barycenter) {
            for(int k = 0; k < numberOfInputs_; k++) {
              if(inv_clustering_[i] == inv_clustering_[k]) {
                Good g = Good(
                  b.x_, b.y_, false, centroids_with_price_max_[k].size());
                g.setPrice(initial_off_diagonal_prices[2][k]);
                g.SetCriticalCoordinates(
                  b.coords_[0], b.coords_[1], b.coords_[2]);
                centroids_with_price_max_[k].emplace_back(g);
              }
            }
            Good g = Good(
              b.x_, b.y_, false, centroids_max_[inv_clustering_[i]].size());
            g.SetCriticalCoordinates(b.coords_[0], b.coords_[1], b.coords_[2]);
            centroids_max_[inv_clustering_[i]].emplace_back(g);
          }
        }
        compteur_for_adding_points++;
      }
    }
  }

  return new_min_persistence;
}

void ttk::PDClustering::initializeBarycenterComputers(
  std::vector<double> & /*min_persistence*/) {

  if(do_min_) {
    barycenter_computer_min_.resize(k_);
    for(int c = 0; c < k_; c++) {
      std::vector<BidderDiagram> diagrams_c;
      for(int const idx : clustering_[c]) {
        diagrams_c.emplace_back(current_bidder_diagrams_min_[idx]);
      }
      barycenter_computer_min_[c] = {};
      barycenter_computer_min_[c].setThreadNumber(threadNumber_);
      barycenter_computer_min_[c].setWasserstein(wasserstein_);
      barycenter_computer_min_[c].setDiagramType(0);
      barycenter_computer_min_[c].setUseProgressive(false);
      barycenter_computer_min_[c].setDeterministic(true);
      barycenter_computer_min_[c].setGeometricalFactor(geometrical_factor_);
      barycenter_computer_min_[c].setDebugLevel(debugLevel_);
      barycenter_computer_min_[c].setNumberOfInputs(diagrams_c.size());
      barycenter_computer_min_[c].setCurrentBidders(diagrams_c);
    }
  }
  if(do_sad_) {
    barycenter_computer_sad_.resize(k_);
    for(int c = 0; c < k_; c++) {
      std::vector<BidderDiagram> diagrams_c;
      for(int const idx : clustering_[c]) {
        diagrams_c.emplace_back(current_bidder_diagrams_saddle_[idx]);
      }
      barycenter_computer_sad_[c] = {};
      barycenter_computer_sad_[c].setThreadNumber(threadNumber_);
      barycenter_computer_sad_[c].setWasserstein(wasserstein_);
      barycenter_computer_sad_[c].setDiagramType(1);
      barycenter_computer_sad_[c].setUseProgressive(false);
      barycenter_computer_sad_[c].setDeterministic(true);
      barycenter_computer_sad_[c].setGeometricalFactor(geometrical_factor_);
      barycenter_computer_sad_[c].setDebugLevel(debugLevel_);
      barycenter_computer_sad_[c].setNumberOfInputs(diagrams_c.size());
      barycenter_computer_sad_[c].setCurrentBidders(diagrams_c);

      std::vector<GoodDiagram> barycenter_goods(clustering_[c].size());
      for(size_t i_diagram = 0; i_diagram < clustering_[c].size();
          i_diagram++) {
        barycenter_goods[i_diagram]
          = centroids_with_price_saddle_[clustering_[c][i_diagram]];
      }
      barycenter_computer_sad_[c].setCurrentBarycenter(barycenter_goods);
    }
  }
  if(do_max_) {
    barycenter_computer_max_.resize(k_);
    for(int c = 0; c < k_; c++) {
      std::vector<BidderDiagram> diagrams_c;
      for(int const idx : clustering_[c]) {
        diagrams_c.emplace_back(current_bidder_diagrams_max_[idx]);
      }
      barycenter_computer_max_[c] = {};
      barycenter_computer_max_[c].setDiagrams(inputDiagramsMax_);
      barycenter_computer_max_[c].setThreadNumber(threadNumber_);
      barycenter_computer_max_[c].setWasserstein(wasserstein_);
      barycenter_computer_max_[c].setDiagramType(2);
      barycenter_computer_max_[c].setUseProgressive(false);
      barycenter_computer_max_[c].setDeterministic(true);
      barycenter_computer_max_[c].setGeometricalFactor(geometrical_factor_);
      barycenter_computer_max_[c].setDebugLevel(debugLevel_);
      barycenter_computer_max_[c].setNumberOfInputs(diagrams_c.size());
      barycenter_computer_max_[c].setCurrentBidders(diagrams_c);

      std::vector<GoodDiagram> barycenter_goods(clustering_[c].size());
      for(size_t i_diagram = 0; i_diagram < clustering_[c].size();
          i_diagram++) {
        barycenter_goods[i_diagram]
          = centroids_with_price_max_[clustering_[c][i_diagram]];
      }
      barycenter_computer_max_[c].setCurrentBarycenter(barycenter_goods);
    }
  }
}

void ttk::PDClustering::printDistancesToFile() {
  std::ofstream ufile("u_vec.txt");
  std::ofstream lfile("l_mat.txt");
  std::ofstream approx_file("a_mat.txt");
  if(ufile.is_open() && lfile.is_open()) {
    for(int i = 0; i < numberOfInputs_; i++) {
      ufile << u_[i] << " ";
      for(int j = 0; j < k_; j++) {
        lfile << l_[i][j] << " ";
      }
      lfile << "\n";
    }
  }

  for(int c = 0; c < k_; c++) {
    for(int const i : clustering_[c]) {
      approx_file << (u_[i] + l_[i][c]) / 2 << " ";
    }
    approx_file << "\n";
  }
  lfile.close();
  ufile.close();
  approx_file.close();
}

void ttk::PDClustering::printRealDistancesToFile() {
  std::cout << "Computing real distances to every clusters" << std::endl;
  std::ofstream file("a_real_mat.txt");
  if(file.is_open()) {
    for(int c = 0; c < k_; c++) {
      for(int const i : clustering_[c]) {
        file << distanceToCentroid_[i] << " ";
      }
      file << "\n";
    }
    file.close();
  } else {
    std::cout << "file not open" << std::endl;
  }
}

void ttk::PDClustering::printPricesToFile(int iteration) {
  std::ofstream file(
    "prices_evolution.txt", std::ofstream::out | std::ofstream::app);
  if(file.is_open()) {
    file << "\nITERATION " << iteration << "\n" << std::endl;
    for(int i = 0; i < k_; i++) {
      file << "\ncentroid " << i << std::endl;

      for(size_t j = 0; j < centroids_with_price_max_[i].size(); j++) {
        Good const g = centroids_with_price_max_[i].at(j);
        file << g.getPrice() << " ";
      }
    }
  }
  file.close();
}

double ttk::PDClustering::computeRealCost() {
  double total_real_cost_min = 0;
  double total_real_cost_max = 0;
  double total_real_cost_sad = 0;
  double sq_distance;
  if(original_dos[0]) {
    for(int c = 0; c < k_; c++) {
      double real_cost_cluster = 0;
      for(int i = 0; i < numberOfInputs_; i++) {
        GoodDiagram const current_barycenter
          = centroidWithZeroPrices(centroids_min_[c]);
        sq_distance
          = computeDistance(bidder_diagrams_min_[i], current_barycenter, 0.01);
        real_cost_cluster += sq_distance;
      }
      total_real_cost_min += real_cost_cluster;
    }
  }
  if(original_dos[1]) {
    for(int c = 0; c < k_; c++) {
      double real_cost_cluster = 0;
      for(int i = 0; i < numberOfInputs_; i++) {
        GoodDiagram const current_barycenter
          = centroidWithZeroPrices(centroids_saddle_[c]);
        sq_distance = computeDistance(
          bidder_diagrams_saddle_[i], current_barycenter, 0.01);
        real_cost_cluster += sq_distance;
      }
      total_real_cost_sad += real_cost_cluster;
    }
  }
  if(original_dos[2]) {
    for(int c = 0; c < k_; c++) {
      double real_cost_cluster = 0;
      for(int i = 0; i < numberOfInputs_; i++) {
        GoodDiagram const current_barycenter
          = centroidWithZeroPrices(centroids_max_[c]);
        sq_distance
          = computeDistance(bidder_diagrams_max_[i], current_barycenter, 0.01);
        real_cost_cluster += sq_distance;
      }
      total_real_cost_max += real_cost_cluster;
    }
  }
  return total_real_cost_min + total_real_cost_sad + total_real_cost_max;
}

void ttk::PDClustering::computeBarycenterForTwoGlobal(
  std::vector<std::vector<std::vector<std::vector<MatchingType>>>>
    &all_matchings_per_type_and_cluster) {

  if(do_min_) {
    computeBarycenterForTwo(
      all_matchings_per_type_and_cluster[0][0], current_bidder_ids_min_,
      current_bidder_diagrams_min_, bidder_diagrams_min_, centroids_min_[0]);
  }
  if(do_sad_) {
    computeBarycenterForTwo(all_matchings_per_type_and_cluster[0][1],
                            current_bidder_ids_sad_,
                            current_bidder_diagrams_saddle_,
                            bidder_diagrams_saddle_, centroids_saddle_[0]);
  }
  if(do_max_) {
    computeBarycenterForTwo(
      all_matchings_per_type_and_cluster[0][2], current_bidder_ids_max_,
      current_bidder_diagrams_max_, bidder_diagrams_max_, centroids_max_[0]);
  }
}

void ttk::PDClustering::computeBarycenterForTwo(
  std::vector<std::vector<MatchingType>> &matchings,
  std::vector<std::vector<int>> &bidders_ids,
  std::vector<BidderDiagram> &current_bidder_diagrams,
  std::vector<BidderDiagram> &bidder_diagrams,
  GoodDiagram &barycenter) {

  auto &matchings0 = matchings[0];
  auto &diagram0 = bidder_diagrams[0];
  auto &current_diagram0 = current_bidder_diagrams[0];
  auto &ids0 = bidders_ids[0];

  auto &matchings1 = matchings[1];
  auto &diagram1 = bidder_diagrams[1];
  auto &current_diagram1 = current_bidder_diagrams[1];
  auto &ids1 = bidders_ids[1];

  std::vector<int> new_to_old_id(current_diagram1.size());
  // 1. Invert the current_bidder_ids_ vector
  for(size_t j = 0; j < ids1.size(); j++) {
    int const new_id = ids1[j];
    if(new_id >= 0) {
      new_to_old_id[new_id] = j;
    }
  }

  std::vector<MatchingType> matching_to_add(0);
  std::vector<MatchingType> matching_to_add2(0);

  for(size_t i = 0; i < matchings1.size(); i++) {
    MatchingType &t = matchings1[i];
    int const bidderId = std::get<0>(t);
    int const goodId = std::get<1>(t);

    if(bidderId >= 0) {
      Bidder const b = diagram1.at(new_to_old_id[bidderId]);
      double const bx = b.x_;
      double const by = b.y_;
      if(goodId >= 0) {
        Good const g = barycenter.at(goodId);
        double const gx = g.x_;
        double const gy = g.y_;
        barycenter.at(goodId).x_ = (bx + gx) / 2;
        barycenter.at(goodId).y_ = (by + gy) / 2;
        // divide by 4 in order to display the cost of the half matching
        // i.e. the cost of matching to the barycenter
        std::get<2>(t) /= 4;

      } else {
        double gx = (bx + by) / 2;
        double gy = (bx + by) / 2;
        gx = (gx + bx) / 2;
        gy = (gy + by) / 2;
        double const cost = Geometry::pow((gx - bx), wasserstein_)
                            + Geometry::pow((gy - by), wasserstein_);
        MatchingType const t2
          = std::make_tuple(bidderId, barycenter.size(), cost);
        MatchingType const t3 = std::make_tuple(-1, barycenter.size(), cost);
        Good const g = Good(gx, gy, false, barycenter.size());
        // g.SetCriticalCoordinates(b.coords_[0], b.coords_[1], b.coords_[2]);
        barycenter.emplace_back(g);
        matching_to_add.emplace_back(t2);
        matching_to_add2.emplace_back(t3);
      }
    } else {
      if(goodId >= 0) {
        Good const g = barycenter.at(goodId);
        double const gx = (g.x_ + g.y_) / 2;
        double const gy = (g.x_ + g.y_) / 2;
        barycenter.at(goodId).x_ = (gx + g.x_) / 2;
        barycenter.at(goodId).y_ = (gy + g.y_) / 2;
        std::get<2>(t) /= 4;
      }
    }
  }
  for(size_t j = 0; j < matching_to_add.size(); j++) {
    matchings1.emplace_back(matching_to_add[j]);
  }
  for(size_t j = 0; j < matching_to_add2.size(); j++) {
    matchings0.emplace_back(matching_to_add2[j]);
  }

  // correct the costs of matchings in diagram 0
  // costs are initially allzeros because the barycenter is identical
  // to diagram 0
  std::vector<int> new_to_old_id2(current_diagram0.size());
  // 1. Invert the current_bidder_ids_ vector
  for(size_t j = 0; j < ids0.size(); j++) {
    int const new_id = ids0[j];
    if(new_id >= 0) {
      new_to_old_id2[new_id] = j;
    }
  }

  for(size_t i = 0; i < matchings0.size(); i++) {
    MatchingType &t = matchings0[i];
    int const bidderId = std::get<0>(t);
    int const goodId = std::get<1>(t);

    if(bidderId >= 0 and goodId >= 0) {
      Bidder const b = diagram0.at(new_to_old_id2[bidderId]);
      double const bx = b.x_;
      double const by = b.y_;
      Good const g = barycenter.at(goodId);
      double const gx = g.x_;
      double const gy = g.y_;
      double const cost = Geometry::pow((gx - bx), wasserstein_)
                          + Geometry::pow((gy - by), wasserstein_);
      std::get<2>(t) = cost;
    }
  }
}
