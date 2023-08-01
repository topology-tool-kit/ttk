#include <PersistenceDiagramBarycenter.h>

void ttk::PersistenceDiagramBarycenter::execute(
  std::vector<DiagramType> &intermediateDiagrams,
  DiagramType &barycenter,
  std::vector<std::vector<std::vector<MatchingType>>> &all_matchings) {

  Timer tm;

  printMsg("Computing Barycenter of " + std::to_string(numberOfInputs_)
           + " diagrams.");

  std::vector<DiagramType> data_min(numberOfInputs_);
  std::vector<DiagramType> data_sad(numberOfInputs_);
  std::vector<DiagramType> data_max(numberOfInputs_);

  std::vector<std::vector<int>> data_min_idx(numberOfInputs_);
  std::vector<std::vector<int>> data_sad_idx(numberOfInputs_);
  std::vector<std::vector<int>> data_max_idx(numberOfInputs_);

  bool do_min = false;
  bool do_sad = false;
  bool do_max = false;

  // Create diagrams for min, saddle and max persistence pairs
  for(int i = 0; i < numberOfInputs_; i++) {
    DiagramType &CTDiagram = intermediateDiagrams[i];

    for(size_t j = 0; j < CTDiagram.size(); ++j) {
      PersistencePair const &t = CTDiagram[j];

      ttk::CriticalType const nt1 = t.birth.type;
      ttk::CriticalType const nt2 = t.death.type;

      double const dt = t.persistence();
      // if (abs<double>(dt) < zeroThresh) continue;
      if(dt > 0) {
        if(nt1 == ttk::CriticalType::Local_minimum
           && nt2 == ttk::CriticalType::Local_maximum) {
          data_max[i].push_back(t);
          data_max_idx[i].push_back(j);
          do_max = true;
        } else {
          if(nt1 == ttk::CriticalType::Local_maximum
             || nt2 == ttk::CriticalType::Local_maximum) {
            data_max[i].push_back(t);
            data_max_idx[i].push_back(j);
            do_max = true;
          }
          if(nt1 == ttk::CriticalType::Local_minimum
             || nt2 == ttk::CriticalType::Local_minimum) {
            data_min[i].push_back(t);
            data_min_idx[i].push_back(j);
            do_min = true;
          }
          if((nt1 == ttk::CriticalType::Saddle1
              && nt2 == ttk::CriticalType::Saddle2)
             || (nt1 == ttk::CriticalType::Saddle2
                 && nt2 == ttk::CriticalType::Saddle1)) {
            data_sad[i].push_back(t);
            data_sad_idx[i].push_back(j);
            do_sad = true;
          }
        }
      }
    }
  }

  DiagramType barycenter_min;
  DiagramType barycenter_sad;
  DiagramType barycenter_max;

  std::vector<std::vector<MatchingType>> matching_min, matching_sad,
    matching_max;

  double total_cost = 0, min_cost = 0, sad_cost = 0, max_cost = 0;
  /*omp_set_num_threads(1);
  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel sections
  #endif
  {
    #ifdef TTK_ENABLE_OPENMP
    #pragma omp section
    #endif
    {*/
  if(do_min) {
    printMsg("Computing Minima barycenter...");
    PDBarycenter bary_min{};
    bary_min.setThreadNumber(threadNumber_);
    bary_min.setWasserstein(wasserstein_);
    bary_min.setNumberOfInputs(numberOfInputs_);
    bary_min.setDiagramType(0);
    bary_min.setUseProgressive(use_progressive_);
    bary_min.setGeometricalFactor(alpha_);
    bary_min.setDebugLevel(debugLevel_);
    bary_min.setDeterministic(deterministic_);
    bary_min.setLambda(lambda_);
    bary_min.setMethod(method_);
    bary_min.setEarlyStoppage(early_stoppage_);
    bary_min.setEpsilonDecreases(epsilon_decreases_);
    bary_min.setReinitPrices(reinit_prices_);
    bary_min.setDiagrams(&data_min);
    matching_min = bary_min.execute(barycenter_min);
    min_cost = bary_min.getCost();
    total_cost += min_cost;
  }
  /*}

  #ifdef TTK_ENABLE_OPENMP
  #pragma omp section
  #endif
  {*/
  if(do_sad) {
    printMsg("Computing Saddles barycenter...");
    PDBarycenter bary_sad{};
    bary_sad.setThreadNumber(threadNumber_);
    bary_sad.setWasserstein(wasserstein_);
    bary_sad.setNumberOfInputs(numberOfInputs_);
    bary_sad.setDiagramType(1);
    bary_sad.setUseProgressive(use_progressive_);
    bary_sad.setGeometricalFactor(alpha_);
    bary_sad.setLambda(lambda_);
    bary_sad.setDebugLevel(debugLevel_);
    bary_sad.setMethod(method_);
    bary_sad.setEarlyStoppage(early_stoppage_);
    bary_sad.setEpsilonDecreases(epsilon_decreases_);
    bary_sad.setDeterministic(deterministic_);
    bary_sad.setReinitPrices(reinit_prices_);
    bary_sad.setDiagrams(&data_sad);
    matching_sad = bary_sad.execute(barycenter_sad);
    sad_cost = bary_sad.getCost();
    total_cost += sad_cost;
  }
  /*}

  #ifdef TTK_ENABLE_OPENMP
  #pragma omp section
  #endif
  {*/
  if(do_max) {
    printMsg("Computing Maxima barycenter...");
    PDBarycenter bary_max{};
    bary_max.setThreadNumber(threadNumber_);
    bary_max.setWasserstein(wasserstein_);
    bary_max.setNumberOfInputs(numberOfInputs_);
    bary_max.setDiagramType(2);
    bary_max.setUseProgressive(use_progressive_);
    bary_max.setGeometricalFactor(alpha_);
    bary_max.setLambda(lambda_);
    bary_max.setMethod(method_);
    bary_max.setDebugLevel(debugLevel_);
    bary_max.setEarlyStoppage(early_stoppage_);
    bary_max.setDeterministic(deterministic_);
    bary_max.setEpsilonDecreases(epsilon_decreases_);
    bary_max.setReinitPrices(reinit_prices_);
    bary_max.setDiagrams(&data_max);
    matching_max = bary_max.execute(barycenter_max);
    max_cost = bary_max.getCost();
    total_cost += max_cost;
  }
  //}
  //}

  // Reconstruct matchings
  all_matchings.resize(1);
  all_matchings[0].resize(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; i++) {

    if(do_min) {
      for(size_t j = 0; j < matching_min[i].size(); j++) {
        MatchingType t = matching_min[i][j];
        int const bidder_id = std::get<0>(t);
        std::get<0>(t) = data_min_idx[i][bidder_id];
        if(std::get<1>(t) < 0) {
          std::get<1>(t) = -1;
        }
        all_matchings[0][i].push_back(t);
      }
    }

    if(do_sad) {
      for(size_t j = 0; j < matching_sad[i].size(); j++) {
        MatchingType t = matching_sad[i][j];
        int const bidder_id = std::get<0>(t);
        std::get<0>(t) = data_sad_idx[i][bidder_id];
        if(std::get<1>(t) >= 0) {
          std::get<1>(t) = std::get<1>(t) + barycenter_min.size();
        } else {
          std::get<1>(t) = -1;
        }
        all_matchings[0][i].push_back(t);
      }
    }

    if(do_max) {
      for(size_t j = 0; j < matching_max[i].size(); j++) {
        MatchingType t = matching_max[i][j];
        int const bidder_id = std::get<0>(t);
        std::get<0>(t) = data_max_idx[i][bidder_id];
        if(std::get<1>(t) >= 0) {
          std::get<1>(t)
            = std::get<1>(t) + barycenter_min.size() + barycenter_sad.size();
        } else {
          std::get<1>(t) = -1;
        }
        all_matchings[0][i].push_back(t);
      }
    }
  }
  // Reconstruct barcenter
  for(size_t j = 0; j < barycenter_min.size(); j++) {
    const auto &dt = barycenter_min[j];
    barycenter.push_back(dt);
  }
  for(size_t j = 0; j < barycenter_sad.size(); j++) {
    const auto &dt = barycenter_sad[j];
    barycenter.push_back(dt);
  }
  for(size_t j = 0; j < barycenter_max.size(); j++) {
    const auto &dt = barycenter_max[j];
    barycenter.push_back(dt);
  }

  // Recreate 3D critical coordinates of barycentric points
  std::vector<int> number_of_matchings_for_point(barycenter.size());
  std::vector<float> cords_x1(barycenter.size());
  std::vector<float> cords_y1(barycenter.size());
  std::vector<float> cords_z1(barycenter.size());
  std::vector<float> cords_x2(barycenter.size());
  std::vector<float> cords_y2(barycenter.size());
  std::vector<float> cords_z2(barycenter.size());
  for(unsigned i = 0; i < barycenter.size(); i++) {
    number_of_matchings_for_point[i] = 0;
    cords_x1[i] = 0;
    cords_y1[i] = 0;
    cords_z1[i] = 0;
    cords_x2[i] = 0;
    cords_y2[i] = 0;
    cords_z2[i] = 0;
  }

  for(unsigned i = 0; i < all_matchings[0].size(); i++) {
    DiagramType &CTDiagram = intermediateDiagrams[i];
    for(unsigned j = 0; j < all_matchings[0][i].size(); j++) {
      MatchingType t = all_matchings[0][i][j];
      int const bidder_id = std::get<0>(t);
      int const bary_id = std::get<1>(t);

      const auto &bidder = CTDiagram[bidder_id];
      number_of_matchings_for_point[bary_id] += 1;
      cords_x1[bary_id] += bidder.birth.coords[0];
      cords_y1[bary_id] += bidder.birth.coords[1];
      cords_z1[bary_id] += bidder.birth.coords[2];
      cords_x2[bary_id] += bidder.death.coords[0];
      cords_y2[bary_id] += bidder.death.coords[1];
      cords_z2[bary_id] += bidder.death.coords[2];
    }
  }

  for(unsigned i = 0; i < barycenter.size(); i++) {
    if(number_of_matchings_for_point[i] > 0) {
      barycenter[i].birth.coords[0]
        = cords_x1[i] / number_of_matchings_for_point[i];
      barycenter[i].birth.coords[1]
        = cords_y1[i] / number_of_matchings_for_point[i];
      barycenter[i].birth.coords[2]
        = cords_z1[i] / number_of_matchings_for_point[i];
      barycenter[i].death.coords[0]
        = cords_x2[i] / number_of_matchings_for_point[i];
      barycenter[i].death.coords[1]
        = cords_y2[i] / number_of_matchings_for_point[i];
      barycenter[i].death.coords[2]
        = cords_z2[i] / number_of_matchings_for_point[i];
    }
  }

  printMsg("Min-saddle cost    : " + std::to_string(min_cost));
  printMsg("Saddle-saddle cost : " + std::to_string(sad_cost));
  printMsg("Saddle-max cost    : " + std::to_string(max_cost));
  printMsg("Total cost         : " + std::to_string(total_cost));
  printMsg("Complete", 1, tm.getElapsedTime(), threadNumber_);
}
