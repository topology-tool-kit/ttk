#include <AssignmentAuction.h>
#include <AssignmentMunkres.h>
#include <TrackingFromCriticalPoints.h>
#include <algorithm>

// Compute L_p distance betweem (p,f(p)) and (q,f(q)) where p and q are critical
// points

double ttk::TrackingFromCriticalPoints::criticalPointDistance(
  const std::array<float, 3> &coords_p1,
  const double &sfValue_p1,
  const std::array<float, 3> &coords_p2,
  const double &sfValue_p2,
  const int &p = 2) {
  return std::pow(xWeight_ * std::pow(coords_p1[0] - coords_p2[0], p)
                    + yWeight_ * std::pow(coords_p1[1] - coords_p2[1], p)
                    + zWeight_ * std::pow(coords_p1[2] - coords_p2[2], p)
                    + fWeight_ * std::pow(sfValue_p1 - sfValue_p2, p),
                  1.0 / p);
}

// Sort the critical points by types

void ttk::TrackingFromCriticalPoints::sortCriticalPoint(
  const DiagramType &d,
  const double minimumRelevantPersistence,
  std::vector<std::array<float, 3>> &maxCoords,
  std::vector<std::array<float, 3>> &sad_1Coords,
  std::vector<std::array<float, 3>> &sad_2Coords,
  std::vector<std::array<float, 3>> &minCoords,
  std::vector<double> &maxScalar,
  std::vector<double> &sad_1Scalar,
  std::vector<double> &sad_2Scalar,
  std::vector<double> &minScalar,
  std::vector<SimplexId> &mapMax,
  std::vector<SimplexId> &mapSad_1,
  std::vector<SimplexId> &mapSad_2,
  std::vector<SimplexId> &mapMin) {
  for(unsigned int i = 0; i < d.size(); i++) {
    std::array<float, 3> birthCoords = d[i].birth.coords;
    std::array<float, 3> deathCoords = d[i].death.coords;
    if(std::abs(d[i].persistence()) > minimumRelevantPersistence) {
      switch(d[i].dim) {
        case 0:
          minCoords.push_back(birthCoords);
          minScalar.push_back(d[i].birth.sfValue);
          mapMin.push_back(i);

          sad_1Coords.push_back(deathCoords);
          sad_1Scalar.push_back(d[i].death.sfValue);
          mapSad_1.push_back(i);
          break;
        case 1:
          sad_2Coords.push_back(birthCoords);
          sad_2Scalar.push_back(d[i].birth.sfValue);
          mapSad_2.push_back(i);

          maxCoords.push_back(deathCoords);
          maxScalar.push_back(d[i].death.sfValue);
          mapMax.push_back(i);
          break;
        case 2:
          sad_1Coords.push_back(birthCoords);
          sad_1Scalar.push_back(d[i].birth.sfValue);
          mapSad_1.push_back(i);

          sad_2Coords.push_back(deathCoords);
          sad_2Scalar.push_back(d[i].death.sfValue);
          mapSad_2.push_back(i);
          break;
      }
    }
  }
}

void ttk::TrackingFromCriticalPoints::buildCostMatrix(
  const std::vector<std::array<float, 3>> &coords_1,
  const std::vector<double> &sfValues_1,
  const std::vector<std::array<float, 3>> &coords_2,
  const std::vector<double> &sfValues_2,
  const float &costDeathBirth,
  std::vector<std::vector<double>> &matrix) {
  int size_1 = coords_1.size();
  int size_2 = coords_2.size();
  int matrix_size = (size_1 > 0 && size_2 > 0) ? size_1 + size_2 : 0;
  for(int i = 0; i < size_1; i++) {
    for(int j = 0; j < size_2; j++) {
      matrix[i][j] = criticalPointDistance(
        coords_1[i], sfValues_1[i], coords_2[j], sfValues_2[j]);
    }
  }
  if(!adaptiveDeathBirthCost_) {
    for(int i = size_1; i < matrix_size; i++) {
      for(int j = 0; j < size_2; j++) {
        matrix[i][j] = costDeathBirth;
      }
    }
    for(int i = 0; i < size_1; i++) {
      for(int j = size_2; j < matrix_size; j++) {
        matrix[i][j] = costDeathBirth;
      }
    }
  } else if(adaptiveDeathBirthCost_) {
    for(int j = 0; j < size_2; j++) {
      double c = matrix[0][j];
      for(int i = 1; i < size_1; i++) {
        c = matrix[i][j] < c ? matrix[i][j] : c;
      }
      for(int i = size_1; i < matrix_size; i++) {
        matrix[i][j] = c / (epsilonAdapt_);
      }
    }

    for(int i = 0; i < size_1; i++) {
      double d = matrix[i][0];
      for(int j = 1; j < size_2; j++) {
        d = matrix[i][j] < d ? matrix[i][j] : d;
      }
      for(int j = size_2; j < matrix_size; j++) {
        matrix[i][j] = d / (epsilonAdapt_);
      }
    }
  }
}

void ttk::TrackingFromCriticalPoints::performMatchings(
  const std::vector<DiagramType> &persistenceDiagrams,
  std::vector<std::vector<MatchingType>> &maximaMatchings,
  std::vector<std::vector<MatchingType>> &sad_1_Matchings,
  std::vector<std::vector<MatchingType>> &sad_2_Matchings,
  std::vector<std::vector<MatchingType>> &minimaMatchings,
  std::vector<std::vector<SimplexId>> &maxMap,
  std::vector<std::vector<SimplexId>> &sad_1Map,
  std::vector<std::vector<SimplexId>> &sad_2Map,
  std::vector<std::vector<SimplexId>> &minMap) {

  int fieldNumber = persistenceDiagrams.size();

  std::vector<std::vector<std::array<float, 3>>> maxCoords(fieldNumber);
  std::vector<std::vector<std::array<float, 3>>> sad_1Coords(fieldNumber);
  std::vector<std::vector<std::array<float, 3>>> sad_2Coords(fieldNumber);
  std::vector<std::vector<std::array<float, 3>>> minCoords(fieldNumber);

  std::vector<std::vector<double>> maxScalar(fieldNumber);
  std::vector<std::vector<double>> sad_1Scalar(fieldNumber);
  std::vector<std::vector<double>> sad_2Scalar(fieldNumber);
  std::vector<std::vector<double>> minScalar(fieldNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < fieldNumber; i++) {

    double minimumRelevantPersistence{};
    if(i < fieldNumber - 1) {
      minimumRelevantPersistence
        = ttk::TrackingFromCriticalPoints::computeRelevantPersistence(
          persistenceDiagrams[i], persistenceDiagrams[i + 1]);
    }
    if(i == fieldNumber - 1) {
      minimumRelevantPersistence
        = ttk::TrackingFromCriticalPoints::computeRelevantPersistence(
          persistenceDiagrams[i - 1], persistenceDiagrams[i]);
    }
    sortCriticalPoint(persistenceDiagrams[i], minimumRelevantPersistence,
                      maxCoords[i], sad_1Coords[i], sad_2Coords[i],
                      minCoords[i], maxScalar[i], sad_1Scalar[i],
                      sad_2Scalar[i], minScalar[i], maxMap[i], sad_1Map[i],
                      sad_2Map[i], minMap[i]);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < fieldNumber - 1; i++) {

    float costDeathBirth
      = epsilonConstant_
        * computeBoundingBoxRadius(
          persistenceDiagrams[i], persistenceDiagrams[i + 1]);
    int maxSize = (maxCoords[i].size() > 0 && maxCoords[i + 1].size() > 0)
                    ? maxCoords[i].size() + maxCoords[i + 1].size()
                    : 0;
    int sad_1Size = (sad_1Coords[i].size() > 0 && sad_1Coords[i + 1].size() > 0)
                      ? sad_1Coords[i].size() + sad_1Coords[i + 1].size()
                      : 0;
    int sad_2Size = (sad_2Coords[i].size() > 0 && sad_2Coords[i + 1].size() > 0)
                      ? sad_2Coords[i].size() + sad_2Coords[i + 1].size()
                      : 0;
    int minSize = (minCoords[i].size() > 0 && minCoords[i + 1].size() > 0)
                    ? minCoords[i].size() + minCoords[i + 1].size()
                    : 0;

    std::vector<std::vector<double>> maxMatrix(
      maxSize, std::vector<double>(maxSize, 0));
    std::vector<std::vector<double>> sad_1Matrix(
      sad_1Size, std::vector<double>(sad_1Size, 0));
    std::vector<std::vector<double>> sad_2Matrix(
      sad_2Size, std::vector<double>(sad_2Size, 0));
    std::vector<std::vector<double>> minMatrix(
      minSize, std::vector<double>(minSize, 0));

    buildCostMatrix(maxCoords[i], maxScalar[i], maxCoords[i + 1],
                    maxScalar[i + 1], costDeathBirth, maxMatrix);
    buildCostMatrix(sad_1Coords[i], sad_1Scalar[i], sad_1Coords[i + 1],
                    sad_1Scalar[i + 1], costDeathBirth, sad_1Matrix);
    buildCostMatrix(sad_2Coords[i], sad_2Scalar[i], sad_2Coords[i + 1],
                    sad_2Scalar[i + 1], costDeathBirth, sad_2Matrix);
    buildCostMatrix(minCoords[i], minScalar[i], minCoords[i + 1],
                    minScalar[i + 1], costDeathBirth, minMatrix);

    assignmentSolver(maxMatrix, maximaMatchings[i]);
    assignmentSolver(sad_1Matrix, sad_1_Matchings[i]);
    assignmentSolver(sad_2Matrix, sad_2_Matchings[i]);
    assignmentSolver(minMatrix, minimaMatchings[i]);
  }
}

int ttk::TrackingFromCriticalPoints::computeGlobalId(
  const DiagramType &persistenceDiagram,
  const CriticalType &type,
  const SimplexId &id) {

  switch(persistenceDiagram[id].dim) {
    case 0:
      return (type == CriticalType::Local_minimum
                ? persistenceDiagram[id].birth.id
                : persistenceDiagram[id].death.id);
    case 1:
      return (type == CriticalType::Saddle2 ? persistenceDiagram[id].birth.id
                                            : persistenceDiagram[id].death.id);
    case 2:
      return (type == CriticalType::Saddle1 ? persistenceDiagram[id].birth.id
                                            : persistenceDiagram[id].death.id);
  }
  return -1;
}

void ttk::TrackingFromCriticalPoints::performTrackingForOneType(
  const std::vector<DiagramType> &persistenceDiagrams,
  const std::vector<std::vector<MatchingType>> &matchings,
  const std::vector<std::vector<SimplexId>> &map,
  const CriticalType &currentType,
  std::vector<trackingTuple> &trackings,
  std::vector<std::vector<double>> &trackingCosts,
  std::vector<double> &trackingPersistence) {
  int fieldNumber = matchings.size() + 1;

  SimplexId deathLimitId = map[1].size();
  SimplexId birthLimitId = map[0].size();

  std::vector<int> previousStepMap(deathLimitId, -1);
  std::vector<int> sw;

  for(unsigned int i = 0; i < matchings[0].size(); i++) {
    SimplexId endLocalId = std::get<1>(matchings[0][i]);
    SimplexId startLocalId = std::get<0>(matchings[0][i]);
    if(endLocalId < deathLimitId && startLocalId < birthLimitId) {

      SimplexId startId = computeGlobalId(
        persistenceDiagrams[0], currentType, map[0][startLocalId]);

      SimplexId endId = computeGlobalId(
        persistenceDiagrams[1], currentType, map[1][endLocalId]);

      std::vector<SimplexId> chain = {startId, endId};
      std::tuple<int, int, std::vector<SimplexId>> tt
        = std::make_tuple(0, 1, chain);
      trackings.push_back(tt);
      trackingPersistence.push_back(
        persistenceDiagrams[0][map[0][startLocalId]].persistence()
        + persistenceDiagrams[1][map[1][endLocalId]].persistence());
      std::vector<double> newCostEntry;
      newCostEntry.push_back(std::get<2>(matchings[0][i]));

      trackingCosts.push_back(newCostEntry);

      previousStepMap[std::get<1>(matchings[0][i])] = trackings.size() - 1;
    }
  }
  for(int i = 1; i < fieldNumber - 1; i++) {
    birthLimitId = map[i].size();
    deathLimitId = map[i + 1].size();
    sw.resize(deathLimitId, -1);
    for(unsigned int j = 0; j < matchings[i].size(); j++) {
      SimplexId startLocalId = std::get<0>(matchings[i][j]);
      SimplexId endLocalId = std::get<1>(matchings[i][j]);

      bool wasPreviouslyMatched = false;
      if(startLocalId < birthLimitId)
        wasPreviouslyMatched = previousStepMap[startLocalId] != -1;

      bool validPair = startLocalId < birthLimitId && endLocalId < deathLimitId;

      if(validPair && wasPreviouslyMatched) {
        int trackingId = previousStepMap[startLocalId];
        SimplexId endId = computeGlobalId(
          persistenceDiagrams[i + 1], currentType, map[i + 1][endLocalId]);
        std::get<1>(trackings[trackingId])++;
        std::get<2>(trackings[trackingId]).push_back(endId);
        trackingCosts[trackingId].push_back(std::get<2>(matchings[i][j]));
        trackingPersistence[trackingId]
          += persistenceDiagrams[i + 1][map[i + 1][endLocalId]].persistence();
        sw[endLocalId] = trackingId;
        previousStepMap[startLocalId] = -1;
      }

      else if(validPair && !wasPreviouslyMatched) {

        SimplexId startId = computeGlobalId(
          persistenceDiagrams[i], currentType, map[i][startLocalId]);
        SimplexId endId = computeGlobalId(
          persistenceDiagrams[i + 1], currentType, map[i + 1][endLocalId]);
        std::vector<ttk::SimplexId> chain = {startId, endId};
        std::tuple<int, int, std::vector<SimplexId>> tt
          = std::make_tuple(i, i + 1, chain);
        trackings.push_back(tt);
        std::vector<double> newCostEntry;
        newCostEntry.push_back(std::get<2>(matchings[i][j]));
        trackingCosts.push_back(newCostEntry);
        trackingPersistence.push_back(
          persistenceDiagrams[i][map[i][startLocalId]].persistence()
          + persistenceDiagrams[i + 1][map[i + 1][endLocalId]].persistence());
        sw[endLocalId] = trackings.size() - 1;
      }
    }
    previousStepMap = sw;
    sw.clear();
  }

  for(unsigned int i = 0; i < trackings.size(); i++) {
    trackingPersistence[i]
      /= (std::get<1>(trackings[i]) - std::get<0>(trackings[i]) + 1);
  }
}

void ttk::TrackingFromCriticalPoints::performTrackings(
  const std::vector<DiagramType> &persistenceDiagrams,
  const std::vector<std::vector<MatchingType>> &maximaMatchings,
  const std::vector<std::vector<MatchingType>> &sad_1_Matchings,
  const std::vector<std::vector<MatchingType>> &sad_2_Matchings,
  const std::vector<std::vector<MatchingType>> &minimaMatchings,
  const std::vector<std::vector<SimplexId>> &maxMap,
  const std::vector<std::vector<SimplexId>> &sad_1Map,
  const std::vector<std::vector<SimplexId>> &sad_2Map,
  const std::vector<std::vector<SimplexId>> &minMap,
  std::vector<trackingTuple> &allTrackings,
  std::vector<std::vector<double>> &allTrackingsCosts,
  std::vector<double> &allTrackingsMeanPersistences,
  unsigned int (&typesArrayLimits)[3]) {

  std::vector<ttk::trackingTuple> trackingsMax;
  std::vector<ttk::trackingTuple> trackingsSad_1;
  std::vector<ttk::trackingTuple> trackingsSad_2;
  std::vector<ttk::trackingTuple> trackingsMin;

  std::vector<std::vector<double>> maxTrackingCost;
  std::vector<std::vector<double>> sad_1_TrackingCost;
  std::vector<std::vector<double>> sad_2_TrackingCost;
  std::vector<std::vector<double>> minTrackingCost;

  std::vector<double> trackingsPersistenceMax;
  std::vector<double> trackingsPersistenceSad_1;
  std::vector<double> trackingsPersistenceSad_2;
  std::vector<double> trackingsPersistenceMin;

  performTrackingForOneType(persistenceDiagrams, maximaMatchings, maxMap,
                            CriticalType::Local_maximum, trackingsMax,
                            maxTrackingCost, trackingsPersistenceMax);
  allTrackings.insert(
    allTrackings.end(), trackingsMax.begin(), trackingsMax.end());
  allTrackingsMeanPersistences.insert(allTrackingsMeanPersistences.end(),
                                      trackingsPersistenceMax.begin(),
                                      trackingsPersistenceMax.end());
  allTrackingsCosts.insert(
    allTrackingsCosts.end(), maxTrackingCost.begin(), maxTrackingCost.end());
  typesArrayLimits[0] = allTrackings.size();

  performTrackingForOneType(persistenceDiagrams, sad_1_Matchings, sad_1Map,
                            CriticalType::Saddle1, trackingsSad_1,
                            sad_1_TrackingCost, trackingsPersistenceSad_1);
  allTrackings.insert(
    allTrackings.end(), trackingsSad_1.begin(), trackingsSad_1.end());
  allTrackingsMeanPersistences.insert(allTrackingsMeanPersistences.end(),
                                      trackingsPersistenceSad_1.begin(),
                                      trackingsPersistenceSad_1.end());
  allTrackingsCosts.insert(allTrackingsCosts.end(), sad_1_TrackingCost.begin(),
                           sad_1_TrackingCost.end());
  typesArrayLimits[1] = allTrackings.size();

  performTrackingForOneType(persistenceDiagrams, sad_2_Matchings, sad_2Map,
                            CriticalType::Saddle2, trackingsSad_2,
                            sad_2_TrackingCost, trackingsPersistenceSad_2);
  allTrackings.insert(
    allTrackings.end(), trackingsSad_2.begin(), trackingsSad_2.end());
  allTrackingsMeanPersistences.insert(allTrackingsMeanPersistences.end(),
                                      trackingsPersistenceSad_2.begin(),
                                      trackingsPersistenceSad_2.end());
  allTrackingsCosts.insert(allTrackingsCosts.end(), sad_2_TrackingCost.begin(),
                           sad_2_TrackingCost.end());
  typesArrayLimits[2] = allTrackings.size();

  performTrackingForOneType(persistenceDiagrams, minimaMatchings, minMap,
                            CriticalType::Local_minimum, trackingsMin,
                            minTrackingCost, trackingsPersistenceMin);
  allTrackings.insert(
    allTrackings.end(), trackingsMin.begin(), trackingsMin.end());
  allTrackingsMeanPersistences.insert(allTrackingsMeanPersistences.end(),
                                      trackingsPersistenceMin.begin(),
                                      trackingsPersistenceMin.end());
  allTrackingsCosts.insert(
    allTrackingsCosts.end(), minTrackingCost.begin(), minTrackingCost.end());
}

void ttk::TrackingFromCriticalPoints::assignmentSolver(
  std::vector<std::vector<double>> &costMatrix,
  std::vector<ttk::MatchingType> &matching) {
  if(costMatrix.size() > 0) {
    if(assignmentMethod_ == 0) {
      ttk::AssignmentAuction<double> solver;
      solver.setInput(costMatrix);
      solver.run(matching);
      solver.clearMatrix();
    } else if(assignmentMethod_ == 1) {
      ttk::AssignmentMunkres<double> solver;
      solver.setInput(costMatrix);
      solver.run(matching);
      solver.clearMatrix();
    }
  }
}
