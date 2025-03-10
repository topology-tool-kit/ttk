#include <AssignmentAuction.h>
#include <SeparatrixStability.h>

#include <cmath>
ttk::SeparatrixStability::SeparatrixStability() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("SeparatrixStability");
}

void ttk::SeparatrixStability::assignmentSolver(
  std::vector<std::vector<double>> &costMatrix,
  std::vector<ttk::MatchingType> &matching) {
  if(costMatrix.size() > 0) {
    ttk::AssignmentAuction<double> solver;
    solver.setInput(costMatrix);
    solver.run(matching);
    solver.clearMatrix();
  }
}

void ttk::SeparatrixStability::buildCostMatrix(
  const std::vector<std::array<double, 3>> &coords1,
  const std::vector<std::array<double, 3>> &coords2,
  std::vector<std::vector<double>> &costMatrix) {
  int size = coords1.size();
  costMatrix.resize(size);
  for(int i = 0; i < size; i++) {
    costMatrix[i].resize(size);
  }
  for(int i = 0; i < size; i++) {
    for(int j = 0; j < size; j++) {
      costMatrix[i][j]
        = std::sqrt(std::pow(coords1[i][0] - coords2[j][0], 2)
                    + std::pow(coords1[i][1] - coords2[j][1], 2)
                    + std::pow(coords1[i][2] - coords2[j][2], 2));
    }
  }
}

int ttk::SeparatrixStability::buildMatchingsWithOtherBlocks(
  const std::vector<std::vector<std::array<double, 3>>> &coords,
  const int &block_id,
  std::vector<std::vector<MatchingType>> &matchings) {

  int n_blocks = coords.size();

  for(int i = 0; i < n_blocks; i++) {
    if(i == block_id)
      continue;
    std::vector<std::vector<double>> costMatrix;
    int blockIdInMatchingVector = i < block_id ? i : i - 1;
    buildCostMatrix(coords[block_id], coords[i], costMatrix);
    assignmentSolver(costMatrix, matchings[blockIdInMatchingVector]);
  }
  return 1;
}

int ttk::SeparatrixStability::buildOccurenceArraysFull(
  const std::vector<GraphMatrixFull> &ajacencyMatrices,
  const int &n_separatrices,
  const std::vector<std::vector<std::array<double, 3>>> &coordsSource,
  const std::vector<std::vector<std::array<double, 3>>> &coordsDestination,
  const int &block_id,
  std::vector<int> &edgeOccurences,
  std::vector<bool> &isIsomorphicWith,
  std::vector<std::vector<int>> &matchingArraySource,
  std::vector<std::vector<int>> &matchingArrayDestination,
  std::vector<std::vector<int>> &matchingArraySeparatrix){

  int n_blocks = ajacencyMatrices.size();
  for(int i = 0; i < n_blocks - 1; i++) {
    if(coordsSource[i].size() != coordsSource[i + 1].size()
       || coordsDestination[i].size() != coordsDestination[i + 1].size()) {
      printErr("All blocks must have same number of vertex.");
      return 0;
    }
  }

  int n_source = ajacencyMatrices[0].size();
  int n_destination = ajacencyMatrices[0][0].size();
  isIsomorphicWith.resize(n_blocks, true);
  matchingArraySource.resize(n_blocks, std::vector<int>(n_source));
  matchingArrayDestination.resize(n_blocks, std::vector<int>(n_destination));



  std::vector<std::vector<MatchingType>> matchingsSource(n_blocks - 1);
  std::vector<std::vector<MatchingType>> matchingsDestination(n_blocks - 1);

  buildMatchingsWithOtherBlocks(coordsSource, block_id, matchingsSource);
  buildMatchingsWithOtherBlocks(
    coordsDestination, block_id, matchingsDestination);

  edgeOccurences.resize(n_separatrices, 1);
  matchingArraySeparatrix.resize(n_blocks, std::vector<int>(n_separatrices, -1));

  for(int k = 0; k < n_blocks; k++) {
    if(k == block_id){
      for (int i = 0 ; i < n_destination; i++){
        matchingArrayDestination[k][i]=i;
      }
      for (int i = 0 ; i < n_source; i++){
        matchingArraySource[k][i]=i;
      }
      for (int i = 0 ; i < n_separatrices; i++){
        matchingArraySeparatrix[k][i]=i;
      }
      continue;
    }
    int otherBlockdIdMatchingsVector = k < block_id ? k : k - 1;
    for(int i = 0; i < n_source; i++) {
      for(int j = 0; j < n_destination; j++) {
        int thisBlockSourceId
          = std::get<0>(matchingsSource[otherBlockdIdMatchingsVector][i]);
        int thisBlockDestinationId
          = std::get<0>(matchingsDestination[otherBlockdIdMatchingsVector][j]);
        int otherBlockSourceId
          = std::get<1>(matchingsSource[otherBlockdIdMatchingsVector][i]);
        int otherBlockDestinationId = std::get<1>(
          matchingsDestination[otherBlockdIdMatchingsVector][j]);
        matchingArraySource[k][otherBlockSourceId]=thisBlockSourceId;
        matchingArrayDestination[k][otherBlockDestinationId]=thisBlockDestinationId;
        if(ajacencyMatrices[block_id][thisBlockSourceId]
                                [thisBlockDestinationId]
          != -1 &&
          ajacencyMatrices[k][otherBlockSourceId]
                                 [otherBlockDestinationId]
            != -1) {          
          int separatriceId
              = ajacencyMatrices[block_id][thisBlockSourceId]
                                     [thisBlockDestinationId];
          int separatriceIdInOtherBlock = ajacencyMatrices[k][otherBlockSourceId]
          [otherBlockDestinationId];
          edgeOccurences[separatriceId]++;
          matchingArraySeparatrix[k][separatriceIdInOtherBlock]= separatriceId;
        }else if((ajacencyMatrices[block_id][thisBlockSourceId]
                                [thisBlockDestinationId]
          != -1) ^
          (ajacencyMatrices[k][otherBlockSourceId]
                                 [otherBlockDestinationId]
            != -1)){
              isIsomorphicWith[k]=false;
            }
      }
    }
  }
  //re-index the matching id for the source point so that they are different from the destination point id
  for (int i = 0 ; i < n_blocks; i++){
    int sourceSize = matchingArraySource[i].size();
    int destinationSize = matchingArrayDestination[i].size();
    for (int j = 0 ; j < sourceSize; j++){
      matchingArraySource[i][j]+=destinationSize;
    }
  }
  return 1;
}

void ttk::SeparatrixStability::computeGraphMinor(
  const GraphMatrixFull &adjacencyMatrixFull,
  GraphMatrixMinor &adjacencyMatrix) {
  int n_row = adjacencyMatrixFull.size();
  int n_col = adjacencyMatrixFull[0].size();
  adjacencyMatrix.resize(n_col);
  for(int i = 0; i < n_col; i++) {
    adjacencyMatrix[i].resize(n_col);
  }

  for(int i = 0; i < n_row; i++) {
    for(int j = 0; j < n_col; j++) {
      if(adjacencyMatrixFull[i][j] != -1) {
        for(int k = 0; k < j; k++) {
          if(adjacencyMatrixFull[i][k] != -1) {
            std::pair<int, int> newEdge = std::make_pair(
              adjacencyMatrixFull[i][j], adjacencyMatrixFull[i][k]);
            adjacencyMatrix[k][j].push_back(newEdge);
          }
        }
      }
    }
  }
}

int ttk::SeparatrixStability::buildOccurenceArraysMinor(
  const std::vector<GraphMatrixFull> &adjacencyMatricesFull,
  const int &n_separatrices,
  const std::vector<std::vector<std::array<double, 3>>> &coords,
  const int &block_id,
  std::vector<int> &edgeOccurences,
  std::vector<bool> &isIsomorphicWith,
  std::vector<std::vector<int>> &matchingArray,
  std::vector<std::vector<int>> &matchingArraySeparatrix) {

  int n_blocks = adjacencyMatricesFull.size();
  isIsomorphicWith.resize(n_blocks, true);

  std::vector<GraphMatrixMinor> adjacencyMatricesMinor(n_blocks);

  for (int i = 0 ; i < n_blocks; i++){
    computeGraphMinor(adjacencyMatricesFull[i], adjacencyMatricesMinor[i]);
  }

  int n_points = adjacencyMatricesMinor[0].size();
  matchingArray.resize(n_blocks, std::vector<int>(n_points));

  for(int i = 0; i < n_blocks - 1; i++) {
    if(coords[i].size() != coords[i + 1].size()) {
      printErr("All blocks must have same number of vertex.\n");
      return 0;
    }
  }
  std::vector<std::vector<MatchingType>> matchings(n_blocks - 1);
  buildMatchingsWithOtherBlocks(coords, block_id, matchings);

  edgeOccurences.resize(n_separatrices);

  for(int i = 0; i < n_points; i++) {
    for(int j = 0; j < n_points; j++) {
      if(!adjacencyMatricesMinor[block_id][i][j].empty()) {;
        for(auto edge : adjacencyMatricesMinor[block_id][i][j]) {
          edgeOccurences[edge.first] = 1;
          edgeOccurences[edge.second] = 1;
        }
      }
    }
  }

  for(int k = 0; k < n_blocks; k++) {

    if(k == block_id){
      for (int i = 0 ; i < n_points ; i++){
        matchingArray[k][i]=i;
      }
      for (int i = 0 ; i < n_separatrices; i++){
        matchingArraySeparatrix[k][i]=i;
      }
      continue;
    }

    int otherBlockdIdMatchingsVector = k < block_id ? k : k - 1;

    for (int i = 0; i < n_points; i++){
      int thisBlockPointId = std::get<0>(matchings[otherBlockdIdMatchingsVector][i]);
      int otherBlockPointId = std::get<1>(matchings[otherBlockdIdMatchingsVector][i]);
      matchingArray[k][otherBlockPointId]=thisBlockPointId;
    }

    for(int i = 0; i < n_points; i++) {
      for(int j = i + 1; j < n_points; j++) {
        int tmp2 = std::get<0>(matchings[otherBlockdIdMatchingsVector][j]);
        int tmp1 = std::get<0>(matchings[otherBlockdIdMatchingsVector][i]);
        int thisBlockVertex1 = std::min(tmp1, tmp2);
        int thisBlockVertex2 = std::max(tmp1, tmp2);
        tmp1 = std::get<1>(matchings[otherBlockdIdMatchingsVector][i]);
        tmp2 = std::get<1>(matchings[otherBlockdIdMatchingsVector][j]);
        int otherBlockVertex1 = std::min(tmp1, tmp2);
        int otherBlockVertex2 = std::max(tmp1, tmp2);
        bool existsInOtherBlock
          = !adjacencyMatricesMinor[k][otherBlockVertex1][otherBlockVertex2]
              .empty();

        if(!adjacencyMatricesMinor[block_id][thisBlockVertex1][thisBlockVertex2]
              .empty() && existsInOtherBlock > 0) {
          for(auto edge : adjacencyMatricesMinor[block_id][thisBlockVertex1]
                                                [thisBlockVertex2]) {
            std::pair<int, int> edgeInOtherBlock = adjacencyMatricesMinor[k][otherBlockVertex1][otherBlockVertex2][0];
            edgeOccurences[edge.first] += existsInOtherBlock;
            edgeOccurences[edge.second] += existsInOtherBlock;
            matchingArraySeparatrix[k][edgeInOtherBlock.first]=edge.first;
            matchingArraySeparatrix[k][edgeInOtherBlock.second]=edge.second;
          }
        }
        else if((!adjacencyMatricesMinor[block_id][thisBlockVertex1][thisBlockVertex2]
              .empty()) ^ existsInOtherBlock){
                isIsomorphicWith[k]=false;
              }
      }
    }
  }
  return 1;
}

int ttk::SeparatrixStability::buildOccurenceArrays(
      const std::vector<GraphMatrixFull> &adjacencyMatrices,
      const std::vector<int> &separatrixCountForEachBlock,
      const std::vector<std::vector<std::array<double, 3>>> &coordsSource,
      const std::vector<std::vector<std::array<double, 3>>> &coordsDestination,
      const bool &mergeEdgesOnSaddles,
      std::vector<std::vector<int>> &edgesOccurencesForEachBlock,
      std::vector<std::vector<bool>> &isomorphismForEachBlock,
      std::vector<std::vector<std::vector<int>>> &matchingArrayForEachBlockSource,
      std::vector<std::vector<std::vector<int>>> &matchingArrayForEachBlockDestination,
      std::vector<std::vector<std::vector<int>>> &matchingArraySeparatrixForEachBlock){

  int n_blocks = adjacencyMatrices.size();
  int status;
  for (int i = 0 ; i < n_blocks; i++){
    matchingArraySeparatrixForEachBlock[i].resize(n_blocks);
    for (int j = 0 ; j < n_blocks; j++){
      matchingArraySeparatrixForEachBlock[i][j].resize(separatrixCountForEachBlock[j], -1);
    }
  }

  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(threadNumber_)
  #endif // TTK_ENABLE_OPENMP

  for(int i = 0; i < n_blocks; i++) {
    if(!mergeEdgesOnSaddles) {
      status = this->buildOccurenceArraysFull(
        adjacencyMatrices, separatrixCountForEachBlock[i], coordsSource,
        coordsDestination, i, edgesOccurencesForEachBlock[i], isomorphismForEachBlock[i],
        matchingArrayForEachBlockSource[i],
        matchingArrayForEachBlockDestination[i],
        matchingArraySeparatrixForEachBlock[i]);
    } else {
      status = this->buildOccurenceArraysMinor(
        adjacencyMatrices, separatrixCountForEachBlock[i],
        coordsDestination, i, edgesOccurencesForEachBlock[i], isomorphismForEachBlock[i],
        matchingArrayForEachBlockDestination[i],
        matchingArraySeparatrixForEachBlock[i]);
    }
  }

  return status;
}