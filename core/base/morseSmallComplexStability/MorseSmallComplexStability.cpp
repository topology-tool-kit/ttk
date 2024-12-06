#include <MorseSmallComplexStability.h>
#include <AssignmentAuction.h>

#include <cmath>
ttk::MorseSmallComplexStability::MorseSmallComplexStability() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("MorseSmallComplexStability");
}

void ttk::MorseSmallComplexStability::assignmentSolver(
  std::vector<std::vector<double>> &costMatrix,
  std::vector<ttk::MatchingType> &matching) {
  if(costMatrix.size() > 0) {
    ttk::AssignmentAuction<double> solver;
    solver.setInput(costMatrix);
    solver.run(matching);
    solver.clearMatrix();
  }
}

void ttk::MorseSmallComplexStability::buildCostMatrix(const std::vector<std::array<double, 3>> &coords1,
                                                        const std::vector<std::array<double, 3>> &coords2,
                                                        std::vector<std::vector<double>> &costMatrix){
  int size = coords1.size();
  costMatrix.resize(size);
  for (int i = 0 ; i < size; i++){
    costMatrix[i].resize(size);
  }
  for (int i = 0 ; i < size ; i++){
    for (int j = 0 ; j < size ; j++){
      costMatrix[i][j]=std::sqrt(std::pow(coords1[i][0] - coords2[j][0], 2)
                        + std::pow(coords1[i][1] - coords2[j][1], 2)
                        + std::pow(coords1[i][2] - coords2[j][2], 2));
    }
  }
}

int ttk::MorseSmallComplexStability::buildMatchingsWithOtherBlocks(const std::vector<std::vector<std::array<double, 3>>> &coords, 
                                                                    const int &block_id,
                                                                    std::vector<std::vector<MatchingType>> &matchings){

  int n_blocks = coords.size();

  for (int i = 0 ; i < n_blocks ; i++){
    if(i==block_id)continue;
    std::vector<std::vector<double>> costMatrix;
    int blockIdInMatchingVector = i < block_id ? i : i-1; 
    buildCostMatrix(coords[block_id], coords[i], costMatrix);
    assignmentSolver(costMatrix, matchings[blockIdInMatchingVector]);
  }
  return 1;
}

int ttk::MorseSmallComplexStability::buildOccurenceArraysFull(
                                  const std::vector<GraphMatrixFull> &adjacencyMatricesFull, 
                                  const int &n_separatrices, 
                                  const std::vector<std::vector<std::array<double, 3>>> &coordsSource,
                                  const std::vector<std::vector<std::array<double, 3>>> &coordsDestination,
                                  const int &block_id,
                                  std::vector<int> &edgeOccurences){

  int n_blocks = adjacencyMatricesFull.size();

  for (int i = 0 ; i < n_blocks - 1 ; i++){
    if (coordsSource[i].size() != coordsSource[i+1].size() 
          || coordsDestination[i].size() != coordsDestination[i+1].size()){
      printErr("All blocks must have same number of vertex.");
      return 0;
    }
  }

  int n_source = adjacencyMatricesFull[0].size();
  int n_destination = adjacencyMatricesFull[0][0].size();

  std::vector<std::vector<MatchingType>> matchingsSource(n_blocks - 1);
  std::vector<std::vector<MatchingType>> matchingsDestination(n_blocks - 1);
  
  buildMatchingsWithOtherBlocks(coordsSource, block_id, matchingsSource);
  buildMatchingsWithOtherBlocks(coordsDestination, block_id, matchingsDestination);

  edgeOccurences.resize(n_separatrices, 1);

  for (int k = 0 ; k < n_blocks; k++){
  if(k==block_id)continue;
  int otherBlockdIdMatchingsVector = k < block_id ? k : k-1; 
    for (int i = 0 ; i < n_source ; i++){
      for (int j = 0 ; j < n_destination; j++){
        int thisBlockSourceId = std::get<0>(matchingsSource[otherBlockdIdMatchingsVector][i]);
        int thisBlockDestinationId = std::get<0>(matchingsDestination[otherBlockdIdMatchingsVector][j]);
        if(adjacencyMatricesFull[block_id][thisBlockSourceId][thisBlockDestinationId]!=0){
          
          int otherBlockSourceId = std::get<1>(matchingsSource[otherBlockdIdMatchingsVector][i]);
          int otherBlockDestinationId = std::get<1>(matchingsDestination[otherBlockdIdMatchingsVector][j]);

          if(adjacencyMatricesFull[k][otherBlockSourceId][otherBlockDestinationId]!=0){
            int separatriceId = adjacencyMatricesFull[block_id][thisBlockSourceId][thisBlockDestinationId];
            edgeOccurences[separatriceId]++;
          }
        }
      }
    }
  }
  return 1;
}


int ttk::MorseSmallComplexStability::buildOccurenceArraysMinor(const std::vector<GraphMatrixMinor> &adjacencyMatricesMinor, 
                                                            const int &n_separatrices,
                                                            const std::vector<std::vector<std::array<double, 3>>> &coords,
                                                            const int &block_id,
                                                            std::vector<int> &edgeOccurences){
  int n_blocks = adjacencyMatricesMinor.size();
  int n_points = adjacencyMatricesMinor[0].size();

  for (int i = 0 ; i < n_blocks - 1 ; i++){
    if (coords[i].size() != coords[i+1].size()){
      printErr("All blocks must have same number of vertex.\n");
      return 0;
    }
  }
  std::vector<std::vector<MatchingType>> matchings(n_blocks-1);
  buildMatchingsWithOtherBlocks(coords, block_id, matchings);

  edgeOccurences.resize(n_separatrices);

  for (int i = 0 ; i < n_points ; i ++){
    for (int j = 0 ; j < n_points ; j++){
      int n_edge_ij = adjacencyMatricesMinor[block_id][i][j].size(); 
      for (auto edge : adjacencyMatricesMinor[block_id][i][j]){
        edgeOccurences[edge.first] = n_edge_ij; 
        edgeOccurences[edge.second] = n_edge_ij; 
        }
      }
    }
  
  for (int k = 0 ; k < n_blocks ; k++){

    if(k==block_id)continue;

    int otherBlockdIdMatchingsVector = k < block_id ? k : k-1; 

    for (int i = 0 ; i < n_points ; i++){
      for (int j = 0 ; j < n_points ; j++){

        int tmp1 = std::get<0>(matchings[otherBlockdIdMatchingsVector][i]);
        int tmp2 = std::get<0>(matchings[otherBlockdIdMatchingsVector][j]);
        int thisBlockVertex1 = std::min(tmp1, tmp2);
        int thisBlockVertex2 = std::max(tmp1, tmp2); 

        if(!adjacencyMatricesMinor[block_id][thisBlockVertex1][thisBlockVertex2].empty()){

          tmp1 = std::get<1>(matchings[otherBlockdIdMatchingsVector][i]);
          tmp2 = std::get<1>(matchings[otherBlockdIdMatchingsVector][j]);
          int otherBlockVertex1 = std::min(tmp1, tmp2);
          int otherBlockVertex2 = std::max(tmp1, tmp2); 

          int occurenceInOtherBlock=adjacencyMatricesMinor[k][otherBlockVertex1][otherBlockVertex2].size();

          for (auto edge : adjacencyMatricesMinor[block_id][thisBlockVertex1][thisBlockVertex2]){
            edgeOccurences[edge.first] += occurenceInOtherBlock;
            edgeOccurences[edge.second] += occurenceInOtherBlock;
          }

        }
      }
    } 
  }
  return 1;
} 