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

void ttk::MorseSmallComplexStability::makePartition(const std::vector<std::vector<MatchingType>> &matchings, 
                                                    const int &n_blocks,
                                                    const int &n_points,
                                                    std::vector<std::vector<int>> &classIdToVertexIds){
  classIdToVertexIds.resize(n_points);
  std::vector<int> previousPositionOfId(n_points);
  std::vector<int> switchArray(n_points);
  for (int i = 0 ; i < n_points; i++){
    classIdToVertexIds[i].push_back(std::get<0>(matchings[0][i]));
    classIdToVertexIds[i].push_back(std::get<1>(matchings[0][i]));
    previousPositionOfId[std::get<1>(matchings[0][i])]=i;
  }
  for (int i = 1 ; i < n_blocks - 1 ; i++){
    for (int j = 0 ; j < n_points ; j++){
      int backId = std::get<0>(matchings[i][j]);
      int frontId = std::get<1>(matchings[i][j]);
      classIdToVertexIds[previousPositionOfId[backId]].push_back(frontId);
      switchArray[frontId]=previousPositionOfId[backId];
    }
    previousPositionOfId=switchArray;
    for (int k = 0 ; k < n_points ; k++){
      std::cout<<previousPositionOfId[k]<<"   ";
    }
    std::cout<<std::endl;
  }
}

int ttk::MorseSmallComplexStability::buildVertexEquivalenceClasses(
                                  const std::vector<std::vector<std::array<double, 3>>> &coords, 
                                  std::vector<std::vector<int>> &classIdToVertexIds){

  bool allSameSize = true;
  int n_blocks = coords.size();

  for (int i = 0 ; i < n_blocks - 1 ; i++){
    if (coords[i].size() != coords[i+1].size())allSameSize=false;
    std::cout<<"coords size for block "<<i<<" = "<<coords[i].size()<<std::endl;
    std::cout<<"coords size for block "<<i+1<<" = "<<coords[i+1].size()<<std::endl;
  }
  if(!allSameSize){
    return 0;
  }
  int n_points = coords[0].size();

  std::vector<std::vector<MatchingType>> matchings(n_blocks - 1);

  for (int i = 0 ; i < n_blocks - 1; i++){
    std::vector<std::vector<double>> costMatrix;
    buildCostMatrix(coords[i], coords[i+1], costMatrix);
    assignmentSolver(costMatrix, matchings[i]);
  }
  for (int i = 0 ; i < n_blocks - 1 ; i++){

    std::cout<<"matchings from block "<<i<<"to block "<<i+1<<std::endl;
    for (int j = 0; j < matchings[i].size(); j++){
      std::cout<<std::get<0>(matchings[i][j])<<", "<<std::get<1>(matchings[i][j])<<std::endl;
    }
  }
  makePartition(matchings, n_blocks, n_points, classIdToVertexIds);
  return 1;
}

int ttk::MorseSmallComplexStability::buildOccurenceArrays(const std::vector<GraphMatrix> &adjacencyMatrices, 
                                                            const std::vector<int> &separatrixCountForEachBlock,
                                                            const std::vector<std::vector<std::array<double, 3>>> &coords,
                                                            std::vector<std::vector<int>> &edgeOccurenceForEachBlock){
  int n_blocks = adjacencyMatrices.size();
  int n_points = adjacencyMatrices[0].size();
  
  std::vector<std::vector<int>> classIdToVertexIds(n_points);
  if(!buildVertexEquivalenceClasses(coords, classIdToVertexIds)){
    printErr("All blocks must have same number of vertex.");
    return 0;
  }

  for (int i  = 0; i < n_points ; i++){
    std::cout<<"equivalence class "<<i<<std::endl;
    for (int j = 0 ; j < classIdToVertexIds[i].size(); j++){
      std::cout<<classIdToVertexIds[i][j]<<", ";
    }
    std::cout<<std::endl;
  }
  std::cout<<"balise 1"<<std::endl;

  std::vector<std::vector<int>> occurenceMatrix(n_points);
  for (int i = 0 ; i < n_points ; i++){
    for (int j = 0 ; j < n_points ; j++){
      int occurence=0;
      for (int k = 0 ; k < n_blocks ; k++){
        int vertexId1 = classIdToVertexIds[i][k];
        int vertexId2 = classIdToVertexIds[j][k];
        if(adjacencyMatrices[k][vertexId1][vertexId2].has_value())occurence++;
      }
      occurenceMatrix[i].push_back(occurence);
    }
  }

  for (int i = 0 ; i < n_points ; i++){
    for (int j = 0 ; j < n_points ; j++){
      std::cout<<occurenceMatrix[i][j]<<"  ";
    }
    std::cout<<std::endl;
  }

  for (int i = 0 ; i < n_blocks ; i++){
    edgeOccurenceForEachBlock[i].resize(separatrixCountForEachBlock[i]);
    for (int j = 0 ; j < n_points ; j++){
      for (int k = 0; k < n_points ; k++){
        int vertexId1 = classIdToVertexIds[j][i];
        int vertexId2 = classIdToVertexIds[k][i];
        if (adjacencyMatrices[i][vertexId1][vertexId2].has_value()){
          int separatrixId1 = adjacencyMatrices[i][vertexId1][vertexId2].value().first;
          int separatrixId2 = adjacencyMatrices[i][vertexId1][vertexId2].value().second;
          edgeOccurenceForEachBlock[i][separatrixId1]=occurenceMatrix[j][k];
          edgeOccurenceForEachBlock[i][separatrixId2]=occurenceMatrix[j][k];
        }
      }
    }
  }

  for (int i = 0 ; i < n_blocks ; i++){
    for (int j = 0 ; j < edgeOccurenceForEachBlock[i].size() ; j++){
      std::cout<<edgeOccurenceForEachBlock[i][j]<<"  ";
    }
    std::cout<<std::endl;
  }
  return 1;
} 