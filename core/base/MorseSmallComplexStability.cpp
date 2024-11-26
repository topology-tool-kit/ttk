#include <MorseSmallComplexStability.h>

ttk::MorseSmallComplexStability::MorseSmallComplexStability() {
  this->setDebugMsgPrefix("MorseSmallComplexStability");
}

int ttk::MorseSmallComplexStability::computeGraphData(std::vector<std::vector<std::array<double, 3>>> &blockCoords,
                                           std::vector<std::vector<double>> &blockSfValues,
                                           std::vector<std::vector<std::pair<int, int>>> &blockEdges,
                                           std::vector<std::vector<int>> &edgeCount,
                                           unsigned int n_block){

  return 0;
}
