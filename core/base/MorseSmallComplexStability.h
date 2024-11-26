#pragma once

// ttk common includes
#include <Debug.h>
#include<array>
#include <vector>
#include<utility>


namespace ttk {

  class MorseSmallComplexStability : virtual public Debug {

  public:
    MorseSmallComplexStability();
    
    int computeGraphData ( std::vector<std::vector<std::array<double, 3>>> &blockCoords,
                  std::vector<std::vector<double>> &blockSfValues,
                  std::vector<std::vector<std::pair<int, int>>> &blockEdges,
                  std::vector<std::vector<int>> &edgeCount,
                  unsigned int n_block);


  };
} // namespace ttk

