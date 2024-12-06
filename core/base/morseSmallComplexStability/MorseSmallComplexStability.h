#pragma once

#include <Debug.h>
#include <Triangulation.h>
#include <optional>
#include <utility>
#include <vector>

namespace ttk {

  using MatchingType = std::tuple<int, int, double>;

  class MorseSmallComplexStability : virtual public Debug {
  
  public:

    using GraphMatrixMinor = std::vector<std::vector<std::vector<std::pair<int, int>>>>;
    using GraphMatrixFull = std::vector<std::vector<int>>;

    MorseSmallComplexStability();

    int buildOccurenceArraysMinor(const std::vector<GraphMatrixMinor> &adjacencyMatrices, 
                                  const int &n_separatrices,
                                  const std::vector<std::vector<std::array<double, 3>>> &coords,
                                  const int &block_id,
                                  std::vector<int> &edgeOccurences);

    int buildOccurenceArraysFull(const std::vector<GraphMatrixFull> &adjacencyMatricesFull, 
                                  const int &n_separatrices, 
                                  const std::vector<std::vector<std::array<double, 3>>> &coordsSource,
                                  const std::vector<std::vector<std::array<double, 3>>> &coordsDestination,
                                  const int &block_id,
                                  std::vector<int> &edgesOccurences);

  private:

    int buildMatchingsWithOtherBlocks(const std::vector<std::vector<std::array<double, 3>>> &coords, 
                                      const int &block_id,
                                      std::vector<std::vector<MatchingType>> &matchings);

    void buildCostMatrix(const std::vector<std::array<double, 3>> &coords1,
                          const std::vector<std::array<double, 3>> &coords2,
                          std::vector<std::vector<double>> &matrix);

    void assignmentSolver(std::vector<std::vector<double>> &costMatrix,
                          std::vector<ttk::MatchingType> &matching);

  };

}
