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

    using GraphMatrix = std::vector<std::vector<std::vector<std::pair<int, int>>>>;

    MorseSmallComplexStability();

    int buildOccurenceArrays(const std::vector<GraphMatrix> &adjacencyMatrices, 
                              const std::vector<int> &separatrixCountForEachBlock,
                              const std::vector<std::vector<std::array<double, 3>>> &coords,
                              std::vector<std::vector<int>> &occurencesMatrix);

  private:

    int buildVertexEquivalenceClasses(const std::vector<std::vector<std::array<double, 3>>> &coords, 
                                  std::vector<std::vector<int>> &classToVertexId);

    void buildCostMatrix(const std::vector<std::array<double, 3>> &coords1,
                          const std::vector<std::array<double, 3>> &coords2,
                          std::vector<std::vector<double>> &matrix);

    void assignmentSolver(std::vector<std::vector<double>> &costMatrix,
                          std::vector<ttk::MatchingType> &matching);

    void makePartition(const std::vector<std::vector<MatchingType>> &matchings, 
                        const int &n_blocks,
                        const int &n_points,
                        std::vector<std::vector<int>> &partitions);

  };

}
