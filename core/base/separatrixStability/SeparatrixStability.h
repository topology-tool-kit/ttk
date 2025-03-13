#pragma once

#include <Debug.h>
#include <Triangulation.h>
#include <optional>
#include <utility>
#include <vector>

namespace ttk {

  using MatchingType = std::tuple<int, int, double>;

  class SeparatrixStability : virtual public Debug {

  public:
    using GraphMatrixMinor
      = std::vector<std::vector<std::vector<std::pair<int, int>>>>;
    using GraphMatrixFull = std::vector<std::vector<int>>;

    SeparatrixStability();

    int buildOccurenceArrays(
      const std::vector<GraphMatrixFull> &adjacencyMatrices,
      const std::vector<int> &separatrixCountForEachBlock ,
      const std::vector<std::vector<std::array<double, 3>>> &coordsSource,
      const std::vector<std::vector<std::array<double, 3>>> &coordsDestination,
      const bool &mergeEdgesOnSaddles, 
      std::vector<std::vector<int>> &edgesOccurencesForEachBlock,
      std::vector<std::vector<bool>> &isomorphismForEachBlock,
      std::vector<std::vector<std::vector<int>>> &matchingArrayForEachBlockSource,
      std::vector<std::vector<std::vector<int>>> &matchingArrayForEachBlockDestination,
      std::vector<std::vector<std::vector<int>>> &matchingArraySeparatrixForEachBlock);

    inline void setEpsilon(double e){
      epsilon=e;
    }

  private:

    int buildMatchingsWithOtherBlocks(
      const std::vector<std::vector<std::array<double, 3>>> &coords,
      const int &block_id,
      std::vector<std::vector<MatchingType>> &matchings);

    void computeGraphMinor(const GraphMatrixFull &adjacencyMatrixFull,
                            GraphMatrixMinor &adjacencyMatrix);

    void buildCostMatrix(const std::vector<std::array<double, 3>> &coords1,
                         const std::vector<std::array<double, 3>> &coords2,
                         std::vector<std::vector<double>> &matrix);

    void assignmentSolver(std::vector<std::vector<double>> &costMatrix,
                          std::vector<ttk::MatchingType> &matching);
    

    int buildOccurenceArraysMinor(
      const std::vector<GraphMatrixFull> &adjacencyMatrices,
      const int &n_separatrices,
      const std::vector<std::vector<std::array<double, 3>>> &coords,
      const int &block_id,
      std::vector<int> &edgeOccurences,
      std::vector<bool> &isIsomorphicWith,
      std::vector<std::vector<int>> &matchingArray,
      std::vector<std::vector<int>> &matchingArraySeparatrix);

    int buildOccurenceArraysFull(
      const std::vector<GraphMatrixFull> &adjacencyMatrices,
      const int &n_separatrices,
      const std::vector<std::vector<std::array<double, 3>>> &coordsSource,
      const std::vector<std::vector<std::array<double, 3>>> &coordsDestination,
      const int &block_id,
      std::vector<int> &edgesOccurences,
      std::vector<bool> &isIsomorphicWith,
      std::vector<std::vector<int>> &matchingArraySource,
      std::vector<std::vector<int>> &matchingArrayDestination,
      std::vector<std::vector<int>> &matchingArraySeparatrix);

    double epsilon{};
  };


} // namespace ttk
