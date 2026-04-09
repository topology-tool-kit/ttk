/// \ingroup vtk
/// \class ttkSeparatrixStability
/// \author Thomas Daniel <thomas.daniel124@gmail.com>
/// \date July 2025.
///
/// \brief Compute the occurrence rate of separatrices in an ensemble.
///
/// This module takes as an input a set of 1-dimensional separatrices of the 
/// Morse-Smale complex for an ensemble dataset. It computes as an output, for 
/// each separatrix, its rate of occurrence in the ensemble (based on partial 
/// isomorphism computations).
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/molecularVibration/">Molecular
///   Vibration example</a>
///   
/// \b Related \b publication: \n
/// "BondMatcher: H-Bond Stability Analysis in Molecular Systems" \n
/// Thomas Daniel, Malgorzata Olejniczak, Julien Tierny \n
/// IEEE Transactions on Visualization and Computer Graphics \n
/// Proc. of IEEE VIS 2025.

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

    int buildOccurrenceArrays(
      const std::vector<GraphMatrixFull> &adjacencyMatrices,
      const std::vector<int> &separatrixCountForEachBlock,
      const std::vector<std::vector<std::array<double, 3>>> &coordsSource,
      const std::vector<std::vector<std::array<double, 3>>> &coordsDestination,
      const std::vector<std::vector<double>> &scalarsSource,
      const std::vector<std::vector<double>> &scalarsDestination,
      const bool &mergeEdgesOnSaddles,
      std::vector<std::vector<int>> &edgesOccurrencesForEachBlock,
      std::vector<std::vector<bool>> &isomorphismForEachBlock,
      std::vector<std::vector<std::vector<int>>>
        &matchingArrayForEachBlockSource,
      std::vector<std::vector<std::vector<int>>>
        &matchingArrayForEachBlockDestination,
      std::vector<std::vector<std::vector<int>>>
        &matchingArraySeparatrixForEachBlock);

    inline void setEpsilon(double e) {
      epsilon = e;
    }

    inline void setWeights(const double &px,
                           const double &py,
                           const double &pz,
                           const double &pf) {
      Px = px;
      Py = py;
      Pz = pz;
      Pf = pf;
    }

  private:
    int buildMatchingsWithOtherBlocks(
      const std::vector<std::vector<std::array<double, 3>>> &coords,
      const std::vector<std::vector<double>> &scalars,
      const int &block_id,
      std::vector<std::vector<MatchingType>> &matchings);

    void computeGraphMinor(const GraphMatrixFull &adjacencyMatrixFull,
                           GraphMatrixMinor &adjacencyMatrix);

    void buildCostMatrix(const std::vector<std::array<double, 3>> &coords1,
                         const std::vector<std::array<double, 3>> &coords2,
                         const std::vector<double> &scalars1,
                         const std::vector<double> &scalars2,
                         std::vector<std::vector<double>> &matrix);

    void assignmentSolver(std::vector<std::vector<double>> &costMatrix,
                          std::vector<ttk::MatchingType> &matching);

    int buildOccurrenceArraysMinor(
      const std::vector<GraphMatrixFull> &adjacencyMatrices,
      const int &n_separatrices,
      const std::vector<std::vector<std::array<double, 3>>> &coords,
      const std::vector<std::vector<double>> &scalars,
      const int &block_id,
      std::vector<int> &edgeOccurrences,
      std::vector<bool> &isIsomorphicWith,
      std::vector<std::vector<int>> &matchingArray,
      std::vector<std::vector<int>> &matchingArraySeparatrix);

    int buildOccurrenceArraysFull(
      const std::vector<GraphMatrixFull> &adjacencyMatrices,
      const int &n_separatrices,
      const std::vector<std::vector<std::array<double, 3>>> &coordsSource,
      const std::vector<std::vector<std::array<double, 3>>> &coordsDestination,
      const std::vector<std::vector<double>> &scalarsSource,
      const std::vector<std::vector<double>> &scalarsDestination,
      const int &block_id,
      std::vector<int> &edgesOccurrences,
      std::vector<bool> &isIsomorphicWith,
      std::vector<std::vector<int>> &matchingArraySource,
      std::vector<std::vector<int>> &matchingArrayDestination,
      std::vector<std::vector<int>> &matchingArraySeparatrix);

    double epsilon{};
    double Px{1};
    double Py{1};
    double Pz{1};
    double Pf{1};
  };

} // namespace ttk
