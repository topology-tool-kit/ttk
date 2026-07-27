/// \ingroup base
/// \class ttk::ClusteringMetrics
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date January 2022
///
/// This module defines the %ClusteringMetrics class that computes both the NMI
/// (normalized mutual information) and the ARI (adjusted rand index) scores
/// between two clustering of the same points.
///
/// \b Related \b publication: \n
/// 'Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)'
/// Mathieu Pont, Jules Vidal and Julien Tierny.
/// IEEE Transactions on Visualization and Computer Graphics, 2023.
///
/// \sa ClusteringMetrics

#pragma once

#include <vector>
// ttk common includes
#include <Debug.h>

namespace ttk {

  /**
   * The ClusteringMetrics class provides methods to compute two scores
   * comparing two clustering of the same points. It can compute the NMI score
   * and the ARI score.
   */
  class ClusteringMetrics : virtual public Debug {

  public:
    ClusteringMetrics();

    int execute(const int *clustering1,
                const int *clustering2,
                const size_t n,
                double &nmiValue,
                double &ariValue) const;

    int
      computeContingencyTables(const int *clust1,
                               const int *clust2,
                               const size_t nPoint,
                               std::vector<std::vector<int>> &contingencyMatrix,
                               std::vector<int> &sumLin,
                               std::vector<int> &sumCol) const;

    // For single ARI or single NMI computing, the contingency matrix, as well
    // as the sums of each lines and each column must be provided as arguments.
    int computeARI(const std::vector<std::vector<int>> &contingencyMatrix,
                   const std::vector<int> &sumLin,
                   const std::vector<int> &sumCol,
                   const size_t nPoint,
                   double &ariValue) const;
    int computeNMI(const std::vector<std::vector<int>> &contingencyMatrix,
                   const std::vector<int> &sumLin,
                   const std::vector<int> &sumCol,
                   const size_t nPoint,
                   double &nmiValue) const;
  }; // ClusteringMetrics class

} // namespace ttk
