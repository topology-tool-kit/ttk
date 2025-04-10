/// \ingroup base
/// \class ttk::DimensionReductionMetrics
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date September 2024.
///
/// \brief TTK base class that computes different scores for the quality
/// of a dimension reduction
///
/// This module defines the %DimensionReductionMetrics class that computes
/// different scores for the quality of a dimension reduction
///
/// \sa ttkDimensionReductionMetrics.cpp %for a usage example.

#pragma once

// ttk common includes
#include <Debug.h>

namespace ttk {

  class DimensionReductionMetrics : virtual public Debug {

  public:
    DimensionReductionMetrics();

    /**
     * @brief Main entry point
     *
     * @param[in] input First point cloud (typically the HD input)
     * @param[in] latent Second point cloud (typically the LD representation)
     */
    void execute(std::vector<std::vector<double>> const &input,
                 std::vector<std::vector<double>> const &latent);

    struct Metrics {
      double w0, w1, ta, lc, rmse, trust, cont, lcmc, mrreh, mrrel;
    };
    Metrics get() const {
      return {w0_, w1_, ta_, lc_, rmse_, trust_, cont_, lcmc_, mrreh_, mrrel_};
    }

  protected:
    /** "p" parameter used for the p-Wasserstein distances */
    double Wasserstein{2.};

    /** size of the triplet sample for the triplet accuracy computation; if set
     * to -1, all the N(N-1)(N-2)/6 triplets are used */
    int SampleSize{-1};

    /** size "K" of the K-neighborhood used for the rank-based measures
     * (trustworthiness, continuity, LCMC, MRRE) */
    unsigned NeighborhoodSize{10};

    /** p-Wasserstein distance between the 0-dimensional persistence diagrams in
     * both space (to the power p) */
    double w0_{0.};

    /** p-Wasserstein distance between the 1-dimensional persistence diagrams in
     * both space (to the power p) */
    double w1_{0.};

    /** Triplet accuracy between the input and the representation, i.e. the
     * percentage of triplets whose distances in both spaces have the same
     * relative order */
    double ta_{0.};

    /** Linear correlation of pairwise distances between the input and the
     * representation */
    double lc_{0.};

    /** Root mean squared error between distance matrices of the input and the
     * representation */
    double rmse_{0.};

    /** Trustworthiness is penalized when neighbors in the representation are
     * not neighbors in the input */
    double trust_{0.};

    /** Continuity is penalized when neighbors in the input are not neighbors
     * in the representation */
    double cont_{0.};

    /** Local continuity meta criterion translates the similarity of
     * neighborhoods in the input and the representation */
    double lcmc_{0.};

    /** Mean relative rank error with respect to the ranks in the
     * representation */
    double mrreh_{0.};

    /** Mean relative rank error with respect to the ranks in the input */
    double mrrel_{0.};

  private:
    unsigned n_;
    unsigned dimHigh_;
    unsigned dimLow_;
    std::vector<double> inputCompressedDistanceMatrix_;
    std::vector<double> latentCompressedDistanceMatrix_;

    inline double inputDM(unsigned i, unsigned j) const {
      if(i == j)
        return 0.;
      else
        return inputCompressedDistanceMatrix_[std::max(i, j)
                                                * (std::max(i, j) - 1) / 2
                                              + std::min(i, j)];
    }

    inline double latentDM(unsigned i, unsigned j) const {
      if(i == j)
        return 0.;
      else
        return latentCompressedDistanceMatrix_[std::max(i, j)
                                                 * (std::max(i, j) - 1) / 2
                                               + std::min(i, j)];
    }

    void computeTopologicalMetrics();
    void computeTripletAccuracy();
    void computePairwiseDistanceBasedMetrics();
    void computeRankBasedMetrics();

    bool tripletOrderPreserved(int i, int j, int k) const;

  }; // DimensionReductionMetrics class

} // namespace ttk
