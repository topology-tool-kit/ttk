/// \ingroup base
/// \class ttk::PersistenceDiagramClustering
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date September 2019
///
/// \brief TTK processing package for the computation of Wasserstein barycenters
/// and K-Means clusterings of a set of persistence diagrams.
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa ttkPersistenceDiagramClustering

#pragma once

#ifndef diagramTuple
#define diagramTuple                                                       \
  std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId,            \
             ttk::CriticalType, dataType, ttk::SimplexId, dataType, float, \
             float, float, dataType, float, float, float>
#endif

#ifndef BNodeType
#define BNodeType ttk::CriticalType
#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1 ttk::CriticalType::Saddle1
#define BSaddle2 ttk::CriticalType::Saddle2
#define BIdVertex ttk::SimplexId
#endif

// base code includes
#include <PDClustering.h>
#include <PersistenceDiagramBarycenter.h>

namespace ttk {

  class PersistenceDiagramClustering : virtual public Debug {

  public:
    PersistenceDiagramClustering() {
      this->setDebugMsgPrefix("PersistenceDiagramClustering");
    }

    ~PersistenceDiagramClustering() = default;

    template <class dataType>
    std::vector<int> execute(
      std::vector<std::vector<diagramTuple>> &intermediateDiagrams,
      std::vector<std::vector<diagramTuple>> &centroids,
      std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings);

    std::array<double, 3> getDistances() const {
      return this->distances;
    }

    template <class dataType>
    static dataType abs(const dataType var) {
      return (var >= 0) ? var : -var;
    }

  protected:
    // Critical pairs used for clustering
    // 0:min-saddles ; 1:saddles-saddles ; 2:sad-max ; else : all

    // distance results per pair type
    std::array<double, 3> distances{};

    int DistanceWritingOptions{0};
    int PairTypeClustering{-1};
    bool Deterministic{true};
    int WassersteinMetric{2};

    bool UseProgressive{true};

    bool ForceUseOfAlgorithm{false};
    bool UseInterruptible{true};
    double Alpha{1.0};
    bool UseAdditionalPrecision{false};
    double DeltaLim{0.01};
    double Lambda{1.0};
    double TimeLimit{999999};

    int NumberOfClusters{1};
    bool UseAccelerated{false};
    bool UseKmeansppInit{false};

    int points_added_;
    int points_deleted_;

    // std::vector<BidderDiagram<void>> bidder_diagrams_;
    // std::vector<GoodDiagram<void>> barycenter_goods_;
  };

  template <class dataType>
  std::vector<int> PersistenceDiagramClustering::execute(
    std::vector<std::vector<diagramTuple>> &intermediateDiagrams,
    std::vector<std::vector<diagramTuple>> &final_centroids,
    std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings) {

    const int numberOfInputs_ = intermediateDiagrams.size();
    Timer tm;
    {
      printMsg("Clustering " + std::to_string(numberOfInputs_) + " diagrams in "
               + std::to_string(NumberOfClusters) + " cluster(s).");
    }

    std::vector<std::vector<diagramTuple>> data_min(numberOfInputs_);
    std::vector<std::vector<diagramTuple>> data_sad(numberOfInputs_);
    std::vector<std::vector<diagramTuple>> data_max(numberOfInputs_);

    std::vector<std::vector<int>> data_min_idx(numberOfInputs_);
    std::vector<std::vector<int>> data_sad_idx(numberOfInputs_);
    std::vector<std::vector<int>> data_max_idx(numberOfInputs_);

    std::vector<int> inv_clustering(numberOfInputs_);

    bool do_min = false;
    bool do_sad = false;
    bool do_max = false;

    // Create diagrams for min, saddle and max persistence pairs
    for(int i = 0; i < numberOfInputs_; i++) {
      std::vector<diagramTuple> &CTDiagram = intermediateDiagrams[i];

      for(size_t j = 0; j < CTDiagram.size(); ++j) {
        diagramTuple t = CTDiagram[j];

        BNodeType nt1 = std::get<1>(t);
        BNodeType nt2 = std::get<3>(t);

        dataType dt = std::get<4>(t);
        // if (abs<dataType>(dt) < zeroThresh) continue;
        if(dt > 0) {
          if(nt1 == BLocalMin && nt2 == BLocalMax) {
            if(PairTypeClustering == 2) {
              data_max[i].push_back(t);
              data_max_idx[i].push_back(j);
              do_max = true;
            } else {
              data_min[i].push_back(t);
              data_min_idx[i].push_back(j);
              do_min = true;
            }
          } else {
            if(nt1 == BLocalMax || nt2 == BLocalMax) {
              data_max[i].push_back(t);
              data_max_idx[i].push_back(j);
              do_max = true;
            }
            if(nt1 == BLocalMin || nt2 == BLocalMin) {
              data_min[i].push_back(t);
              data_min_idx[i].push_back(j);
              do_min = true;
            }
            if((nt1 == BSaddle1 && nt2 == BSaddle2)
               || (nt1 == BSaddle2 && nt2 == BSaddle1)) {
              data_sad[i].push_back(t);
              data_sad_idx[i].push_back(j);
              do_sad = true;
            }
          }
        }
      }
    }

    {
      std::stringstream msg;
      switch(PairTypeClustering) {
        case(0):
          msg << "Only MIN-SAD Pairs";
          do_max = false;
          do_sad = false;
          break;
        case(1):
          msg << "Only SAD-SAD Pairs";
          do_max = false;
          do_min = false;
          break;
        case(2):
          msg << "Only SAD-MAX Pairs";
          do_min = false;
          do_sad = false;
          break;
        default:
          msg << "All critical pairs: "
                 "global clustering";
          break;
      }
      printMsg(msg.str());
    }

    std::vector<std::vector<std::vector<std::vector<matchingTuple>>>>
      all_matchings_per_type_and_cluster;
    PDClustering<dataType> KMeans = PDClustering<dataType>();
    KMeans.setNumberOfInputs(numberOfInputs_);
    KMeans.setWasserstein(WassersteinMetric);
    KMeans.setUseProgressive(UseProgressive);
    KMeans.setAccelerated(UseAccelerated);
    KMeans.setUseKDTree(true);
    KMeans.setTimeLimit(TimeLimit);
    KMeans.setGeometricalFactor(Alpha);
    KMeans.setLambda(Lambda);
    KMeans.setDeterministic(Deterministic);
    KMeans.setForceUseOfAlgorithm(ForceUseOfAlgorithm);
    KMeans.setDebugLevel(debugLevel_);
    KMeans.setDeltaLim(DeltaLim);
    KMeans.setUseDeltaLim(UseAdditionalPrecision);
    KMeans.setDistanceWritingOptions(DistanceWritingOptions);
    KMeans.setKMeanspp(UseKmeansppInit);
    KMeans.setK(NumberOfClusters);
    KMeans.setDiagrams(&data_min, &data_sad, &data_max);
    KMeans.setDos(do_min, do_sad, do_max);
    inv_clustering
      = KMeans.execute(final_centroids, all_matchings_per_type_and_cluster);
    std::vector<std::vector<int>> centroids_sizes
      = KMeans.get_centroids_sizes();

    this->distances = KMeans.getDistances();

    /// Reconstruct matchings
    //
    std::vector<int> cluster_size;
    std::vector<int> idxInCluster(numberOfInputs_);

    for(int j = 0; j < numberOfInputs_; ++j) {
      unsigned int c = inv_clustering[j];
      if(c + 1 > cluster_size.size()) {
        cluster_size.resize(c + 1);
        cluster_size[c] = 1;
        idxInCluster[j] = 0;
      } else {
        cluster_size[c]++;
        idxInCluster[j] = cluster_size[c] - 1;
      }
    }

    bool removeDuplicateGlobalPair = false;
    if(NumberOfClusters > 1 and do_min and do_max) {
      removeDuplicateGlobalPair = true;
    }

    all_matchings.resize(NumberOfClusters);
    for(int c = 0; c < NumberOfClusters; c++) {
      all_matchings[c].resize(numberOfInputs_);
      if(removeDuplicateGlobalPair) {
        centroids_sizes[c][0] -= 1;
      }
    }
    for(int i = 0; i < numberOfInputs_; i++) {
      unsigned int c = inv_clustering[i];

      if(do_min) {
        for(unsigned int j = 0;
            j
            < all_matchings_per_type_and_cluster[c][0][idxInCluster[i]].size();
            j++) {
          matchingTuple t
            = all_matchings_per_type_and_cluster[c][0][idxInCluster[i]][j];
          int bidder_id = std::get<0>(t);
          if(bidder_id < (int)data_min[i].size()) {
            if(bidder_id < 0) { // matching with a diagonal point
              std::get<0>(t) = -1;
            } else {
              std::get<0>(t) = data_min_idx[i][bidder_id];
            }
            // cout<<" IDS :  "<<bidder_id<<" "<<std::get<0>(t)<<endl;
            if(std::get<1>(t) < 0) {
              std::get<1>(t) = -1;
            }
            all_matchings[inv_clustering[i]][i].push_back(t);
          }
        }
      }

      if(do_sad) {
        for(unsigned int j = 0;
            j
            < all_matchings_per_type_and_cluster[c][1][idxInCluster[i]].size();
            j++) {
          matchingTuple t
            = all_matchings_per_type_and_cluster[c][1][idxInCluster[i]][j];
          int bidder_id = std::get<0>(t);
          if(bidder_id < (int)data_sad[i].size()) {
            if(bidder_id < 0) { // matching with a diagonal point
              std::get<0>(t) = -1;
            } else {
              std::get<0>(t) = data_sad_idx[i][bidder_id];
            }
            if(std::get<1>(t) >= 0) {
              std::get<1>(t) = std::get<1>(t) + centroids_sizes[c][0];
            } else {
              std::get<1>(t) = -1;
            }
            all_matchings[inv_clustering[i]][i].push_back(t);
          }
        }
      }

      if(do_max) {
        for(unsigned int j = 0;
            j
            < all_matchings_per_type_and_cluster[c][2][idxInCluster[i]].size();
            j++) {
          matchingTuple t
            = all_matchings_per_type_and_cluster[c][2][idxInCluster[i]][j];
          int bidder_id = std::get<0>(t);
          if(bidder_id < (int)data_max[i].size()) {
            if(bidder_id < 0) { // matching with a diagonal point
              std::get<0>(t) = -1;
            } else {
              std::get<0>(t) = data_max_idx[i][bidder_id];
            }
            // std::get<0>(t) = data_max_idx[i][bidder_id];
            if(std::get<1>(t) > 0) {
              std::get<1>(t) = std::get<1>(t) + centroids_sizes[c][0]
                               + centroids_sizes[c][1];
            } else if(std::get<1>(t) == 0) {
              if(!removeDuplicateGlobalPair) {
                std::get<1>(t) = std::get<1>(t) + centroids_sizes[c][0]
                                 + centroids_sizes[c][1];
              }
            } else {
              std::get<1>(t) = -1;
            }
            all_matchings[inv_clustering[i]][i].push_back(t);
          }
        }
      }
    }

    printMsg("Complete", 1, tm.getElapsedTime(), threadNumber_);
    return inv_clustering;
  }
} // namespace ttk
