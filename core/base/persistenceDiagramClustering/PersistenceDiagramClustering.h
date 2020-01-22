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

#ifndef _PERSISTENCEDIAGRAMSCLUSTERING_H
#define _PERSISTENCEDIAGRAMSCLUSTERING_H

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
//
#include <Wrapper.h>
//
#include <PersistenceDiagram.h>
//
#include <PersistenceDiagramBarycenter.h>
//
#include <limits>
//
#include <PDClustering.h>
//

using namespace std;
using namespace ttk;

namespace ttk {
  template <typename dataType>
  class PersistenceDiagramClustering : public Debug {

  public:
    PersistenceDiagramClustering() {
      wasserstein_ = 2;
      use_progressive_ = 1;
      use_kmeanspp_ = 0;
      use_accelerated_ = 0;
      numberOfInputs_ = 0;
      threadNumber_ = 1;
    };

    ~PersistenceDiagramClustering(){};

    std::vector<int> execute(
      std::vector<std::vector<diagramTuple>> &intermediateDiagrams,
      std::vector<std::vector<diagramTuple>> &centroids,
      std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings);

    inline int setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
      // 			if(inputData_)
      // 			free(inputData_);
      // 			inputData_ = (void **) malloc(numberOfInputs*sizeof(void *));
      // 			for(int i=0 ; i<numberOfInputs ; i++){
      // 			inputData_[i] = NULL;
      // 			}
      return 0;
    }

    inline void setWasserstein(const std::string &wasserstein) {
      wasserstein_ = (wasserstein == "inf") ? -1 : stoi(wasserstein);
    }

    inline void setThreadNumber(const int &ThreadNumber) {
      threadNumber_ = ThreadNumber;
    }

    inline void setUseProgressive(const bool use_progressive) {
      use_progressive_ = use_progressive;
    }

    inline void setAlpha(const double alpha) {
      alpha_ = alpha;
    }
    inline void setLambda(const double lambda) {
      lambda_ = lambda;
    }

    inline void setTimeLimit(const double time_limit) {
      time_limit_ = time_limit;
    }

    inline void setUseKmeansppInit(const bool UseKmeansppInit) {
      use_kmeanspp_ = UseKmeansppInit;
    }

    inline void setUseAccelerated(const bool UseAccelerated) {
      use_accelerated_ = UseAccelerated;
    }

    inline void setNumberOfClusters(const int NumberOfClusters) {
      n_clusters_ = NumberOfClusters;
    }
    inline void setForceUseOfAlgorithm(const bool forceUseOfAlgorithm) {
      forceUseOfAlgorithm_ = forceUseOfAlgorithm;
    }
    inline void setDeterministic(const bool deterministic) {
      deterministic_ = deterministic;
    }
    inline void setPairTypeClustering(const int pairTypeClustering) {
      pairTypeClustering_ = pairTypeClustering;
    }

    inline void setUseDeltaLim(const bool useDeltaLim) {
      useDeltaLim_ = useDeltaLim;
    }

    inline void setDistanceWritingOptions(const int distanceWritingOptions) {
      distanceWritingOptions_ = distanceWritingOptions;
    }

    inline void setDeltaLim(const double deltaLim) {
      deltaLim_ = deltaLim;
    }
    inline void setOutputDistanceMatrix(const bool arg) {
      outputDistanceMatrix_ = arg;
    }
    inline void setUseFullDiagrams(const bool arg) {
      useFullDiagrams_ = arg;
    }
    inline void setPerClusterDistanceMatrix(const bool arg) {
      perClusterDistanceMatrix_ = arg;
    }
    template <typename type>
    static type abs(const type var) {
      return (var >= 0) ? var : -var;
    }

    inline const std::vector<std::vector<double>> &&getDiagramsDistMat() {
      return std::move(diagramsDistMat_);
    }
    inline const std::vector<std::vector<double>> &&getCentroidsDistMat() {
      return std::move(centroidsDistMat_);
    }
    inline const std::vector<double> &&getDistanceToCentroid() {
      return std::move(distanceToCentroid_);
    }

  protected:
    // Critical pairs used for clustering
    // 0:min-saddles ; 1:saddles-saddles ; 2:sad-max ; else : all

    double deltaLim_;
    bool useDeltaLim_;
    int distanceWritingOptions_;
    int pairTypeClustering_;
    bool forceUseOfAlgorithm_;
    bool deterministic_;
    int wasserstein_;
    int n_clusters_;

    int numberOfInputs_;
    int threadNumber_;
    bool use_progressive_;
    bool use_accelerated_;
    bool use_kmeanspp_;
    double alpha_;
    double lambda_;
    double time_limit_;

    int points_added_;
    int points_deleted_;

    std::vector<BidderDiagram<dataType>> bidder_diagrams_;
    std::vector<GoodDiagram<dataType>> barycenter_goods_;

    bool outputDistanceMatrix_{false};
    bool useFullDiagrams_{false};
    bool perClusterDistanceMatrix_{false};
    std::vector<std::vector<double>> centroidsDistMat_{};
    std::vector<std::vector<double>> diagramsDistMat_{};
    std::vector<double> distanceToCentroid_{};
  };

  template <typename dataType>
  std::vector<int> PersistenceDiagramClustering<dataType>::execute(
    std::vector<std::vector<diagramTuple>> &intermediateDiagrams,
    std::vector<std::vector<diagramTuple>> &final_centroids,
    std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings) {

    Timer tm;
    {
      std::stringstream msg;
      msg << "[PersistenceDiagramClustering] Clustering " << numberOfInputs_
          << " diagrams in " << n_clusters_ << " cluster(s)." << std::endl;
      dMsg(std::cout, msg.str(), infoMsg);
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
            data_max[i].push_back(t);
            data_max_idx[i].push_back(j);
            do_max = true;
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
      switch(pairTypeClustering_) {
        case(0):
          msg << "[PersistenceDiagramClustering] Only MIN-SAD Pairs";
          do_max = false;
          do_sad = false;
          break;
        case(1):
          msg << "[PersistenceDiagramClustering] Only SAD-SAD Pairs";
          do_max = false;
          do_min = false;
          break;
        case(2):
          msg << "[PersistenceDiagramClustering] Only SAD-MAX Pairs";
          do_min = false;
          do_sad = false;
          break;
        default:
          msg << "[PersistenceDiagramClustering] All critical pairs: "
                 "global clustering";
          break;
      }
      msg << std::endl;
      dMsg(std::cout, msg.str(), advancedInfoMsg);
    }

    vector<vector<vector<vector<matchingTuple>>>>
      all_matchings_per_type_and_cluster;
    PDClustering<dataType> KMeans = PDClustering<dataType>();
    KMeans.setWasserstein(wasserstein_);
    KMeans.setThreadNumber(threadNumber_);
    KMeans.setNumberOfInputs(numberOfInputs_);
    KMeans.setUseProgressive(use_progressive_);
    KMeans.setAccelerated(use_accelerated_);
    KMeans.setUseKDTree(true);
    KMeans.setTimeLimit(time_limit_);
    KMeans.setGeometricalFactor(alpha_);
    KMeans.setLambda(lambda_);
    KMeans.setDeterministic(deterministic_);
    KMeans.setForceUseOfAlgorithm(forceUseOfAlgorithm_);
    KMeans.setDebugLevel(debugLevel_);
    KMeans.setDeltaLim(deltaLim_);
    KMeans.setUseDeltaLim(useDeltaLim_);
    KMeans.setDistanceWritingOptions(distanceWritingOptions_);
    KMeans.setKMeanspp(use_kmeanspp_);
    KMeans.setK(n_clusters_);
    KMeans.setDiagrams(&data_min, &data_sad, &data_max);
    KMeans.setDos(do_min, do_sad, do_max);
    KMeans.setOutputDistanceMatrix(outputDistanceMatrix_);
    KMeans.setUseFullDiagrams(useFullDiagrams_);
    KMeans.setPerClusterDistanceMatrix(perClusterDistanceMatrix_);
    inv_clustering
      = KMeans.execute(final_centroids, all_matchings_per_type_and_cluster);
    vector<vector<int>> centroids_sizes = KMeans.get_centroids_sizes();

    centroidsDistMat_ = KMeans.getCentroidsDistanceMatrix();
    diagramsDistMat_ = KMeans.getDiagramsDistanceMatrix();
    distanceToCentroid_ = KMeans.getDistanceToCentroid();

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
      if(debugLevel_ > 20) {
        cout << "id in cluster " << idxInCluster[j] << endl;
      }
    }

    all_matchings.resize(n_clusters_);
    for(int c = 0; c < n_clusters_; c++) {
      all_matchings[c].resize(numberOfInputs_);
    }
    for(int i = 0; i < numberOfInputs_; i++) {
      unsigned int c = inv_clustering[i];
      // cout<<"CLUSTER : "<<inv_clustering[i]<<endl;
      // LARGEST PAIR MUST BE FIRST
      if(do_max) {
        matchingTuple t
          = all_matchings_per_type_and_cluster[c][2][idxInCluster[i]][0];
        int bidder_id = std::get<0>(t);
        if(bidder_id >= 0 && bidder_id < (int)data_max[i].size()) {
          std::get<0>(t) = data_max_idx[i][bidder_id];
          if(std::get<1>(t) >= 0) {
            std::get<1>(t)
              = std::get<1>(t) + centroids_sizes[c][0] + centroids_sizes[c][1];
          } else {
            std::get<1>(t) = -1;
          }
          all_matchings[inv_clustering[i]][i].push_back(t);
        }
      }
      if(do_min) {
        for(unsigned int j = 0;
            j
            < all_matchings_per_type_and_cluster[c][0][idxInCluster[i]].size();
            j++) {
          matchingTuple t
            = all_matchings_per_type_and_cluster[c][0][idxInCluster[i]][j];
          int bidder_id = std::get<0>(t);
          if(bidder_id >= 0 && bidder_id < (int)data_min[i].size()) {
            std::get<0>(t) = data_min_idx[i][bidder_id];
            // cout<<" IDS :  "<<bidder_id<<" "<<std::get<0>(t)<<endl;
            if(std::get<1>(t) < 0) {
              std::get<1>(t) = -1;
            }
            all_matchings[inv_clustering[i]][i].push_back(t);
          }
        }
      }

      if(do_sad) {
        for(unsigned int j = 1;
            j
            < all_matchings_per_type_and_cluster[c][1][idxInCluster[i]].size();
            j++) {
          matchingTuple t
            = all_matchings_per_type_and_cluster[c][1][idxInCluster[i]][j];
          int bidder_id = std::get<0>(t);
          if(bidder_id >= 0 && bidder_id < (int)data_sad[i].size()) {
            std::get<0>(t) = data_sad_idx[i][bidder_id];
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
          if(bidder_id >= 0 && bidder_id < (int)data_max[i].size()) {
            std::get<0>(t) = data_max_idx[i][bidder_id];
            if(std::get<1>(t) >= 0) {
              std::get<1>(t) = std::get<1>(t) + centroids_sizes[c][0]
                               + centroids_sizes[c][1];
            } else {
              std::get<1>(t) = -1;
            }
            all_matchings[inv_clustering[i]][i].push_back(t);
          }
        }
      }
    }

    {
      std::stringstream msg;
      msg << "[PersistenceDiagramClustering] Processed in "
          << tm.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
          << std::endl;
      dMsg(std::cout, msg.str(), infoMsg);
    }
    return inv_clustering;
  }
} // namespace ttk

#endif
