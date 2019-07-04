#ifndef _PDCLUSTERING_H
#define _PDCLUSTERING_H

#include 				  <Auction.h>
#include 				  <KDTree.h>
#include 				  <limits>
#include 				  <PersistenceDiagramsClustering.h>

using namespace std;
using namespace ttk;

namespace ttk{
  template<typename dataType>
  class PDClustering : public Debug{

	public:

		PDClustering(){
			wasserstein_ = 2;
			geometrical_factor_ = 1;
			threadNumber_ = 1;
			use_progressive_ = true;
            deterministic_ = true;
			time_limit_ = std::numeric_limits<double>::max();
			epsilon_min_ = 1e-5;
			epsilon_.resize(3);
			precision_criterion_ = false;
			precision_min_ = false;
			precision_sad_ = false;
			precision_max_ = false;
			cost_min_=0;
			cost_max_=0;
			cost_sad_=0;
		};

		~PDClustering(){};


        std::vector<int> execute(std::vector<std::vector<diagramTuple>>& final_centroids);

		dataType getMostPersistent(int type=-1);
		dataType getLessPersistent(int type=-1);
		std::vector<std::vector<dataType>> getMinDiagonalPrices();
		std::vector<std::vector<dataType>> getMinPrices();

		dataType computeDistance(BidderDiagram<dataType>& D1, BidderDiagram<dataType>& D2, dataType delta_lim);
		dataType computeDistance(BidderDiagram<dataType> D1, GoodDiagram<dataType> D2, dataType delta_lim);
		dataType computeDistance(BidderDiagram<dataType>* D1, GoodDiagram<dataType>* D2, dataType delta_lim);
		dataType computeDistance(GoodDiagram<dataType>& D1, GoodDiagram<dataType>& D2, dataType delta_lim);

		GoodDiagram<dataType> centroidWithZeroPrices(GoodDiagram<dataType> centroid);
		BidderDiagram<dataType> centroidToDiagram(GoodDiagram<dataType> centroid);
		GoodDiagram<dataType> diagramToCentroid(BidderDiagram<dataType> diagram);
		BidderDiagram<dataType> diagramWithZeroPrices(BidderDiagram<dataType> diagram);

		void setBidderDiagrams();
		void initializeEmptyClusters();
		void initializeCentroids();
		void initializeCentroidsKMeanspp();
		void initializeAcceleratedKMeans();
		void initializeBarycenterComputers();
		void printDistancesToFile();
		void printRealDistancesToFile();
		void printPricesToFile(int);
		dataType computeRealCost();

        std::vector<dataType> enrichCurrentBidderDiagrams(std::vector<dataType> previous_min_persistence, 
		        std::vector<dataType> min_persistence, 
		        std::vector<std::vector<dataType>> initial_diagonal_prices, 
		        std::vector<std::vector<dataType>> initial_off_diagonal_points, 
		        std::vector<int> min_points_to_add, 
		        bool add_points_to_barycenter);

		std::vector<std::vector<dataType>> getDistanceMatrix();
		void getCentroidDistanceMatrix();

		void updateClusters();
		void invertClusters();
		void invertInverseClusters();

		void acceleratedUpdateClusters();
                std::vector<dataType> updateCentroidsPosition(std::vector<std::vector<dataType>>* min_price, std::vector<std::vector<dataType>>* min_diag_price, std::vector<std::vector<std::vector<matchingTuple>>>& all_matchings);

        inline void resetDosToOriginalValues(){
            do_min_ = original_dos[0];
            do_sad_ = original_dos[1];
            do_max_ = original_dos[2];
        }
		inline int setDiagrams(std::vector<std::vector<diagramTuple> > *data_min, std::vector<std::vector<diagramTuple> > *data_saddle, std::vector<std::vector<diagramTuple> > *data_max){
			inputDiagramsMin_ = data_min;
			inputDiagramsSaddle_ = data_saddle;
			inputDiagramsMax_ = data_max;
			return 0;
		}

		inline int setDos(bool doMin, bool doSad, bool doMax){
			do_min_ = doMin;
			do_sad_ = doSad;
			do_max_ = doMax;
			
			original_dos.resize(3);

			original_dos[0]=do_min_;
			original_dos[1]=do_sad_;
			original_dos[2]=do_max_;
			return 0;
		}

		inline int setNumberOfInputs(int numberOfInputs){
			numberOfInputs_ = numberOfInputs;
			return 0;
		}

		inline int setK(const int k){
			k_ = k;
			return 0;
		}


		inline void setWasserstein(const int &wasserstein){
			wasserstein_ = wasserstein;
		}

		inline void setThreadNumber(const int &threadNumber){
			threadNumber_ = threadNumber;
		}

		inline void setUseProgressive(const bool use_progressive){
			use_progressive_ = use_progressive;
		}

		inline void setKMeanspp(const bool use_kmeanspp){
			use_kmeanspp_ = use_kmeanspp;
		}

		inline void setUseKDTree(const bool use_kdtree){
			use_kdtree_ = use_kdtree;
		}

		inline void setAccelerated(const bool use_accelerated){
			use_accelerated_ = use_accelerated;
		}

		inline void setTimeLimit(const double time_limit){
			time_limit_ = time_limit;
		}

		inline void setGeometricalFactor(const double geometrical_factor){
			geometrical_factor_ = geometrical_factor;
		}
    inline void setLambda(const double lambda){
			lambda_ = lambda;
		}
    inline void setDeterministic(const bool deterministic){
			deterministic_ = deterministic;
		}
    inline void setDebugLevel(const int debugLevel){
      debugLevel_ = debugLevel;
    }
  
    inline void printClustering(){
		for(int c=0; c<k_; ++c){
			std::cout<<"Cluster "<< c << " : [";
			for(unsigned int idx=0; idx<clustering_[c].size(); ++idx){
				if(idx==clustering_[c].size()-1){
					std::cout<< clustering_[c][idx]<< "]" <<std::endl;
				}
				else{
					std::cout<< clustering_[c][idx]<< ", ";
				}
			}
		}
    }

    inline void printOldClustering(){
		for(int c=0; c<k_; ++c){
			std::cout<<"Cluster "<< c << " : [";
			for(unsigned int idx=0; idx<old_clustering_[c].size(); ++idx){
				if(idx==old_clustering_[c].size()-1){
					std::cout<< old_clustering_[c][idx]<< "]" <<std::endl;
				}
				else{
					std::cout<< old_clustering_[c][idx]<< ", ";
				}
			}
		}
	}


		template<typename type>
		static type abs(const type var) {
			return (var >= 0) ? var : -var;
		}



    protected:

    std::vector<PDBarycenter<dataType>> barycenter_computer_min_;
    std::vector<PDBarycenter<dataType>> barycenter_computer_sad_;
    std::vector<PDBarycenter<dataType>> barycenter_computer_max_;

    bool 	barycenter_inputs_reset_flag;
    bool	precision_criterion_;
    bool	precision_max_;
    bool	precision_min_;
    bool	precision_sad_;
    bool           deterministic_;
	  int 					wasserstein_;
	  double                geometrical_factor_;
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of the persistence paired
    // lambda = 1 : extremum (min if pair min-sad, max if pair sad-max)
    // lambda = 0 : saddle (bad stability)
    // lambda = 1/2 : middle of the 2 critical points of the pair
    double                lambda_;

	  int 					k_;
      int debugLevel_;
      int                   numberOfInputs_;
      int                   threadNumber_;
	  bool                  use_progressive_;
	  bool                  use_accelerated_;
	  bool                  use_kmeanspp_;
	  bool                  use_kdtree_;
	  double                time_limit_;

	  dataType              epsilon_min_;
      std::vector<double>              epsilon_;
	  dataType              cost_;
	  dataType              cost_min_;
	  dataType              cost_sad_;
	  dataType              cost_max_;

	  std::vector<std::vector<diagramTuple>> *inputDiagramsMin_;
	  std::vector<std::vector<diagramTuple>> *inputDiagramsSaddle_;
	  std::vector<std::vector<diagramTuple>> *inputDiagramsMax_;

      std::vector<bool>                     original_dos;

	  bool                                    do_min_;
      std::vector<BidderDiagram<dataType>>    bidder_diagrams_min_;
	  std::vector<BidderDiagram<dataType>>    current_bidder_diagrams_min_;
      std::vector<GoodDiagram<dataType>>	  centroids_min_;
	  std::vector<GoodDiagram<dataType>>	  centroids_with_price_min_;

	  bool                                    do_sad_;
	  std::vector<BidderDiagram<dataType>>    bidder_diagrams_saddle_;
	  std::vector<BidderDiagram<dataType>>    current_bidder_diagrams_saddle_;
      std::vector<GoodDiagram<dataType>>	  centroids_saddle_;
	  std::vector<GoodDiagram<dataType>>	  centroids_with_price_saddle_;

	  bool                                    do_max_;
	  std::vector<BidderDiagram<dataType>>    bidder_diagrams_max_;
	  std::vector<BidderDiagram<dataType>>    current_bidder_diagrams_max_;
      std::vector<GoodDiagram<dataType>>	  centroids_max_;
	  std::vector<GoodDiagram<dataType>>	  centroids_with_price_max_;

	  std::vector<std::vector<int>>           clustering_;
	  std::vector<std::vector<int>>           old_clustering_;
	  std::vector<int>                        inv_clustering_;

	  std::vector<bool>                       r_;
	  std::vector<dataType>                   u_;
	  std::vector<std::vector<dataType>>      l_;
	  std::vector<std::vector<dataType>>      d_;
	  int                                     n_iterations_;
  };
}


#include <PDClusteringImpl.h>
#endif
