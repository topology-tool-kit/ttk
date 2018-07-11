#ifndef _PDBARYCENTER_H
#define _PDBARYCENTER_H

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
			time_limit_ = std::numeric_limits<double>::max();
			epsilon_min_ = 1e-5;
		};

		~PDClustering(){};


		int execute();
		
		dataType getMostPersistent(int id_of_diagram, int diagram_type);
		
		dataType computeDistance(BidderDiagram<dataType> D1, BidderDiagram<dataType> D2, dataType delta_lim=0.0001);
		dataType computeDistance(BidderDiagram<dataType> D1, GoodDiagram<dataType> D2, dataType delta_lim=0.0001);
		dataType computeDistance(GoodDiagram<dataType> D1, GoodDiagram<dataType> D2, dataType delta_lim=0.0001);
			
		GoodDiagram<dataType> centroidWithZeroPrices(GoodDiagram<dataType> centroid);
		BidderDiagram<dataType> centroidToDiagram(GoodDiagram<dataType> centroid);
		GoodDiagram<dataType> diagramToCentroid(BidderDiagram<dataType> diagram);
		
		void initializeEmptyClusters();
		void initializeCentroids();
		void initializeCentroidsKMeanspp();
		void initializeAcceleratedKMeans();
		
		std::vector<std::vector<dataType>> getDistanceMatrix();
		std::vector<std::vector<dataType>> getCentroidDistanceMatrix();
		
		void updateClusters();
		void invertClusters();
		void invertInverseClusters();
		
		void acceleratedUpdateClusters();
		void updateCentroidsPosition();
		
		
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
		}

		inline int setNumberOfInputs(int numberOfInputs){
			numberOfInputs_ = numberOfInputs;
			return 0;
		}
		
		inline setK(const int k){
			k_ = k;
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
		
		inline void setAccelerated(const bool use_accelerated){
			use_accelerated_ = use_accelerated;
		}
		
		inline void setTimeLimit(const double time_limit){
			time_limit_ = time_limit;
		}
		
		inline void setGeometricalFactor(const double geometrical_factor){
			geometrical_factor_ = geometrical_factor;
		}
		
	
		template<typename type>
		static type abs(const type var) {
			return (var >= 0) ? var : -var;
		}



    protected:
	  int 					wasserstein_;
	  double                geometrical_factor_;
	  int 					k_;
	  
      int                   numberOfInputs_;
      int                   threadNumber_;
	  bool                  use_progressive_;
	  bool                  use_accelerated_;
	  double                time_limit_;
	  
	  dataType              epsilon_min_;
	  dataType              epsilon_;
	  
	  std::vector<std::vector<diagramTuple>> *inputDiagramsMin_;
	  std::vector<std::vector<diagramTuple>> *inputDiagramsSaddle_;
	  std::vector<std::vector<diagramTuple>> *inputDiagramsMax_;
      
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
	  std::vector<int>                        inv_clustering_;
	  
	  std::vector<bool>                       r_;
	  std::vector<dataType>                   u_;
	  std::vector<std::vector<dataType>>      l_;
	  std::vector<std::vector<dataType>>      d_;
  };
}


#include <PDClusteringImpl.h>
#endif
