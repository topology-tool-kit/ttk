#ifndef _PDBARYCENTER_H
#define _PDBARYCENTER_H

#include 				  <Auction.h>
#include 				  <KDTree.h>
#include 				  <limits>
#include 				  <PersistenceDiagramsBarycenter.h>

using namespace std;
using namespace ttk;

namespace ttk{
  template<typename dataType>
  class PDBarycenter : public Debug{

	public:

		PDBarycenter(){
			wasserstein_ = 2;
			geometrical_factor_ = 1;
			threadNumber_ = 1;
			use_progressive_ = true;
			time_limit_ = std::numeric_limits<double>::max();
			epsilon_min_ = 1e-5;
		};

		~PDBarycenter(){};


		std::vector<std::vector<matchingTuple>> execute(std::vector<diagramTuple>& barycenter);
			
		void setBidderDiagrams();
		dataType enrichCurrentBidderDiagrams(dataType previous_min_persistence, dataType min_persistence, std::vector<dataType> initial_diagonal_prices, std::vector<dataType> initial_prices, int min_points_to_add, bool add_points_to_barycenter=true);
		void setInitialBarycenter(dataType min_persistence);
		dataType getMaxPersistence();
		dataType getLowestPersistence();
		dataType getMinimalPrice(int i);
		std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> getKDTree();
		void runMatching(dataType* total_cost, dataType epsilon, std::vector<int> sizes, KDTree<dataType>* kdt, std::vector<KDTree<dataType>*>* correspondance_kdt_map, std::vector<dataType>* min_diag_price, 	std::vector<dataType>* min_price, std::vector<std::vector<matchingTuple>>* all_matchings);
		dataType updateBarycenter(std::vector<std::vector<matchingTuple>>& matchings);
		bool hasBarycenterConverged(std::vector<std::vector<matchingTuple>>& matchings, std::vector<std::vector<matchingTuple>>& previous_matchings);
		std::vector<std::vector<matchingTuple>> correctMatchings(std::vector<std::vector<matchingTuple>> previous_matchings);
		
		
		bool is_matching_stable();
			
		dataType getEpsilon(dataType rho);
		dataType getRho(dataType epsilon);
			
			
// 		inline int setDiagram(int idx, void* data){
// 			if(idx < numberOfInputs_){
// 			inputData_[idx] = data;
// 			}
// 			else{
// 			return -1;
// 			}
// 			return 0;
// 		}
		inline int setDiagrams(std::vector<std::vector<diagramTuple> > *data){
			inputDiagrams_ = data;
			return 0;
		}

		inline int setNumberOfInputs(int numberOfInputs){
			numberOfInputs_ = numberOfInputs;
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
		
		inline void setTimeLimit(const double time_limit){
			time_limit_ = time_limit;
		}
		
		inline void setGeometricalFactor(const double geometrical_factor){
			geometrical_factor_ = geometrical_factor;
		}
		
		inline void setCurrentBidders(std::vector<BidderDiagram<dataType>>& diagrams){
			current_bidder_diagrams_ = diagrams;
		}
		
		inline void setCurrentBarycenter(std::vector<GoodDiagram<dataType>>& barycenters){
			barycenter_goods_ = barycenters;
		}
		
		inline std::vector<BidderDiagram<dataType>>& getCurrentBidders(){
			return current_bidder_diagrams_;
		}
		
		inline std::vector<GoodDiagram<dataType>>& getCurrentBarycenter(){
			return barycenter_goods_;
		}
		
		inline void setDiagramType(const int &diagramType){
			diagramType_ = diagramType;
			if(diagramType_==0){
				nt1_ = BLocalMin;
				nt2_ = BSaddle1;
			}
			else if(diagramType_==1){
				nt1_ = BSaddle1;
				nt2_ = BSaddle2;
			}
			else{
				nt1_ = BSaddle2;
				nt2_ = BLocalMax;
			}
		}
		
		template<typename type>
		static type abs(const type var) {
			return (var >= 0) ? var : -var;
		}



    protected:
	  int 					wasserstein_;
	  double                geometrical_factor_;
	  int					diagramType_;
	  BNodeType 			nt1_;
	  BNodeType 			nt2_;
	  
      int                   numberOfInputs_;
      int                   threadNumber_;
	  bool                  use_progressive_;
	  double                time_limit_;
	  float                 epsilon_min_;
	  std::vector<std::vector<diagramTuple>> *inputDiagrams_;
      
      int points_added_;
	  int points_deleted_;
      
      std::vector<std::vector<dataType>>      all_matchings_;
 	  std::vector<std::vector<dataType>>      all_old_matchings_;
      std::vector<BidderDiagram<dataType>>    bidder_diagrams_;
	  std::vector<BidderDiagram<dataType>>    current_bidder_diagrams_;
	  std::vector<std::vector<int>>           current_bidder_ids_;
      std::vector<GoodDiagram<dataType>>	  barycenter_goods_;
  };
}


#include <PDBarycenterImpl.h>
#endif
