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
		};

		~PDBarycenter(){};


		std::vector<std::vector<matchingTuple>> execute(std::vector<diagramTuple>& barycenter);
			
		void setBidderDiagrams();
		dataType enrichCurrentBidderDiagrams(dataType previous_min_persistence, dataType min_persistence, std::vector<dataType> initial_diagonal_prices);
		void setInitialBarycenter(dataType min_persistence);
		dataType getMaxPersistence();
		dataType getLowestPersistence();
		std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> getKDTree();
		dataType updateBarycenter(std::vector<std::vector<matchingTuple>>& matchings);
		bool hasBarycenterConverged(std::vector<std::vector<matchingTuple>>& matchings, std::vector<std::vector<matchingTuple>>& previous_matchings);
		
		bool is_matching_stable();
			
		dataType getEpsilon(dataType rho);
		dataType getRho(dataType epsilon);
			
			
		inline int setDiagram(int idx, void* data){
			if(idx < numberOfInputs_){
			inputData_[idx] = data;
			}
			else{
			return -1;
			}
			return 0;
		}

		inline int setNumberOfInputs(int numberOfInputs){
			numberOfInputs_ = numberOfInputs;
			inputData_ = (void **) malloc(numberOfInputs*sizeof(void *));
			for(int i=0 ; i<numberOfInputs ; i++){
				inputData_[i] = NULL;
			}
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
      void**                inputData_; //TODO : std::vector<void*>
      int                   threadNumber_;
	  bool                  use_progressive_;
	  double                time_limit_;
      
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
