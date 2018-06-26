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
		};

		~PDBarycenter(){};


		int execute();
			
		void setBidderDiagrams();
		void setInitialBarycenter();
		dataType getMaxPersistence();
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
		
		template<typename type>
		static type abs(const type var) {
			return (var >= 0) ? var : -var;
		}



    protected:
	  int 					wasserstein_;
	  double                geometrical_factor_;
	  
      int                   numberOfInputs_;
      void**                inputData_; //TODO : std::vector<void*>
      
      int points_added_;
	  int points_deleted_;
      
      std::vector<std::vector<dataType>>      all_matchings_;
 	  std::vector<std::vector<dataType>>      all_old_matchings_;
      std::vector<BidderDiagram<dataType>>    bidder_diagrams_;
      std::vector<GoodDiagram<dataType>>	  barycenter_goods_;
  };
}


#include <PDBarycenterImpl.h>
#endif
