/// \ingroup base
/// \class ttk::PersistenceDiagramsBarycenter
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \date August 2016.
///
/// \brief TTK processing package that takes an input ensemble data set 
/// (represented by a list of scalar fields) and which computes various 
/// vertexwise statistics (PDF estimation, bounds, moments, etc.)
///
/// \sa ttkPersistenceDiagramsBarycenter.cpp %for a usage example.

#ifndef _PERSISTENCEDIAGRAMSBARYCENTER_H
#define _PERSISTENCEDIAGRAMSBARYCENTER_H



#ifndef diagramTuple
#define diagramTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::NodeType, ttk::ftm::idVertex, \
  ttk::ftm::NodeType, dataType, ttk::ftm::idVertex, \
  dataType, float, float, float, dataType, float, float, float>
#endif


#ifndef BNodeType
#define BNodeType ttk::ftm::NodeType
#define BLocalMax ttk::ftm::NodeType::Local_maximum
#define BLocalMin ttk::ftm::NodeType::Local_minimum
#define BSaddle1  ttk::ftm::NodeType::Saddle1
#define BSaddle2  ttk::ftm::NodeType::Saddle2
#define BIdVertex ttk::ftm::idVertex
#endif


// base code includes
#include                  <PersistenceDiagramsBarycenter.cpp>
#include                  <Wrapper.h>
#include                  <PersistenceDiagram.h>
#include 				  <Auction.h>
#include 				  <KDTree.h>

using namespace std;
using namespace ttk;

namespace ttk{
  template<typename dataType>
  class PersistenceDiagramsBarycenter : public Debug{

	public:

		PersistenceDiagramsBarycenter(){
			wasserstein_ = 2;
			geometrical_factor_ = 1;
			inputData_ = NULL;
			numberOfInputs_ = 0;
		};

		~PersistenceDiagramsBarycenter(){
			if(inputData_)
				free(inputData_);
		};


		int execute();
			
		void setBidderDiagrams();
		void setInitialBarycenter();
		std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> getKDTree();
			
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
			if(inputData_)
			free(inputData_);
			inputData_ = (void **) malloc(numberOfInputs*sizeof(void *));
			for(int i=0 ; i<numberOfInputs ; i++){
			inputData_[i] = NULL;
			}
			return 0;
		}




    protected:
	  int 					wasserstein_;
	  double                geometrical_factor_;
	  
      int                   numberOfInputs_;
      void**                inputData_; //TODO : std::vector<void*>
      
      std::vector<std::vector<dataType>>      all_matchings_;
 	  std::vector<std::vector<dataType>>      all_old_matchings_;
      std::vector<BidderDiagram<dataType>>    bidder_diagrams_;
      std::vector<GoodDiagram<dataType>>	  barycenter_goods_;
  };
}

// if the package is a pure template class, uncomment the following line

#include <PDBarycenterImpl.h>
#endif // PERSISTENCEDIAGRAMSBARYCENTER_H
