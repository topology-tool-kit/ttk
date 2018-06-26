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
#include 				  <limits>
#include				  <PDBarycenter.h>

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
			threadNumber_ = 1;
		};

		~PersistenceDiagramsBarycenter(){
			if(inputData_)
				free(inputData_);
		};


		int execute();
			
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
		
		inline void setWasserstein(const std::string &wasserstein){
			wasserstein_ = (wasserstein == "inf") ? -1 : stoi(wasserstein);
		}
		
		inline void setThreadNumber(const int &ThreadNumber){
			threadNumber_ = ThreadNumber;
		}
		
		template<typename type>
		static type abs(const type var) {
			return (var >= 0) ? var : -var;
		}



    protected:
	  int 					wasserstein_;
	  double                geometrical_factor_; // TODO include it in barycenter
	  
      int                   numberOfInputs_;
      void**                inputData_; //TODO : std::vector<void*>
      int 					threadNumber_;
      
      int points_added_;
	  int points_deleted_;
      
      std::vector<std::vector<dataType>>      all_matchings_;
 	  std::vector<std::vector<dataType>>      all_old_matchings_;
      std::vector<BidderDiagram<dataType>>    bidder_diagrams_;
      std::vector<GoodDiagram<dataType>>	  barycenter_goods_;
  };
  
  
template <typename dataType> 
int PersistenceDiagramsBarycenter<dataType>::execute(){
	Timer t;
	{
	
	std::vector<std::vector<diagramTuple>*> data_min(numberOfInputs_);
	std::vector<std::vector<diagramTuple>*> data_sad(numberOfInputs_);
	std::vector<std::vector<diagramTuple>*> data_max(numberOfInputs_);
	
	std::vector<std::vector<int>> data_min_idx(numberOfInputs_);
	std::vector<std::vector<int>> data_sad_idx(numberOfInputs_);
	std::vector<std::vector<int>> data_max_idx(numberOfInputs_);
	
	bool do_min = false;
	bool do_sad = false;
	bool do_max = false;
	
	// Create diagrams for min, saddle and max persistence pairs
	for(int i=0; i<numberOfInputs_; i++){
		data_min[i] = new std::vector<diagramTuple>;
		data_sad[i] = new std::vector<diagramTuple>;
		data_max[i] = new std::vector<diagramTuple>;
		std::vector<diagramTuple>* CTDiagram = static_cast<std::vector<diagramTuple>*>(inputData_[i]);
		
		for(int j=0; j<(int) CTDiagram->size(); ++j){
			diagramTuple t = CTDiagram->at(j);
			
			BNodeType nt1 = std::get<1>(t);
			BNodeType nt2 = std::get<3>(t);
			
			//dataType dt = std::get<4>(t);
			//if (abs<dataType>(dt) < zeroThresh) continue;

			if (nt1 == BLocalMin && nt2 == BLocalMax) {
				data_max[i]->push_back(t);
				data_max_idx[i].push_back(j);
				do_max = true;
			}
			else {
				if (nt1 == BLocalMax || nt2 == BLocalMax) {
					data_max[i]->push_back(t);
					data_max_idx[i].push_back(j);
					do_max = true;
				}
				if (nt1 == BLocalMin || nt2 == BLocalMin) {
					data_min[i]->push_back(t);
					data_min_idx[i].push_back(j);
					do_min = true;
				}
				if ((nt1 == BSaddle1 && nt2 == BSaddle2)
					|| (nt1 == BSaddle2 && nt2 == BSaddle1)) {
					data_sad[i]->push_back(t);
					data_sad_idx[i].push_back(j);
					do_sad = true;
				}
			}
		}
	}
	
	/*omp_set_num_threads(1);
	#ifdef TTK_ENABLE_OPENMP
	#pragma omp parallel sections
	#endif
	{
		#ifdef TTK_ENABLE_OPENMP
		#pragma omp section
		#endif
		{*/
			if(do_min){
				std::cout << "Computing Minima barycenter..."<<std::endl;
				PDBarycenter<dataType> bary_min = PDBarycenter<dataType>();
				bary_min.setThreadNumber(threadNumber_);
				bary_min.setWasserstein(wasserstein_);
				bary_min.setNumberOfInputs(numberOfInputs_);
				for(int i=0; i<numberOfInputs_; i++){
					bary_min.setDiagram(i, data_min[i]);
				}
				bary_min.execute();
			}
		/*}
		
		#ifdef TTK_ENABLE_OPENMP
		#pragma omp section
		#endif
		{*/
			if(do_sad){
				std::cout << "Computing Saddles barycenter..."<<std::endl;
				PDBarycenter<dataType> bary_sad = PDBarycenter<dataType>();
				bary_sad.setThreadNumber(threadNumber_);
				bary_sad.setWasserstein(wasserstein_);
				bary_sad.setNumberOfInputs(numberOfInputs_);
				for(int i=0; i<numberOfInputs_; i++){
					bary_sad.setDiagram(i, data_sad[i]);
				}
				bary_sad.execute();
			}
		/*}
		
		#ifdef TTK_ENABLE_OPENMP
		#pragma omp section
		#endif
		{*/
			if(do_max){
				std::cout << "Computing Maxima barycenter..."<<std::endl;
				PDBarycenter<dataType> bary_max = PDBarycenter<dataType>();
				bary_max.setThreadNumber(threadNumber_);
				bary_max.setWasserstein(wasserstein_);
				bary_max.setNumberOfInputs(numberOfInputs_);
				for(int i=0; i<numberOfInputs_; i++){
					bary_max.setDiagram(i, data_max[i]);
				}
				bary_max.execute();
			}
		//}
	//}
	
	for(int i=0; i<numberOfInputs_; i++){
		delete data_min[i];
		delete data_sad[i];
		delete data_max[i];
	}
	
	std::stringstream msg;
	msg << "[PersistenceDiagramsBarycenter] processed in "
		<< t.getElapsedTime() << " s. (" << threadNumber_
		<< " thread(s))."
		<< std::endl;
	dMsg(std::cout, msg.str(), timeMsg);
	}
	
	
	return 0;
}

}
  


// if the package is a pure template class, uncomment the following line

#include <PDBarycenterImpl.h>
#include <PDBarycenter.h>
#endif 
