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

		~PersistenceDiagramsBarycenter(){};


	std::vector<std::vector<matchingTuple> > 
      execute(std::vector<diagramTuple>* barycenter);

// 		inline int setDiagram(int idx, void* data){
// 			if(idx < numberOfInputs_){
// 			inputData_[idx] = data;
// 			}
// 			else{
// 			return -1;
// 			}
// 			return 0;
// 		}
		inline int setDiagrams(void *data){
		inputData_ = data;
		return 0;
    }

		inline int setNumberOfInputs(int numberOfInputs){
			numberOfInputs_ = numberOfInputs;
// 			if(inputData_)
// 			free(inputData_);
// 			inputData_ = (void **) malloc(numberOfInputs*sizeof(void *));
// 			for(int i=0 ; i<numberOfInputs ; i++){
// 			inputData_[i] = NULL;
// 			}
			return 0;
		}
		
		inline void setWasserstein(const std::string &wasserstein){
			wasserstein_ = (wasserstein == "inf") ? -1 : stoi(wasserstein);
		}
		
		inline void setThreadNumber(const int &ThreadNumber){
			threadNumber_ = ThreadNumber;
		}
		
		inline void setUseProgressive(const bool use_progressive){
			use_progressive_ = use_progressive;
		}
		
		inline void setAlpha(const double alpha){
			alpha_ = alpha;
		}
		
		inline void setTimeLimit(const double time_limit){
			time_limit_ = time_limit;
		}
		
		template<typename type>
		static type abs(const type var) {
			return (var >= 0) ? var : -var;
		}
		
		
		inline void setReinitPrices(const bool reinit_prices){
			reinit_prices_ = reinit_prices;
		}
		
		inline void setEpsilonDecreases(const bool epsilon_decreases){
			epsilon_decreases_ = epsilon_decreases;
		}
		
		inline void setEarlyStoppage(const bool early_stoppage){
			early_stoppage_ = early_stoppage;
		}



    protected:
	  int 					wasserstein_;
	  double                geometrical_factor_; // TODO include it in barycenter
	  
      int                   numberOfInputs_;
      void*                inputData_; //TODO : std::vector<void*>
      int 					threadNumber_;
	  bool                  use_progressive_;
	  double                alpha_;
	  double                time_limit_;
      
      
      int points_added_;
	  int points_deleted_;
      
      std::vector<std::vector<dataType>>      all_matchings_;
 	  std::vector<std::vector<dataType>>      all_old_matchings_;
      std::vector<BidderDiagram<dataType>>    bidder_diagrams_;
      std::vector<GoodDiagram<dataType>>	  barycenter_goods_;
	  
	  bool reinit_prices_;
	  bool epsilon_decreases_;
	  bool early_stoppage_;
  };
  
  
template <typename dataType> 
  std::vector<std::vector<matchingTuple>> 
    PersistenceDiagramsBarycenter<dataType>::execute(
      std::vector<diagramTuple>* barycenter){
	Timer t;
	{
	std::vector<std::vector<diagramTuple> > *intermediateDiagrams = 
		(std::vector<std::vector<diagramTuple> > *) inputData_;
	
	std::vector<std::vector<diagramTuple> > data_min(numberOfInputs_);
	std::vector<std::vector<diagramTuple> > data_sad(numberOfInputs_);
	std::vector<std::vector<diagramTuple> > data_max(numberOfInputs_);
	
	std::vector<std::vector<int>> data_min_idx(numberOfInputs_);
	std::vector<std::vector<int>> data_sad_idx(numberOfInputs_);
	std::vector<std::vector<int>> data_max_idx(numberOfInputs_);
	
	std::vector<std::vector<matchingTuple>> all_matchings(numberOfInputs_);
	
	bool do_min = false;
	bool do_sad = false;
	bool do_max = false;
	
	// Create diagrams for min, saddle and max persistence pairs
	for(int i=0; i<numberOfInputs_; i++){
    std::vector<diagramTuple>* CTDiagram = &((*intermediateDiagrams)[i]);
		
		for(int j=0; j<(int) CTDiagram->size(); ++j){
			diagramTuple t = CTDiagram->at(j);
			
			BNodeType nt1 = std::get<1>(t);
			BNodeType nt2 = std::get<3>(t);
			
			dataType dt = std::get<4>(t);
			//if (abs<dataType>(dt) < zeroThresh) continue;
			if(dt>0){
				if (nt1 == BLocalMin && nt2 == BLocalMax) {
					data_max[i].push_back(t);
					data_max_idx[i].push_back(j);
					do_max = true;
				}
				else {
					if (nt1 == BLocalMax || nt2 == BLocalMax) {
						data_max[i].push_back(t);
						data_max_idx[i].push_back(j);
						do_max = true;
					}
					if (nt1 == BLocalMin || nt2 == BLocalMin) {
						data_min[i].push_back(t);
						data_min_idx[i].push_back(j);
						do_min = true;
					}
					if ((nt1 == BSaddle1 && nt2 == BSaddle2)
						|| (nt1 == BSaddle2 && nt2 == BSaddle1)) {
						data_sad[i].push_back(t);
						data_sad_idx[i].push_back(j);
						do_sad = true;
					}
				}
			}
		}
	}
	
	std::vector<diagramTuple> barycenter_min;
	std::vector<diagramTuple> barycenter_sad;
	std::vector<diagramTuple> barycenter_max;
	
	std::vector<std::vector<matchingTuple>> 
    matching_min, matching_sad, matching_max;
	
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
				bary_min.setDiagramType(0);
				bary_min.setUseProgressive(use_progressive_);
				bary_min.setTimeLimit(time_limit_);
				bary_min.setGeometricalFactor(alpha_);
				bary_min.setEarlyStoppage(early_stoppage_);
				bary_min.setEpsilonDecreases(epsilon_decreases_);
				bary_min.setReinitPrices(reinit_prices_);
        bary_min.setDiagrams(&data_min);
				matching_min = bary_min.execute(barycenter_min);
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
				bary_sad.setDiagramType(1);
				bary_sad.setUseProgressive(use_progressive_);
				bary_sad.setTimeLimit(time_limit_);
				bary_sad.setGeometricalFactor(alpha_);
				bary_sad.setEarlyStoppage(early_stoppage_);
				bary_sad.setEpsilonDecreases(epsilon_decreases_);
				bary_sad.setReinitPrices(reinit_prices_);
        bary_sad.setDiagrams(&data_sad);
				matching_sad = bary_sad.execute(barycenter_sad);
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
				bary_max.setDiagramType(2);
				bary_max.setUseProgressive(use_progressive_);
				bary_max.setTimeLimit(time_limit_);
				bary_max.setGeometricalFactor(alpha_);
				bary_max.setEarlyStoppage(early_stoppage_);
				bary_max.setEpsilonDecreases(epsilon_decreases_);
				bary_max.setReinitPrices(reinit_prices_);
        bary_max.setDiagrams(&data_max);
				matching_max = bary_max.execute(barycenter_max);
			}
		//}
	//}
	
	// Reconstruct matchings
	for(int i=0; i<numberOfInputs_; i++){
		
		if(do_min){
			for(unsigned int j=0; j<matching_min[i].size(); j++){
				matchingTuple t = matching_min[i][j];
				int bidder_id = std::get<0>(t);
				std::get<0>(t) = data_min_idx[i][bidder_id];
				all_matchings[i].push_back(t);
			}
		}
		
		if(do_sad){
			for(unsigned int j=0; j<matching_sad[i].size(); j++){
				matchingTuple t = matching_sad[i][j];
				int bidder_id = std::get<0>(t);
				std::get<0>(t) = data_sad_idx[i][bidder_id];
				std::get<1>(t) = std::get<1>(t) + barycenter_min.size();
				all_matchings[i].push_back(t);
			}
		}
		
		if(do_max){
			for(unsigned int j=0; j<matching_max[i].size(); j++){
				matchingTuple t = matching_max[i][j];
				int bidder_id = std::get<0>(t);
				std::get<0>(t) = data_max_idx[i][bidder_id];
				std::get<1>(t) = std::get<1>(t) + barycenter_min.size() + barycenter_sad.size();
				all_matchings[i].push_back(t);
			}
		}
	}
	// Reconstruct barcenter
	for(unsigned int j=0; j<barycenter_min.size(); j++){
		diagramTuple dt = barycenter_min[j];
		barycenter->push_back(dt);
	}
	for(unsigned int j=0; j<barycenter_sad.size(); j++){
		diagramTuple dt = barycenter_sad[j];
		barycenter->push_back(dt);
	}
	for(unsigned int j=0; j<barycenter_max.size(); j++){
		diagramTuple dt = barycenter_max[j];
		barycenter->push_back(dt);
	}
	
	// Recreate 3D critical coordinates of barycentric points
	std::vector<int> number_of_matchings_for_point(barycenter->size());
	std::vector<float> cords_x1(barycenter->size());
	std::vector<float> cords_y1(barycenter->size());
	std::vector<float> cords_z1(barycenter->size());
	std::vector<float> cords_x2(barycenter->size());
	std::vector<float> cords_y2(barycenter->size());
	std::vector<float> cords_z2(barycenter->size());
	for(unsigned i=0; i<barycenter->size(); i++){
		number_of_matchings_for_point[i] = 0;
		cords_x1[i] = 0;
		cords_y1[i] = 0;
		cords_z1[i] = 0;
		cords_x2[i] = 0;
		cords_y2[i] = 0;
		cords_z2[i] = 0;
	}
	
	for(unsigned i=0; i<all_matchings.size(); i++){
    std::vector<diagramTuple>* CTDiagram = &((*intermediateDiagrams)[i]);
		for(unsigned j=0; j<all_matchings[i].size(); j++){
			matchingTuple t = all_matchings[i][j];
			int bidder_id = std::get<0>(t);
			int bary_id = std::get<1>(t);
			
			diagramTuple &bidder = CTDiagram->at(bidder_id);
			number_of_matchings_for_point[bary_id] +=1;
			cords_x1[bary_id] += std::get<7>(bidder);
			cords_y1[bary_id] += std::get<8>(bidder);
			cords_z1[bary_id] += std::get<9>(bidder);
			cords_x2[bary_id] += std::get<11>(bidder);
			cords_y2[bary_id] += std::get<12>(bidder);
			cords_z2[bary_id] += std::get<13>(bidder);
		}
	}
	
	for(unsigned i=0; i<barycenter->size(); i++){
		if(number_of_matchings_for_point[i]>0){
			std::get<7>(barycenter->at(i)) = cords_x1[i] / number_of_matchings_for_point[i];
			std::get<8>(barycenter->at(i)) = cords_y1[i] / number_of_matchings_for_point[i];
			std::get<9>(barycenter->at(i)) = cords_z1[i] / number_of_matchings_for_point[i];
			std::get<11>(barycenter->at(i)) = cords_x2[i] / number_of_matchings_for_point[i];
			std::get<12>(barycenter->at(i)) = cords_y2[i] / number_of_matchings_for_point[i];
			std::get<13>(barycenter->at(i)) = cords_z2[i] / number_of_matchings_for_point[i];
		}
	}
	
// 	for(int i=0; i<numberOfInputs_; i++){
// 		delete data_min[i];
// 		delete data_sad[i];
// 		delete data_max[i];
// 	}
	
	std::stringstream msg;
	msg << "[PersistenceDiagramsBarycenter] processed in "
		<< t.getElapsedTime() << " s. (" << threadNumber_
		<< " thread(s))."
		<< std::endl;
	dMsg(std::cout, msg.str(), timeMsg);
	return all_matchings;
	}
	
}

}
  


// if the package is a pure template class, uncomment the following line

#include <PDBarycenterImpl.h>
#include <PDBarycenter.h>
#endif 
