/// \ingroup base
/// \class ttk::PersistenceDiagramsClustering
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \date August 2016.
///
/// \brief TTK processing package that takes an input ensemble data set
/// (represented by a list of scalar fields) and which computes various
/// vertexwise statistics (PDF estimation, bounds, moments, etc.)
///
/// \sa ttkPersistenceDiagramsClustering.cpp %for a usage example.

#ifndef _PERSISTENCEDIAGRAMSCLUSTERING_H
#define _PERSISTENCEDIAGRAMSCLUSTERING_H



#ifndef diagramTuple
#define diagramTuple std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId, \
  ttk::CriticalType, dataType, ttk::SimplexId, \
  dataType, float, float, float, dataType, float, float, float>
#endif


#ifndef BNodeType
#define BNodeType ttk::CriticalType
#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1  ttk::CriticalType::Saddle1
#define BSaddle2  ttk::CriticalType::Saddle2
#define BIdVertex ttk::SimplexId
#endif


// base code includes
#include                  <PersistenceDiagramsClustering.cpp>
#include                  <Wrapper.h>
#include                  <PersistenceDiagram.h>
#include 				  <PersistenceDiagramsBarycenter.h>
#include 				  <limits>
#include 				  <PDClusteringImpl.h>
#include				  <PDClustering.h>


using namespace std;
using namespace ttk;

namespace ttk{
  template<typename dataType>
  class PersistenceDiagramsClustering : public Debug{

	public:

		PersistenceDiagramsClustering(){
			wasserstein_ = 2;
			inputData_ = NULL;
			numberOfInputs_ = 0;
			threadNumber_ = 1;
		};

		~PersistenceDiagramsClustering(){};


		std::vector<std::vector<matchingTuple> >
		execute(std::vector<diagramTuple>* barycenter);

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
    inline void setLambda(const double lambda){
      lambda_ = lambda;
    }

		inline void setTimeLimit(const double time_limit){
			time_limit_ = time_limit;
		}

		inline void setUseKmeansppInit(const bool UseKmeansppInit){
			use_kmeanspp_ = UseKmeansppInit;
		}

		inline void setUseAccelerated(const bool UseAccelerated){
			use_accelerated_ = UseAccelerated;
		}

		inline void setNumberOfClusters(const int NumberOfClusters){
			n_clusters_ = NumberOfClusters;
		}
    inline void setDeterministic(const bool deterministic){
			deterministic_ = deterministic;
		}
    inline void setPairTypeClustering(const int pairTypeClustering){
			pairTypeClustering_ = pairTypeClustering;
		}
    inline void setDebugLevel(const int debugLevel){
      debugLevel_ = debugLevel;
    }
		template<typename type>
		static type abs(const type var) {
			return (var >= 0) ? var : -var;
		}



    protected:
      // Critical pairs used for clustering
      // 0:min-saddles ; 1:saddles-saddles ; 2:sad-max ; else : all

      int                 debugLevel_;
      int         pairTypeClustering_;
      bool        deterministic_;
	  int 					wasserstein_;
	  int 					n_clusters_;

      int                   numberOfInputs_;
      void*                 inputData_; //TODO : std::vector<void*>
      int 					threadNumber_;
	  bool                  use_progressive_;
	  bool                  use_accelerated_;
	  bool                  use_kmeanspp_;
	  double                alpha_;
    double                lambda_;
	  double                time_limit_;

      int points_added_;
	  int points_deleted_;

      std::vector<std::vector<dataType>>      all_matchings_;
 	  std::vector<std::vector<dataType>>      all_old_matchings_;
      std::vector<BidderDiagram<dataType>>    bidder_diagrams_;
      std::vector<GoodDiagram<dataType>>	  barycenter_goods_;
  };


template <typename dataType>
  std::vector<std::vector<matchingTuple>>
    PersistenceDiagramsClustering<dataType>::execute(
      std::vector<diagramTuple>* barycenter){

	std::cout<< "Launching execute..." << std::endl;
	Timer t;
	{
	std::vector<std::vector<diagramTuple> > *intermediateDiagrams = (std::vector<std::vector<diagramTuple> > *) inputData_;
	std::cout<< "Number of diagrams : " <<  numberOfInputs_ << std::endl;
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

  switch(pairTypeClustering_){
  case(0):
  std::cout << "[ttkPersistenceDiagramsClustering] Only MIN-SAD Pairs" << '\n';
    do_max = false;
    do_sad = false;
    break;
  case(1):
    std::cout << "[ttkPersistenceDiagramsClustering] Only SAD-SAD Pairs" << '\n';
    do_max = false;
    do_min = false;
    break;
  case(2):
  std::cout << "[ttkPersistenceDiagramsClustering] Only SAD-MAX Pairs" << '\n';
    do_min = false;
    do_sad = false;
    break;
  default:
  std::cout << "[ttkPersistenceDiagramsClustering] All critical pairs : global clustering" << '\n';
  break;
  }

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
    KMeans.setDebugLevel(debugLevel_);
	KMeans.setKMeanspp(use_kmeanspp_);
	KMeans.setK(n_clusters_);

	KMeans.setDiagrams(&data_min, &data_sad, &data_max);
	KMeans.setDos(do_min, do_sad, do_max);
	KMeans.execute();



	std::stringstream msg;
	msg << "[PersistenceDiagramsClustering] processed in "
		<< t.getElapsedTime() << " s. (" << threadNumber_
		<< " thread(s))."
		<< std::endl;
	dMsg(std::cout, msg.str(), timeMsg);
	return all_matchings;
	}

}

}

#endif
