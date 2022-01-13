/// \ingroup base
/// \class ttk::PersistenceDiagramBarycenter
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date September 2019
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa PersistenceDiagramClustering

#pragma once

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
#include <KDTree.h>
#include <PDBarycenter.h>
#include <PersistenceDiagramAuction.h>
#include <Wrapper.h>

namespace ttk {
  template <typename dataType>
  class PersistenceDiagramBarycenter : public Debug {

  public:
    PersistenceDiagramBarycenter() {
      wasserstein_ = 2;
      alpha_ = 1;
      lambda_ = 1;
      numberOfInputs_ = 0;
      threadNumber_ = 1;
      time_limit_ = 1;
      deterministic_ = 1;
      reinit_prices_ = 1;
      epsilon_decreases_ = 1;
      use_progressive_ = 1;
      this->setDebugMsgPrefix("PersistenceDiagramBarycenter");
    }

    ~PersistenceDiagramBarycenter() = default;

    void execute(
      std::vector<std::vector<diagramTuple>> &intermediateDiagrams,
      std::vector<diagramTuple> &barycenter,
      std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings);

    // 		inline int setDiagram(int idx, void* data){
    // 			if(idx < numberOfInputs_){
    // 			inputData_[idx] = data;
    // 			}
    // 			else{
    // 			return -1;
    // 			}
    // 			return 0;
    // 		}

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

    inline void setDeterministic(const bool deterministic) {
      deterministic_ = deterministic;
    }

    inline void setWasserstein(const std::string &wasserstein) {
      wasserstein_ = (wasserstein == "inf") ? -1 : stoi(wasserstein);
    }

    inline void setUseProgressive(const bool use_progressive) {
      if(use_progressive)
        epsilon_decreases_ = true;
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

    template <typename type>
    static type abs(const type var) {
      return (var >= 0) ? var : -var;
    }

    inline void setMethod(const int &method) {
      method_ = method;
    }

    inline void setReinitPrices(const bool reinit_prices) {
      reinit_prices_ = reinit_prices;
    }

    inline void setEpsilonDecreases(const bool epsilon_decreases) {
      if(use_progressive_)
        epsilon_decreases_ = true;
      else
        epsilon_decreases_ = epsilon_decreases;
    }

    inline void setEarlyStoppage(const bool early_stoppage) {
      early_stoppage_ = early_stoppage;
    }

  protected:
    bool deterministic_;
    int method_;
    int wasserstein_;
    int numberOfInputs_;
    bool use_progressive_;
    double alpha_;
    double lambda_;
    double time_limit_;

    int points_added_;
    int points_deleted_;

    std::vector<std::vector<dataType>> all_matchings_;
    std::vector<std::vector<dataType>> all_old_matchings_;
    std::vector<BidderDiagram<dataType>> bidder_diagrams_;
    std::vector<GoodDiagram<dataType>> barycenter_goods_;

    bool reinit_prices_;
    bool epsilon_decreases_;
    bool early_stoppage_;
  };

  template <typename dataType>
  void PersistenceDiagramBarycenter<dataType>::execute(
    std::vector<std::vector<diagramTuple>> &intermediateDiagrams,
    std::vector<diagramTuple> &barycenter,
    std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings) {

    Timer tm;
    {
      printMsg("Computing Barycenter of " + std::to_string(numberOfInputs_)
               + " diagrams.");

      std::vector<std::vector<diagramTuple>> data_min(numberOfInputs_);
      std::vector<std::vector<diagramTuple>> data_sad(numberOfInputs_);
      std::vector<std::vector<diagramTuple>> data_max(numberOfInputs_);

      std::vector<std::vector<int>> data_min_idx(numberOfInputs_);
      std::vector<std::vector<int>> data_sad_idx(numberOfInputs_);
      std::vector<std::vector<int>> data_max_idx(numberOfInputs_);

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

      std::vector<diagramTuple> barycenter_min;
      std::vector<diagramTuple> barycenter_sad;
      std::vector<diagramTuple> barycenter_max;

      std::vector<std::vector<matchingTuple>> matching_min, matching_sad,
        matching_max;

      dataType total_cost = 0;
      if(do_min && do_max) {
        time_limit_ = time_limit_ / 2;
      }
      if(do_sad) {
        time_limit_ = time_limit_ / 3;
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
      if(do_min) {
        printMsg("Computing Minima barycenter...");
        PDBarycenter<dataType> bary_min = PDBarycenter<dataType>();
        bary_min.setThreadNumber(threadNumber_);
        bary_min.setWasserstein(wasserstein_);
        bary_min.setNumberOfInputs(numberOfInputs_);
        bary_min.setDiagramType(0);
        bary_min.setUseProgressive(use_progressive_);
        bary_min.setTimeLimit(time_limit_);
        bary_min.setGeometricalFactor(alpha_);
        bary_min.setDebugLevel(debugLevel_);
        bary_min.setDeterministic(deterministic_);
        bary_min.setLambda(lambda_);
        bary_min.setMethod(method_);
        bary_min.setEarlyStoppage(early_stoppage_);
        bary_min.setEpsilonDecreases(epsilon_decreases_);
        bary_min.setReinitPrices(reinit_prices_);
        bary_min.setDiagrams(&data_min);
        matching_min = bary_min.execute(barycenter_min);
        total_cost += bary_min.getCost();
      }
      /*}

      #ifdef TTK_ENABLE_OPENMP
      #pragma omp section
      #endif
      {*/
      if(do_sad) {
        printMsg("Computing Saddles barycenter...");
        PDBarycenter<dataType> bary_sad = PDBarycenter<dataType>();
        bary_sad.setThreadNumber(threadNumber_);
        bary_sad.setWasserstein(wasserstein_);
        bary_sad.setNumberOfInputs(numberOfInputs_);
        bary_sad.setDiagramType(1);
        bary_sad.setUseProgressive(use_progressive_);
        bary_sad.setTimeLimit(time_limit_);
        bary_sad.setGeometricalFactor(alpha_);
        bary_sad.setLambda(lambda_);
        bary_sad.setDebugLevel(debugLevel_);
        bary_sad.setMethod(method_);
        bary_sad.setEarlyStoppage(early_stoppage_);
        bary_sad.setEpsilonDecreases(epsilon_decreases_);
        bary_sad.setDeterministic(deterministic_);
        bary_sad.setReinitPrices(reinit_prices_);
        bary_sad.setDiagrams(&data_sad);
        matching_sad = bary_sad.execute(barycenter_sad);
        total_cost += bary_sad.getCost();
      }
      /*}

      #ifdef TTK_ENABLE_OPENMP
      #pragma omp section
      #endif
      {*/
      if(do_max) {
        printMsg("Computing Maxima barycenter...");
        PDBarycenter<dataType> bary_max = PDBarycenter<dataType>();
        bary_max.setThreadNumber(threadNumber_);
        bary_max.setWasserstein(wasserstein_);
        bary_max.setNumberOfInputs(numberOfInputs_);
        bary_max.setDiagramType(2);
        bary_max.setUseProgressive(use_progressive_);
        bary_max.setTimeLimit(time_limit_);
        bary_max.setGeometricalFactor(alpha_);
        bary_max.setLambda(lambda_);
        bary_max.setMethod(method_);
        bary_max.setDebugLevel(debugLevel_);
        bary_max.setEarlyStoppage(early_stoppage_);
        bary_max.setDeterministic(deterministic_);
        bary_max.setEpsilonDecreases(epsilon_decreases_);
        bary_max.setReinitPrices(reinit_prices_);
        bary_max.setDiagrams(&data_max);
        matching_max = bary_max.execute(barycenter_max);
        total_cost += bary_max.getCost();
      }
      //}
      //}

      // Reconstruct matchings
      all_matchings.resize(1);
      all_matchings[0].resize(numberOfInputs_);
      for(int i = 0; i < numberOfInputs_; i++) {

        if(do_min) {
          for(unsigned int j = 0; j < matching_min[i].size(); j++) {
            matchingTuple t = matching_min[i][j];
            int bidder_id = std::get<0>(t);
            std::get<0>(t) = data_min_idx[i][bidder_id];
            if(std::get<1>(t) < 0) {
              std::get<1>(t) = -1;
            }
            all_matchings[0][i].push_back(t);
          }
        }

        if(do_sad) {
          for(unsigned int j = 0; j < matching_sad[i].size(); j++) {
            matchingTuple t = matching_sad[i][j];
            int bidder_id = std::get<0>(t);
            std::get<0>(t) = data_sad_idx[i][bidder_id];
            if(std::get<1>(t) >= 0) {
              std::get<1>(t) = std::get<1>(t) + barycenter_min.size();
            } else {
              std::get<1>(t) = -1;
            }
            all_matchings[0][i].push_back(t);
          }
        }

        if(do_max) {
          for(unsigned int j = 0; j < matching_max[i].size(); j++) {
            matchingTuple t = matching_max[i][j];
            int bidder_id = std::get<0>(t);
            std::get<0>(t) = data_max_idx[i][bidder_id];
            if(std::get<1>(t) >= 0) {
              std::get<1>(t) = std::get<1>(t) + barycenter_min.size()
                               + barycenter_sad.size();
            } else {
              std::get<1>(t) = -1;
            }
            all_matchings[0][i].push_back(t);
          }
        }
      }
      // Reconstruct barcenter
      for(unsigned int j = 0; j < barycenter_min.size(); j++) {
        diagramTuple dt = barycenter_min[j];
        barycenter.push_back(dt);
      }
      for(unsigned int j = 0; j < barycenter_sad.size(); j++) {
        diagramTuple dt = barycenter_sad[j];
        barycenter.push_back(dt);
      }
      for(unsigned int j = 0; j < barycenter_max.size(); j++) {
        diagramTuple dt = barycenter_max[j];
        barycenter.push_back(dt);
      }

      // Recreate 3D critical coordinates of barycentric points
      std::vector<int> number_of_matchings_for_point(barycenter.size());
      std::vector<float> cords_x1(barycenter.size());
      std::vector<float> cords_y1(barycenter.size());
      std::vector<float> cords_z1(barycenter.size());
      std::vector<float> cords_x2(barycenter.size());
      std::vector<float> cords_y2(barycenter.size());
      std::vector<float> cords_z2(barycenter.size());
      for(unsigned i = 0; i < barycenter.size(); i++) {
        number_of_matchings_for_point[i] = 0;
        cords_x1[i] = 0;
        cords_y1[i] = 0;
        cords_z1[i] = 0;
        cords_x2[i] = 0;
        cords_y2[i] = 0;
        cords_z2[i] = 0;
      }

      for(unsigned i = 0; i < all_matchings[0].size(); i++) {
        std::vector<diagramTuple> &CTDiagram = intermediateDiagrams[i];
        for(unsigned j = 0; j < all_matchings[0][i].size(); j++) {
          matchingTuple t = all_matchings[0][i][j];
          int bidder_id = std::get<0>(t);
          int bary_id = std::get<1>(t);

          diagramTuple &bidder = CTDiagram[bidder_id];
          number_of_matchings_for_point[bary_id] += 1;
          cords_x1[bary_id] += std::get<7>(bidder);
          cords_y1[bary_id] += std::get<8>(bidder);
          cords_z1[bary_id] += std::get<9>(bidder);
          cords_x2[bary_id] += std::get<11>(bidder);
          cords_y2[bary_id] += std::get<12>(bidder);
          cords_z2[bary_id] += std::get<13>(bidder);
        }
      }

      for(unsigned i = 0; i < barycenter.size(); i++) {
        if(number_of_matchings_for_point[i] > 0) {
          std::get<7>(barycenter[i])
            = cords_x1[i] / number_of_matchings_for_point[i];
          std::get<8>(barycenter[i])
            = cords_y1[i] / number_of_matchings_for_point[i];
          std::get<9>(barycenter[i])
            = cords_z1[i] / number_of_matchings_for_point[i];
          std::get<11>(barycenter[i])
            = cords_x2[i] / number_of_matchings_for_point[i];
          std::get<12>(barycenter[i])
            = cords_y2[i] / number_of_matchings_for_point[i];
          std::get<13>(barycenter[i])
            = cords_z2[i] / number_of_matchings_for_point[i];
        }
      }

      printMsg("Total cost : " + std::to_string(total_cost));
      // std::stringstream msg;
      // msg << "[PersistenceDiagramBarycenter] processed in "
      //     << tm.getElapsedTime() << " s. (" << threadNumber_ << "
      //     thread(s))."
      //     << std::endl;
      // dMsg(std::cout, msg.str(), timeMsg);
      printMsg("Complete", 1, tm.getElapsedTime(), threadNumber_);
    }
  }

} // namespace ttk
