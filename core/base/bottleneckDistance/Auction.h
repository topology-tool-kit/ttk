/// \ingroup base
/// \class ttk::Auction
/// \author Joseph Budin <joseph.budin@polytechnique.edu>

#ifndef _AUCTION_H
#define _AUCTION_H

#ifndef matchingTuple
#define matchingTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::idVertex, dataType>
#endif

#ifndef diagramTuple
#define diagramTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::NodeType, ttk::ftm::idVertex, \
  ttk::ftm::NodeType, dataType, ttk::ftm::idVertex, \
  dataType, float, float, float, dataType, float, float, float>
#endif


#include <cmath>
#include <limits>
#include <iostream>
#include <Debug.h>
#include <PersistenceDiagram.h>
#include <AuctionActor.h>
#include <KDTree.h>
#include <queue>

namespace ttk {
  template<typename dataType>
  struct Compare {
	constexpr bool operator()(std::pair<int, dataType> const & a,
								std::pair<int, dataType> const & b) const noexcept
    { return a.second > b.second;}
  };
	
	
    
  template<typename dataType>
  class Auction : public Debug
  {

    public:
		KDTree<dataType>* kdt_;
		std::vector<KDTree<dataType>*> correspondance_kdt_map_;
		
		Auction(int wasserstein, double geometricalFactor, double delta_lim, bool use_kdTree=true) {
            n_bidders_ = 0;
            n_goods_ = 0;
			epsilon_ = 1;
			wasserstein_ = wasserstein;
			delta_lim_ = delta_lim;
			geometricalFactor_ = geometricalFactor;
			use_kdt_ = use_kdTree;
        };
		
		
		Auction(BidderDiagram<dataType>& bidders, GoodDiagram<dataType>& goods, int wasserstein, double geometricalFactor, double delta_lim, KDTree<dataType>* kdt, std::vector<KDTree<dataType>*>& correspondance_kdt_map, dataType epsilon, bool use_kdTree=true) {
			bidders_ = bidders;
			goods_ = goods;
			
            n_bidders_ = bidders.size();
            n_goods_ = goods.size();
			
			for(int i=0; i < n_bidders_; i++){
				//Add diagonal goods
				Bidder<dataType>& b = bidders_.get(i);
				Good<dataType> g = Good<dataType>(b.x_, b.y_, true, -b.id_-1);
				g.projectOnDiagonal();
				diagonal_goods_.addGood(g);
				std::pair<int, dataType> pair = std::make_pair(i, g.getPrice());
				diagonal_queue_.push(pair);
			}
			for(int i=0; i < n_goods_; i++){
				//Add diagonal bidders
				Good<dataType>& g = goods_.get(i);
				Bidder<dataType> b = Bidder<dataType>(g.x_, g.y_, true, -g.id_-1);
				b.projectOnDiagonal();
				b.setPositionInAuction(bidders_.size());
				bidders_.addBidder(b);
			}
			
			epsilon_ = epsilon;
			wasserstein_ = wasserstein;
			delta_lim_ = delta_lim;
			geometricalFactor_ = geometricalFactor;
			
			use_kdt_ = use_kdTree;
			if(use_kdt_) {
				kdt_ = kdt;
				correspondance_kdt_map_ = correspondance_kdt_map;
			}
        };
		
		
		~Auction() {};

		void runAuctionRound(int& n_biddings, const int kdt_index=0);
		dataType getMatchingsAndDistance(std::vector<matchingTuple> *matchings, bool get_diagonal_matches=false);
		dataType run(std::vector<matchingTuple> *matchings);

		
		void BuildAuctionDiagrams(std::vector<diagramTuple> diagram1, std::vector<diagramTuple> diagram2){
			Timer t;
			n_bidders_ = diagram1.size();
			n_goods_ = diagram2.size();
			this->setBidders(diagram1);
			this->setGoods(diagram2);
			for(int i=0; i < n_bidders_; i++){
				//Add diagonal goods
				Bidder<dataType>& b = bidders_.get(i);
				Good<dataType> g = Good<dataType>(b.x_, b.y_, true, -b.id_-1);
				g.projectOnDiagonal();
				diagonal_goods_.addGood(g);
				std::pair<int, dataType> pair = std::make_pair(i, g.getPrice());
				diagonal_queue_.push(pair);
			}
			for(int i=0; i < n_goods_; i++){
				//Add diagonal bidders
				Good<dataType>& g = goods_.get(i);
				Bidder<dataType> b = Bidder<dataType>(g.x_, g.y_, true, -g.id_-1);
				b.projectOnDiagonal();
				b.setPositionInAuction(bidders_.size());
				bidders_.addBidder(b);
			}
			this->buildKDTree();
			std::cout<<"[Initialize auction] Time elapsed : " << t.getElapsedTime() << " s."<<std::endl;
		}
		
		void setBidders(std::vector<diagramTuple> diagram1){
			int d1Size = (int) diagram1.size();
			
			for(int i=0; i<d1Size; i++){
				//Add bidder to bidders
				Bidder<dataType> b = Bidder<dataType>(diagram1[i], i);
				b.setPositionInAuction(bidders_.size());
				bidders_.addBidder(b);
			}
			n_bidders_ = bidders_.size();
		}

		void setGoods(std::vector<diagramTuple> diagram2){
			int d2Size = (int) diagram2.size();
			
			for(int i=0; i<d2Size; i++){
				//Add bidder to bidders
				Good<dataType> g = Good<dataType>(diagram2[i], i);
				goods_.addGood(g);
			}	
			n_goods_ = goods_.size();
		}
		
		
		void buildKDTree(){
			Timer t;
			kdt_ = new KDTree<dataType>(true, wasserstein_);
			const int dimension = geometricalFactor_ >= 1 ? 2 : 5;
			std::vector<dataType> coordinates;
			for(int i=0; i<goods_.size(); i++){
				Good<dataType>& g = goods_.get(i);
				coordinates.push_back(geometricalFactor_*g.x_);
				coordinates.push_back(geometricalFactor_*g.y_);
				if(geometricalFactor_<1){
					coordinates.push_back((1-geometricalFactor_)*g.coords_x_);
					coordinates.push_back((1-geometricalFactor_)*g.coords_y_);
					coordinates.push_back((1-geometricalFactor_)*g.coords_z_);
				}
			}
			correspondance_kdt_map_ = kdt_->build(coordinates.data(), goods_.size(), dimension);
			std::cout<<"[Building KD-Tree] Time elapsed : " << t.getElapsedTime() << " s."<<std::endl;
		}
		
		void setEpsilon(dataType epsilon){
			epsilon_ = epsilon;
		}
		
		void initializeEpsilon(){
			dataType max_persistence = 0;
			for(int i=0; i<bidders_.size(); i++){
				Bidder<dataType>& b = bidders_.get(i);
				dataType persistence = b.getPersistence();
				if(persistence>max_persistence){
					max_persistence = persistence;
				}
			}
			
			for(int i=0; i<goods_.size(); i++){
				Good<dataType>& g = goods_.get(i);
				dataType persistence = g.getPersistence();
				if(persistence>max_persistence){
					max_persistence = persistence;
				}
			}
			this->epsilon_ = 5/4 * pow(max_persistence, wasserstein_);
		}
		
		
		void buildUnassignedBidders(){
			for(int i=0; i<bidders_.size(); i++){
				Bidder<dataType>& b = bidders_.get(i);
				b.setProperty(NULL);
				unassignedBidders_.push_back(i);
			}
		}
		
		
		void reinitializeGoods(){
			for(int i=0; i<goods_.size(); i++){
				Good<dataType>& g = goods_.get(i);
				g.setOwner(-1);
			}
			for(int i=0; i<diagonal_goods_.size(); i++){
				Good<dataType>& g = diagonal_goods_.get(i);
				g.setOwner(-1);
			}
		}
		
		
		dataType getMatchingDistance(){
			dataType d = 0;
			for(int i=0; i<bidders_.size(); i++){
				Bidder<dataType>& b = bidders_.get(i);
				d += b.cost(b.getProperty(), wasserstein_, geometricalFactor_); 
			}
			return d;
		}
		
		
		dataType getRelativePrecision(){
			dataType d = this->getMatchingDistance();
			if(d<1e-12){
				return 0;
			}
			dataType denominator = d - bidders_.size()*epsilon_;
			if(denominator<=0){
				return 1;
			}
			else{
				return pow(d/denominator, 1/((float)wasserstein_)) - 1;
			}
		}
		
		
		template<typename type>
		static type abs(const type var) {
			return (var >= 0) ? var : -var;
		}
	
    protected:
		int wasserstein_;   // Power in Wassertsein distance (by default set to 2)
		BidderDiagram<dataType>  bidders_;
		GoodDiagram<dataType>  goods_;
		GoodDiagram<dataType> diagonal_goods_;
		std::priority_queue<std::pair<int, dataType>, std::vector<std::pair<int, dataType>>, Compare<dataType>> diagonal_queue_;
		std::list<int> unassignedBidders_;
		
		int n_bidders_;
		int n_goods_;
		
		dataType epsilon_;
		double delta_lim_;
		double geometricalFactor_;
		bool use_kdt_;
		
		//KDTree<dataType>* kdt_;
  };
}

#include <AuctionImpl.h>

#endif

