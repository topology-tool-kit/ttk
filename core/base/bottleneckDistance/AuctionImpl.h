#ifndef _AUCTIONIMPL_H
#define _AUCTIONIMPL_H

#ifndef matchingTuple
#define matchingTuple std::tuple<ftm::idVertex, ftm::idVertex, dataType>
#endif


template <typename dataType>
void ttk::Auction<dataType>::runAuctionRound(int& n_biddings)
{
	while(unassignedBidders_.size()>0){
		n_biddings ++;
		int idx = unassignedBidders_.front();
		Bidder<dataType>& b = bidders_.get(idx);
		unassignedBidders_.pop_front();
		
		GoodDiagram<dataType>& all_goods = b.isDiagonal() ? diagonal_goods_ : goods_;
		Good<dataType>& twin_good = b.id_>=0 ? diagonal_goods_.get(b.id_) : goods_.get(-b.id_-1);
		
		int idx_reassigned;
		if(b.isDiagonal()){
			if(use_kdt_){
				idx_reassigned = b.runDiagonalKDTBidding(all_goods, twin_good, wasserstein_, epsilon_, geometricalFactor_, correspondance_kdt_map_, diagonal_queue_);
			}
			else{
				idx_reassigned = b.runDiagonalBidding(all_goods, twin_good, wasserstein_, epsilon_, geometricalFactor_, diagonal_queue_);
			}
		}
		else{
			if(use_kdt_){
				// We can use the kd-tree to speed up the search
				idx_reassigned = b.runKDTBidding(all_goods, twin_good, wasserstein_, epsilon_, geometricalFactor_, kdt_);
			}
			else{
				idx_reassigned = b.runBidding(all_goods, twin_good, wasserstein_, epsilon_, geometricalFactor_);
			}
		}
		
		if(idx_reassigned>=0){
			Bidder<dataType>& reassigned = bidders_.get(idx_reassigned);
			reassigned.setProperty(NULL);
			unassignedBidders_.push_back(idx_reassigned);
		}
	}
}


template <typename dataType>
dataType ttk::Auction<dataType>::getMatchingsAndDistance(std::vector<matchingTuple> *matchings)
{
	dataType wassersteinDistance = 0;
	for (int i=0; i<bidders_.size(); i++){
		Bidder<dataType>& b = bidders_.get(i);
		if(!b.isDiagonal()){
			int good_id = b.getProperty()->id_;
			dataType cost;
			
			if(good_id>-1){
				// good is not diagonal
				cost = b.cost(b.getProperty(), wasserstein_, geometricalFactor_);
				matchingTuple t = std::make_tuple(i, good_id, cost);
				matchings->push_back(t);
			}
			else{
				cost = pow(abs<dataType>((b.x_-b.y_)/2), wasserstein_);
			}
			wassersteinDistance += cost;
		}
		else{
			// b is diagonal
			Good<dataType> g = *b.getProperty();
			wassersteinDistance += pow(abs<dataType>((g.x_-g.y_)/2), wasserstein_);
		}
	}
	return wassersteinDistance;
}


template <typename dataType>
dataType ttk::Auction<dataType>::run(std::vector<matchingTuple> *matchings)
{	
	int n_biddings = 0;
	dataType delta = 5;
	while(delta>delta_lim_){
		epsilon_ /= 5;
		this->buildUnassignedBidders();
		this->reinitializeGoods();
		this->runAuctionRound(n_biddings);
		delta = this->getRelativePrecision();
	}
	std::cout<<"Number of biddings : " << n_biddings <<std::endl;
	dataType wassersteinDistance = this->getMatchingsAndDistance(matchings);
	return wassersteinDistance;
}
#endif
