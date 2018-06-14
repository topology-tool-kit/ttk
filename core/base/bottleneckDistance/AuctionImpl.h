#ifndef _AUCTIONIMPL_H
#define _AUCTIONIMPL_H

#ifndef matchingTuple
#define matchingTuple std::tuple<ftm::idVertex, ftm::idVertex, dataType>
#endif

template <typename dataType>
dataType ttk::Auction<dataType>::run(std::vector<matchingTuple> *matchings)
{
	dataType delta = 5;
	while(delta>delta_lim_){
		epsilon_ /= 5;
		this->buildUnassignedBidders();
		this->reinitializeGoods();
		while(unassignedBidders_.size()>0){
			int idx = unassignedBidders_.front();
			Bidder<dataType>& b = bidders_.get(idx);
			unassignedBidders_.pop_front();
			
			GoodDiagram<dataType>& all_goods = b.isDiagonal() ? diagonal_goods_ : goods_;
			Good<dataType>& diagonal_good = b.id_>=0 ? diagonal_goods_.get(b.id_) : goods_.get(-b.id_-1);
			
			int idx_reassigned;
			if(b.isDiagonal()){
				idx_reassigned = b.runBidding(all_goods, diagonal_good, wasserstein_, epsilon_, geometricalFactor_, correspondance_kdt_map_);
			}
			else{
				// We can use the kd-tree to speed up the search
				idx_reassigned = b.runBidding(all_goods, diagonal_good, wasserstein_, epsilon_, geometricalFactor_, kdt_);
			}
			
			if(idx_reassigned>=0){
				Bidder<dataType>& reassigned = bidders_.get(idx_reassigned);
				reassigned.setProperty(NULL);
				unassignedBidders_.push_back(idx_reassigned);
			}
		}
		delta = this->getRelativePrecision();
	}
	
	dataType wassersteinDistance = 0;
	for (int i=0; i<bidders_.size(); i++){
		Bidder<dataType> b = bidders_.get(i);
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
#endif
