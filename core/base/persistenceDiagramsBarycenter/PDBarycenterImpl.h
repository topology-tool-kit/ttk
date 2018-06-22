#ifndef _PDBARYCENTERIMPL_H
#define _PDBARYCENTERIMPL_H

#include <stdlib.h>     /* srand, rand */
#include <cmath>

using namespace ttk;

template <typename dataType>
int PersistenceDiagramsBarycenter<dataType>::execute(){

	Timer t;
	{
	this->setBidderDiagrams();
	this->setInitialBarycenter();
	
	std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> pair = this->getKDTree();
	KDTree<dataType>* kdt = pair.first;
	std::vector<KDTree<dataType>*>& correspondance_kdt_map = pair.second;
	
	std::vector<std::vector<matchingTuple>> all_matchings;
	dataType total_costs = 0;
	for(int i=0; i<numberOfInputs_; i++){
		Auction<dataType> auction = Auction<dataType>(bidder_diagrams_[i], barycenter_goods_[i], wasserstein_, geometrical_factor_, 0.01, kdt, correspondance_kdt_map);
		std::cout<< "Barycenter size : "<< barycenter_goods_[i].size() << std::endl;
		int n_biddings = 0;
		auction.buildUnassignedBidders();
		auction.reinitializeGoods();
		auction.runAuctionRound(n_biddings, i);
		std::vector<matchingTuple> matchings;
		dataType cost = auction.getMatchingsAndDistance(&matchings, true);
		all_matchings.push_back(matchings);
		total_costs += cost;
	}
	

	dataType max_shift = 0;
	dataType average_shift =0;
	updateBarycenter(all_matchings, max_shift, average_shift);
	std::cout<< "Barycenter size : "<< barycenter_goods_[0].size() << std::endl;
	delete kdt;
	
	
	std::stringstream msg;
	msg << "[PersistenceDiagramsBarycenter] processed in "
		<< t.getElapsedTime() << " s. (" << threadNumber_
		<< " thread(s))."
		<< std::endl;
	dMsg(std::cout, msg.str(), timeMsg);
	}

	return 0;
}







template <typename dataType>
void PersistenceDiagramsBarycenter<dataType>::updateBarycenter(std::vector<std::vector<matchingTuple>>& matchings, dataType& max_shift, dataType& average_shift){
	
	// 1. Initialize variables used in the sequel
	unsigned int n_goods = barycenter_goods_[0].size();
	unsigned int n_diagrams = bidder_diagrams_.size();
	
	std::vector<int> count_diag_matchings(n_goods);     // Number of diagonal matchings for each point of the barycenter
	std::vector<int> x(n_goods);
	std::vector<int> y(n_goods);
	for(unsigned int i=0; i<n_goods; i++){
		count_diag_matchings[i] = 0;
		x[i] = 0;
		y[i] = 0;
	}
	std::vector<dataType> min_prices(n_diagrams);
	for(unsigned int j=0; j<n_diagrams; j++){
		min_prices[j] = std::numeric_limits<dataType>::max();
	}
	
	std::vector<Bidder<dataType>*> points_to_append;  //Will collect bidders linked to diagonal
	
	// 2. Preprocess the matchings
	for(unsigned int j=0; j<matchings.size(); j++){
		for(unsigned int i=0; i<matchings[j].size(); i++){
			int bidder_id = std::get<0>(matchings[j][i]);
			int good_id = std::get<1>(matchings[j][i]);
			if(good_id<0 && bidder_id>=0){
				// Future new barycenter points
				points_to_append.push_back(&bidder_diagrams_[j].get(bidder_id));
			}
			else if(good_id>=0 && bidder_id>=0){
				// Update coordinates (to be divided by the number of diagrams later on)
				x[good_id] += bidder_diagrams_[j].get(bidder_id).x_;
				y[good_id] += bidder_diagrams_[j].get(bidder_id).y_;
			}
			else if(good_id>=0 && bidder_id<0){
				// Counting the number of times this point is linked to the diagonal
				count_diag_matchings[good_id] = count_diag_matchings[good_id] + 1;
			}
		}
	}
	
	// 3. Update the previous points of the barycenter
	for(unsigned int i=0; i<n_goods; i++){
		if(count_diag_matchings[i]<n_diagrams){
			// Barycenter point i is matched at least to one off-diagonal bidder
			// 3.1 Compute the arithmetic mean of off-diagonal bidders linked to it
			dataType x_bar = x[i]/(n_diagrams - count_diag_matchings[i]);
			dataType y_bar = y[i]/(n_diagrams - count_diag_matchings[i]);
			// 3.2 Compute the new coordinates of the point (the more linked to the diagonal it was, the closer to the diagonal it'll be)
			dataType new_x = ( (n_diagrams - count_diag_matchings[i])*x_bar  + count_diag_matchings[i]*(x_bar+y_bar)/2 )/n_diagrams;
			dataType new_y = ( (n_diagrams - count_diag_matchings[i])*y_bar  + count_diag_matchings[i]*(x_bar+y_bar)/2 )/n_diagrams;
			// 3.3 Compute and store how much the point has shifted
			dataType dx = barycenter_goods_[0].get(i).x_ - new_x;
			dataType dy = barycenter_goods_[0].get(i).y_ - new_y;
			dataType shift = pow(abs(dx), wasserstein_) + pow(abs(dy), wasserstein_);
			average_shift += shift/n_goods;
			if(shift>max_shift){
				max_shift = shift;
			}
			// 3.4 Update the position of the point
			for(unsigned int j=0; j<n_diagrams; j++){
				barycenter_goods_[j].get(i).SetCoordinates(new_x, new_y);
				if(barycenter_goods_[j].get(i).getPrice()<min_prices[j]){
					min_prices[j] = barycenter_goods_[j].get(i).getPrice();
				}
			}
			// TODO Reintitialize/play with prices here if you wish
		}
	}
	
	// 4. Delete off-diagonal barycenter points not linked to any
	// off-diagonal bidder
	for(unsigned int i=0; i<n_goods; i++){
		if(count_diag_matchings[i] == n_diagrams){
			dataType shift = barycenter_goods_[0].get(i).getPersistence() / pow(2, 1./wasserstein_);
			average_shift += shift/n_goods;
			if(shift>max_shift){
				max_shift = shift;
			}
			for(unsigned int j=0; j<n_diagrams; j++){
				barycenter_goods_[j].get(i).id_ = -1;
			}
		}
	}
	
	// 5. Append the new points to the barycenter
	for(unsigned int k=0; k<points_to_append.size(); k++){
		Bidder<dataType>* b = points_to_append[k];
		dataType x = b->x_ + (n_diagrams -1)*(b->x_+b->y_)/(2 * n_diagrams); 
		dataType y = b->y_ + (n_diagrams -1)*(b->x_+b->y_)/(2 * n_diagrams); 
		for(unsigned int j=0; j<n_diagrams; j++){
			Good<dataType> g = Good<dataType>(x, y, false, barycenter_goods_[j].size());
			g.setPrice(min_prices[j]);
			barycenter_goods_[j].addGood(g);
		}
	}
	
	// 6. Finally, recreate barycenter_goods via copy :/
	// TODO avoid all these copies of goods
	for(unsigned int j=0; j<n_diagrams; j++){
		int count = 0;
		GoodDiagram<dataType> new_barycenter;
		for(unsigned int i=0; i<barycenter_goods_[j].size(); i++){
			std::cout<< "Barycenter " << j<< " has size : " <<barycenter_goods_[j].size() <<std::endl;
			Good<dataType>& g = barycenter_goods_[j].get(i);
			if(g.id_!=-1){
				g.id_ = count;
				new_barycenter.addGood(g);
				count ++;
			}
		}
		barycenter_goods_[j] = new_barycenter;
	}
	return;
}



template <typename dataType>
dataType PersistenceDiagramsBarycenter<dataType>::getEpsilon(dataType rho){
	return pow(rho, 2)/8;
}

template <typename dataType>
dataType PersistenceDiagramsBarycenter<dataType>::getRho(dataType epsilon){
	return std::sqrt(8*epsilon);
}



template <typename dataType>
void PersistenceDiagramsBarycenter<dataType>::setBidderDiagrams(){
	for(int i=0; i<numberOfInputs_; i++){
		std::vector<diagramTuple>* CTDiagram = static_cast<std::vector<diagramTuple>*>(inputData_[i]);
		BidderDiagram<dataType> bidders;
		for(unsigned int j=0; j<CTDiagram->size(); j++){
			//Add bidder to bidders
			Bidder<dataType> b = Bidder<dataType>((*CTDiagram)[j], j);
			b.setPositionInAuction(bidders.size());
			bidders.addBidder(b);
		}
		bidder_diagrams_.push_back(bidders);
	}
	return;
}


template <typename dataType>
void PersistenceDiagramsBarycenter<dataType>::setInitialBarycenter(){
	//int random_idx = rand() % numberOfInputs_;
	std::cout << "BEWARE, initial barycenter is not chosen randomly..."<< std::endl;
	int random_idx = 0;
	std::vector<diagramTuple>* CTDiagram = static_cast<std::vector<diagramTuple>*>(inputData_[random_idx]);
	
	for(int i=0; i<numberOfInputs_; i++){
		GoodDiagram<dataType> goods;
		for(unsigned int j=0; j<CTDiagram->size(); j++){
			//Add bidder to bidders
			Good<dataType> g = Good<dataType>((*CTDiagram)[j], j);
			goods.addGood(g);
		}
		barycenter_goods_.push_back(goods);
	}
}


template <typename dataType>
std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> PersistenceDiagramsBarycenter<dataType>::getKDTree(){
	Timer t;
	KDTree<dataType>* kdt = new KDTree<dataType>(true, wasserstein_);
	
	const int dimension = geometrical_factor_ >= 1 ? 2 : 5;
	
	std::vector<dataType> coordinates;
	std::vector<std::vector<dataType>> weights;
	
	for(int i=0; i<barycenter_goods_[0].size(); i++){
		Good<dataType>& g = barycenter_goods_[0].get(i);
		coordinates.push_back(geometrical_factor_*g.x_);
		coordinates.push_back(geometrical_factor_*g.y_);
		if(geometrical_factor_<1){
			coordinates.push_back((1-geometrical_factor_)*g.coords_x_);
			coordinates.push_back((1-geometrical_factor_)*g.coords_y_);
			coordinates.push_back((1-geometrical_factor_)*g.coords_z_);
		}
	}
	
	for(unsigned int idx=0; idx<barycenter_goods_.size(); idx++){
		std::vector<dataType> empty_weights;
		weights.push_back(empty_weights);
		for(int i=0; i<barycenter_goods_[idx].size(); i++){
			Good<dataType>& g = barycenter_goods_[idx].get(i);
			weights[idx].push_back(g.getPrice());
		}
	}
	std::vector<KDTree<dataType>*> correspondance_kdt_map = kdt->build(coordinates.data(), barycenter_goods_[0].size(), dimension, weights, barycenter_goods_.size());
	std::cout<<"[Building KD-Tree] Time elapsed : " << t.getElapsedTime() << " s."<<std::endl;
	return std::make_pair(kdt, correspondance_kdt_map);
}







#endif
