#ifndef _PDBARYCENTERIMPL_H
#define _PDBARYCENTERIMPL_H

#define BLocalMax ttk::ftm::NodeType::Local_maximum
#define BLocalMin ttk::ftm::NodeType::Local_minimum
#define BSaddle1  ttk::ftm::NodeType::Saddle1
#define BSaddle2  ttk::ftm::NodeType::Saddle2

#include <stdlib.h>     /* srand, rand */
#include <cmath>

using namespace ttk;

template <typename dataType>
std::vector<std::vector<matchingTuple>> PDBarycenter<dataType>::execute(std::vector<diagramTuple>& barycenter){

	
	std::vector<std::vector<matchingTuple>> previous_matchings;
		
	this->setBidderDiagrams();
	
	dataType max_persistence = getMaxPersistence();
	dataType epsilon_0 = getEpsilon(max_persistence);
	dataType epsilon = epsilon_0;
	dataType min_persistence = max_persistence/2.;
	
	this->enrichCurrentBidderDiagrams(2*max_persistence, min_persistence);
	this->setInitialBarycenter(min_persistence);
	std::cout<< "Barycenter size : "<< barycenter_goods_[0].size() << std::endl;
	
	
	dataType previous_cost = std::numeric_limits<dataType>::max();
	std::vector<dataType> min_diag_price(numberOfInputs_);
	for(int i=0; i<numberOfInputs_; i++){
		min_diag_price[i] = 0;
	}		
	int n_iterations = 0;
	
	bool converged = false;
	while(!converged){
		n_iterations += 1;
		dataType rho = getRho(epsilon);
		if(n_iterations>1 && min_persistence>rho){
			this->enrichCurrentBidderDiagrams(min_persistence, rho);
			epsilon = getEpsilon(min_persistence);
			min_persistence = rho;
			// TODO Enrich barycenter using median diagonal and off-diagonal prices
		}
		std::cout<< "epsilon : "<< epsilon << std::endl;
		std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> pair = this->getKDTree();
		KDTree<dataType>* kdt = pair.first;
		std::vector<KDTree<dataType>*>& correspondance_kdt_map = pair.second;
		
		std::vector<std::vector<matchingTuple>> all_matchings(numberOfInputs_);
		std::vector<int> sizes(numberOfInputs_);
		for(int i=0; i<numberOfInputs_; i++){
			sizes[i] = current_bidder_diagrams_[i].size();
		}
		
		
		dataType total_cost = 0;
		
		/*#ifdef TTK_ENABLE_OPENMP
		omp_set_num_threads(threadNumber_);
		#pragma omp parallel for schedule(dynamic, 1)
		#endif*/
		for(int i=0; i<numberOfInputs_; i++){
			Auction<dataType> auction = Auction<dataType>(&current_bidder_diagrams_[i], &barycenter_goods_[i], wasserstein_, geometrical_factor_, 0.01, kdt, correspondance_kdt_map, epsilon, min_diag_price[i]);
			int n_biddings = 0;
			auction.buildUnassignedBidders();
			auction.reinitializeGoods();
			auction.runAuctionRound(n_biddings, i);
			auction.updateDiagonalPrices();
			
			min_diag_price[i] = auction.getMinimalDiagonalPrice();
			std::vector<matchingTuple> matchings;
			dataType cost = auction.getMatchingsAndDistance(&matchings, true);
			all_matchings[i] = matchings;
			total_cost += cost;
			std::cout<< "Barycenter cost for diagram " << i <<" : "<< cost << std::endl;
			std::cout<< "Number of biddings : " << n_biddings << std::endl;
			// Resizes the diagram which was enrich during the auction 
			// TODO do this inside the auction !
			current_bidder_diagrams_[i].bidders_.resize(sizes[i]);
		}
		std::cout<< "Barycenter cost : "<< total_cost << std::endl;

		dataType max_shift = updateBarycenter(all_matchings);
		std::cout<< "Barycenter size : "<< barycenter_goods_[0].size() << std::endl;
		delete kdt;
		//epsilon /= 5;
		dataType eps_candidate = getEpsilon(pow(max_shift, 1./wasserstein_));
		dataType eps_candidate_2 = epsilon/5;
		if(eps_candidate<epsilon){
			epsilon = eps_candidate;
		}
		if(eps_candidate_2>epsilon){
			epsilon = eps_candidate_2;
		}
		if(epsilon>epsilon_0/n_iterations){
			epsilon = epsilon_0/n_iterations;
		}
		
		converged =  (epsilon < 0.0001 * epsilon_0) && (hasBarycenterConverged(all_matchings, previous_matchings) || total_cost>=previous_cost );
		
		previous_matchings = std::move(all_matchings);
		// TODO Correct matchings !
		previous_cost = total_cost;
	}

	for(int j=0; j<barycenter_goods_[0].size(); j++){
		Good<dataType>& g = barycenter_goods_[0].get(j);
		diagramTuple t = std::make_tuple(0, nt1_, 0, nt2_, g.getPersistence(), j, g.x_, 0,0,0, g.y_, 0,0,0);
		barycenter.push_back(t);
	}
		
	return previous_matchings;
}



template <typename dataType>
bool PDBarycenter<dataType>::hasBarycenterConverged(std::vector<std::vector<matchingTuple>>& matchings, std::vector<std::vector<matchingTuple>>& previous_matchings){
	if(points_added_>0 || points_deleted_>0 || previous_matchings.size()==0){
		return false;
	}
	
	for(unsigned int j=0; j<matchings.size(); j++){
		for(unsigned int i=0; i<matchings[j].size(); i++){
			matchingTuple t = matchings[j][i];
			matchingTuple previous_t = previous_matchings[j][i];
			
			if(std::get<1>(t) != std::get<1>(previous_t) && (std::get<0>(t)>=0 && std::get<0>(previous_t)>=0)){
				return false;
			}
		}
	}
	return true;
}


template <typename dataType>
dataType PDBarycenter<dataType>::updateBarycenter(std::vector<std::vector<matchingTuple>>& matchings){
	
	// 1. Initialize variables used in the sequel
	unsigned int n_goods = barycenter_goods_[0].size();
	unsigned int n_diagrams = current_bidder_diagrams_.size();
	points_added_ = 0;
	points_deleted_ = 0;
	dataType max_shift = 0;
	
	std::vector<unsigned int> count_diag_matchings(n_goods);     // Number of diagonal matchings for each point of the barycenter
	std::vector<dataType> x(n_goods);
	std::vector<dataType> y(n_goods);
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
			
			//int pos = current_bidder_ids_[j][bidder_id];
			
			if(good_id<0 && bidder_id>=0){
				// Future new barycenter point
				points_to_append.push_back(&current_bidder_diagrams_[j].get(bidder_id));
			}
			else if(good_id>=0 && bidder_id>=0){
				// Update coordinates (to be divided by the number of diagrams later on)
				x[good_id] += current_bidder_diagrams_[j].get(bidder_id).x_;
				y[good_id] += current_bidder_diagrams_[j].get(bidder_id).y_;
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
			dataType x_bar = x[i]/(double)(n_diagrams - count_diag_matchings[i]);
			dataType y_bar = y[i]/(double)(n_diagrams - count_diag_matchings[i]);
			// 3.2 Compute the new coordinates of the point (the more linked to the diagonal it was, the closer to the diagonal it'll be)
			dataType new_x = ( (double)(n_diagrams - count_diag_matchings[i])*x_bar  + (double)count_diag_matchings[i]*(x_bar+y_bar)/2. )/(double)n_diagrams;
			dataType new_y = ( (double)(n_diagrams - count_diag_matchings[i])*y_bar  + (double)count_diag_matchings[i]*(x_bar+y_bar)/2. )/(double)n_diagrams;
			// 3.3 Compute and store how much the point has shifted
			dataType dx = barycenter_goods_[0].get(i).x_ - new_x;
			dataType dy = barycenter_goods_[0].get(i).y_ - new_y;
			dataType shift = pow(abs(dx), wasserstein_) + pow(abs(dy), wasserstein_);
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
			// TODO Reinitialize/play with prices here if you wish
		}
	}
	
	// 4. Delete off-diagonal barycenter points not linked to any
	// off-diagonal bidder
	for(unsigned int i=0; i<n_goods; i++){
		if(count_diag_matchings[i] == n_diagrams){
			points_deleted_ += 1;
			dataType shift = pow(barycenter_goods_[0].get(i).getPersistence() / pow(2, 1./wasserstein_), wasserstein_);
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
		points_added_ += 1;
		Bidder<dataType>* b = points_to_append[k];
		dataType x = (b->x_ + (n_diagrams -1)*(b->x_+b->y_)/2.)/(n_diagrams); 
		dataType y = (b->y_ + (n_diagrams -1)*(b->x_+b->y_)/2.)/(n_diagrams); 
		for(unsigned int j=0; j<n_diagrams; j++){
			Good<dataType> g = Good<dataType>(x, y, false, barycenter_goods_[j].size());
			g.setPrice(min_prices[j]);
			barycenter_goods_[j].addGood(g);
		}
	}	
	
	// 6. Finally, recreate barycenter_goods
	for(unsigned int j=0; j<n_diagrams; j++){
		int count = 0;
		GoodDiagram<dataType> new_barycenter;
		for(int i=0; i<barycenter_goods_[j].size(); i++){
			Good<dataType>& g = barycenter_goods_[j].get(i);
			if(g.id_!=-1){
				g.id_ = count;
				new_barycenter.addGood(g);
				count ++;
			}
		}
		barycenter_goods_[j] = new_barycenter;
	}
	return max_shift;
}



template <typename dataType>
dataType PDBarycenter<dataType>::getEpsilon(dataType rho){
	return pow(rho, 2)/8;
}

template <typename dataType>
dataType PDBarycenter<dataType>::getRho(dataType epsilon){
	return std::sqrt(8*epsilon);
}



template <typename dataType>
void PDBarycenter<dataType>::setBidderDiagrams(){
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
		current_bidder_diagrams_.push_back(BidderDiagram<dataType>());
		std::vector<int> ids(bidders.size());
		current_bidder_ids_.push_back(ids);
	}
	return;
}


template <typename dataType>
void PDBarycenter<dataType>::enrichCurrentBidderDiagrams(dataType previous_min_persistence, dataType min_persistence){
	for(int i=0; i<numberOfInputs_; i++){
		for(int j=0; j<bidder_diagrams_[i].size(); j++){
			Bidder<dataType> b = bidder_diagrams_[i].get(j);
			dataType persistence = b.getPersistence();
			if(persistence>=min_persistence && persistence<previous_min_persistence){
				b.id_ = current_bidder_diagrams_[i].size();
				b.setPositionInAuction(current_bidder_diagrams_[i].size());
				current_bidder_diagrams_[i].addBidder(b);
				
				// b.id_ --> position of b in current_bidder_diagrams_[i]
				current_bidder_ids_[i][j] = current_bidder_diagrams_[i].size()-1;
			}
		}
	}
}



template <typename dataType>
dataType PDBarycenter<dataType>::getMaxPersistence(){
	dataType max_persistence = 0;
	for(int i=0; i<numberOfInputs_; i++){
		BidderDiagram<dataType>& D = bidder_diagrams_[i];
		for(int j=0; j<D.size(); j++){
			//Add bidder to bidders
			Bidder<dataType>& b = D.get(j);
			dataType persistence = b.getPersistence();
			if(persistence>max_persistence){
				max_persistence = persistence;
			}
		}
	}
	return max_persistence;
}


template <typename dataType>
void PDBarycenter<dataType>::setInitialBarycenter(dataType min_persistence){
	//int random_idx = rand() % numberOfInputs_;
	std::cout << "BEWARE, initial barycenter is not chosen randomly..."<< std::endl;
	int random_idx = 0;
	std::vector<diagramTuple>* CTDiagram = static_cast<std::vector<diagramTuple>*>(inputData_[random_idx]);
	
	for(int i=0; i<numberOfInputs_; i++){
		GoodDiagram<dataType> goods;
		int count=0;
		for(unsigned int j=0; j<CTDiagram->size(); j++){
			//Add bidder to bidders
			Good<dataType> g = Good<dataType>((*CTDiagram)[j], count);
			if(g.getPersistence()>=min_persistence){
				goods.addGood(g);
				count ++;
			}
		}
		barycenter_goods_.push_back(goods);
	}
}


template <typename dataType>
std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> PDBarycenter<dataType>::getKDTree(){
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
