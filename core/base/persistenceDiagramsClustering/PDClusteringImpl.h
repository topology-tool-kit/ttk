#ifndef _PDCLUSTERINGIMPL_H
#define _PDCLUSTERINGIMPL_H

#define BLocalMax ttk::ftm::NodeType::Local_maximum
#define BLocalMin ttk::ftm::NodeType::Local_minimum
#define BSaddle1  ttk::ftm::NodeType::Saddle1
#define BSaddle2  ttk::ftm::NodeType::Saddle2

#include <stdlib.h>     /* srand, rand */
#include <cmath>
#include <PDClustering.h>

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>

using namespace ttk;

template <typename dataType>
int PDClustering<dataType>::execute(){
	Timer t;
	{
	bool converged = false;
	bool diagrams_complete= !use_progressive_;
	n_iterations_ = 0;
	double total_time = 0;
	
	//dataType cost = std::numeric_limits<dataType>::max();
	setBidderDiagrams();
	cost_ = std::numeric_limits<dataType>::max();
	dataType min_cost = std::numeric_limits<dataType>::max();
	epsilon_ = pow(getMostPersistent(), 2)/8.;
	dataType epsilon0 = epsilon_;
	
	// Getting current diagrams (with only at most min_points_to_add points)
	dataType max_persistence = getMostPersistent();
	dataType lowest_persistence = getLessPersistent();
	int min_points_to_add = 30;
	dataType min_persistence = 0;
	if(use_progressive_){
		// min_persistence = max_persistence/2.;
		min_persistence = 0;
	}
	else{
		min_points_to_add = std::numeric_limits<int>::max();
	}
	std::vector<std::vector<dataType>> min_diag_price(3);
	for(int c=0; c<3; ++c){
		for(int i=0; i<numberOfInputs_; i++){
			min_diag_price[c].push_back(0);
		}
	}
	min_persistence = enrichCurrentBidderDiagrams(2*max_persistence, min_persistence, min_diag_price, min_points_to_add);
	
	// Initializing centroids and clusters
	if(use_kmeanspp_){
		initializeCentroidsKMeanspp();
	}
	else{
		initializeCentroids();
	}
	initializeEmptyClusters();
	if(use_accelerated_){
		initializeAcceleratedKMeans();
		getCentroidDistanceMatrix();
		acceleratedUpdateClusters();
	}
	else{
		updateClusters();
	}
	printClustering();
	
	while(!converged || (!diagrams_complete && use_progressive_)){
		Timer t_inside;{
			n_iterations_++;
			dataType max_shift = updateCentroidsPosition();
			if(epsilon_<1e-5){
				converged=true;
			}
			dataType epsilon_candidate = std::max(std::min(max_shift/8., epsilon0/pow(n_iterations_, 2)), epsilon_/5.);
			if(epsilon_candidate<epsilon_){
				epsilon_=epsilon_candidate;
			}
			
			std::cout<< "Iteration "<< n_iterations_<<", Epsilon = "<< epsilon_<< std::endl;
			std::cout<< "Max shift : "<< max_shift << std::endl;
			dataType rho = epsilon_>0 ? std::sqrt(8.0*epsilon_) : -1;
			if(use_progressive_ && n_iterations_>1 && min_persistence>rho && !diagrams_complete){
				if(epsilon_<5e-5){
					// Add all remaining points for final convergence.
					min_persistence = 0;
					min_points_to_add = std::numeric_limits<int>::max();
					diagrams_complete = true;
					use_progressive_=false;
				}
				dataType epsilon_candidate = pow(min_persistence, 2)/8.;
				if(epsilon_candidate>epsilon_){
					// Should always be the case except if min_persistence is equal to zero
					epsilon_ = epsilon_candidate;
				}
				min_diag_price = getMinDiagonalPrices();
				min_persistence = enrichCurrentBidderDiagrams(min_persistence, rho, min_diag_price, min_points_to_add);
				if(min_persistence<lowest_persistence){
					use_progressive_=false;
					diagrams_complete = true;
				}
				// TODO Enrich barycenter using median diagonal and off-diagonal prices
			}
			
			if(use_accelerated_){
				acceleratedUpdateClusters();
			}
			else{
				updateClusters();
			}
			
			if(cost_<min_cost && n_iterations_>2 && !use_progressive_){
				min_cost=cost_;
			}
			else if(n_iterations_>2 && epsilon_<epsilon0/500. && !use_progressive_){
				converged = true;
			}
			std::cout<< "Cost = "<< cost_<< std::endl;
			printClustering();
		}
		total_time +=t_inside.getElapsedTime();
		if(total_time>time_limit_){
			converged = true;
			diagrams_complete = true;
		}
	}
	}// End of timer
	return 0;
}

template <typename dataType>
dataType PDClustering<dataType>::getMostPersistent(){
	dataType max_persistence = 0;
	if(do_min_){
		for(unsigned int i=0; i< bidder_diagrams_min_.size(); ++i){
			for(int j=0; j< bidder_diagrams_min_[i].size(); ++j){
				Bidder<dataType> b = bidder_diagrams_min_[i].get(j);
				dataType persistence = b.getPersistence();
				if(persistence>max_persistence){
					max_persistence = persistence; 
				}
			}
		}
	}
	
	if(do_sad_){
		for(unsigned int i=0; i< bidder_diagrams_saddle_.size(); ++i){
			for(int j=0; j< bidder_diagrams_saddle_[i].size(); ++j){
				Bidder<dataType> b = bidder_diagrams_saddle_[i].get(j);
				dataType persistence = b.getPersistence();
				if(persistence>max_persistence){
					max_persistence = persistence; 
				}
			}
		}
	}
	
	if(do_max_){
		for(unsigned int i=0; i< bidder_diagrams_max_.size(); ++i){
			for(int j=0; j< bidder_diagrams_max_[i].size(); ++j){
				Bidder<dataType> b = bidder_diagrams_max_[i].get(j);
				dataType persistence = b.getPersistence();
				if(persistence>max_persistence){
					max_persistence = persistence; 
				}
			}
		}
	}
	return max_persistence;
}


template <typename dataType>
dataType PDClustering<dataType>::getLessPersistent(){
	dataType min_persistence = 0;
	if(do_min_){
		for(unsigned int i=0; i< bidder_diagrams_min_.size(); ++i){
			for(int j=0; j< bidder_diagrams_min_[i].size(); ++j){
				Bidder<dataType> b = bidder_diagrams_min_[i].get(j);
				dataType persistence = b.getPersistence();
				if(persistence<min_persistence){
					min_persistence = persistence; 
				}
			}
		}
	}
	
	if(do_sad_){
		for(unsigned int i=0; i< bidder_diagrams_saddle_.size(); ++i){
			for(int j=0; j< bidder_diagrams_saddle_[i].size(); ++j){
				Bidder<dataType> b = bidder_diagrams_saddle_[i].get(j);
				dataType persistence = b.getPersistence();
				if(persistence<min_persistence){
					min_persistence = persistence; 
				}
			}
		}
	}
	
	if(do_max_){
		for(unsigned int i=0; i< bidder_diagrams_max_.size(); ++i){
			for(int j=0; j< bidder_diagrams_max_[i].size(); ++j){
				Bidder<dataType> b = bidder_diagrams_max_[i].get(j);
				dataType persistence = b.getPersistence();
				if(persistence<min_persistence){
					min_persistence = persistence; 
				}
			}
		}
	}
	return min_persistence;
}

template <typename dataType>
std::vector<std::vector<dataType>> PDClustering<dataType>::getMinDiagonalPrices(){
	std::vector<std::vector<dataType>> min_prices(3);
	if(do_min_){
		for(unsigned int i=0; i< current_bidder_diagrams_min_.size(); ++i){
			min_prices[0].push_back(std::numeric_limits<dataType>::max());
			for(int j=0; j< current_bidder_diagrams_min_[i].size(); ++j){
				Bidder<dataType> b = current_bidder_diagrams_min_[i].get(j);
				dataType price = b.diagonal_price_;
				if(price<min_prices[0][i]){
					min_prices[0][i] = price; 
				}
			}
			if(min_prices[0][i]>=std::numeric_limits<dataType>::max()/2.){
				min_prices[0][i] = 0;
			}
		}
	}
	
	if(do_sad_){
		for(unsigned int i=0; i< bidder_diagrams_saddle_.size(); ++i){
			min_prices[1].push_back(std::numeric_limits<dataType>::max());
			for(int j=0; j< current_bidder_diagrams_saddle_[i].size(); ++j){
				Bidder<dataType> b = current_bidder_diagrams_saddle_[i].get(j);
				dataType price = b.diagonal_price_;
				if(price<min_prices[1][i]){
					min_prices[1][i] = price; 
				}
			}
			if(min_prices[1][i] >= std::numeric_limits<dataType>::max()/2.){
				min_prices[1][i] = 0;
			}
		}
	}
	
	if(do_max_){
		for(unsigned int i=0; i< current_bidder_diagrams_max_.size(); ++i){
			min_prices[2].push_back(std::numeric_limits<dataType>::max());
			for(int j=0; j< current_bidder_diagrams_max_[i].size(); ++j){
				Bidder<dataType> b = current_bidder_diagrams_max_[i].get(j);
				dataType price = b.diagonal_price_;
				if(price<min_prices[2][i]){
					min_prices[2][i] = price; 
				}
			}
			if(min_prices[2][i] >= std::numeric_limits<dataType>::max()/2.){
				min_prices[2][i] = 0;
			}
		}
	}
	return min_prices;
}


template <typename dataType>
dataType PDClustering<dataType>::computeDistance(BidderDiagram<dataType>& D1, BidderDiagram<dataType>& D2, dataType delta_lim){
	GoodDiagram<dataType> D2_bis = diagramToCentroid(D2);
	return computeDistance(D1, D2_bis, delta_lim);
}

template <typename dataType>
dataType PDClustering<dataType>::computeDistance(BidderDiagram<dataType> D1, GoodDiagram<dataType> D2, dataType delta_lim){
	std::vector<matchingTuple> matchings;
	D2 = centroidWithZeroPrices(D2);
	Auction<dataType> auction(wasserstein_, geometrical_factor_, delta_lim, use_kdtree_);
	auction.BuildAuctionDiagrams(&D1, &D2);
	dataType cost = auction.run(&matchings);
	return cost;
}


template <typename dataType>
dataType PDClustering<dataType>::computeDistance(BidderDiagram<dataType>* D1, GoodDiagram<dataType>* D2, dataType delta_lim){
	std::vector<matchingTuple> matchings;
	Auction<dataType> auction(wasserstein_, geometrical_factor_, delta_lim, use_kdtree_);
	int size1 = D1->size();
	auction.BuildAuctionDiagrams(D1, D2);
	dataType cost = auction.run(&matchings);
	// Diagonal Points were added in the original diagram. The following line removes them.
	D1->bidders_.resize(size1);
	return cost;
}


template <typename dataType>
dataType PDClustering<dataType>::computeDistance(GoodDiagram<dataType>& D1, GoodDiagram<dataType>& D2, dataType delta_lim){
	BidderDiagram<dataType> D1_bis = centroidToDiagram(D1);
	return computeDistance(D1_bis, D2, delta_lim);
}


template <typename dataType>
GoodDiagram<dataType> PDClustering<dataType>::centroidWithZeroPrices(GoodDiagram<dataType> centroid){
	GoodDiagram<dataType> GD = GoodDiagram<dataType>();
	for(int i=0; i<centroid.size(); i++){
		Good<dataType> g = centroid.get(i);
		g.setPrice(0);
		GD.addGood(g);
	}
	return GD;
}


template <typename dataType>
BidderDiagram<dataType> PDClustering<dataType>::diagramWithZeroPrices(BidderDiagram<dataType> diagram){
	BidderDiagram<dataType> BD = BidderDiagram<dataType>();
	for(int i=0; i<diagram.size(); i++){
		Bidder<dataType> b = diagram.get(i);
		b.setDiagonalPrice(0);
		BD.addBidder(b);
	}
	return BD;
}


template <typename dataType>
BidderDiagram<dataType> PDClustering<dataType>::centroidToDiagram(GoodDiagram<dataType> centroid){
	BidderDiagram<dataType> BD = BidderDiagram<dataType>();
	for(int i=0; i<centroid.size(); i++){
		Good<dataType> g = centroid.get(i);
		
		Bidder<dataType> b = Bidder<dataType>(g.x_, g.y_, g.isDiagonal(), BD.size());
		b.SetCriticalCoordinates(g.coords_x_, g.coords_y_, g.coords_z_);
		b.setPositionInAuction(BD.size());
		BD.addBidder(b);
	}
	return BD;
}

template <typename dataType>
GoodDiagram<dataType> PDClustering<dataType>::diagramToCentroid(BidderDiagram<dataType> diagram){
	GoodDiagram<dataType> GD = GoodDiagram<dataType>();
	for(int i=0; i<diagram.size(); i++){
		Bidder<dataType> b = diagram.get(i);
		
		Good<dataType> g = Good<dataType>(b.x_, b.y_, b.isDiagonal(), GD.size());
		g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
		GD.addGood(g);
	}
	return GD;
}


template <typename dataType>
void PDClustering<dataType>::initializeEmptyClusters(){
	clustering_ = std::vector<std::vector<int>>(k_);
}

template <typename dataType>
void PDClustering<dataType>::initializeCentroids(){
	std::vector<int> idx(numberOfInputs_);
	// To perform a random draw with replacement, the vector {1, 2, ..., numberOfInputs_} is
	// shuffled, and we consider its k_ first elements to be the initial centroids.
	for(int i=0; i<numberOfInputs_; i++){
		idx[i] = i;
	}
	std::random_shuffle(idx.begin(), idx.end());
	
	for(int c=0; c<k_; c++){
		if(do_min_){
			GoodDiagram<dataType> centroid_min = diagramToCentroid(current_bidder_diagrams_min_[idx[c]]);
			centroids_min_.push_back(centroid_min);
		}
		if(do_sad_){
			GoodDiagram<dataType> centroid_sad = diagramToCentroid(current_bidder_diagrams_saddle_[idx[c]]);
			centroids_saddle_.push_back(centroid_sad);
		}
		if(do_max_){
			GoodDiagram<dataType> centroid_max = diagramToCentroid(current_bidder_diagrams_max_[idx[c]]);
			centroids_max_.push_back(centroid_max);
		}
	}
}

template <typename dataType>
void PDClustering<dataType>::initializeCentroidsKMeanspp(){
	std::vector<int> indexes_clusters;
	int random_idx = rand() % numberOfInputs_;
	indexes_clusters.push_back(random_idx);
	
	if(do_min_){
		GoodDiagram<dataType> centroid_min = diagramToCentroid(current_bidder_diagrams_min_[random_idx]);
		centroids_min_.push_back(centroid_min);
	}
	if(do_sad_){
		GoodDiagram<dataType> centroid_sad = diagramToCentroid(current_bidder_diagrams_saddle_[random_idx]);
		centroids_saddle_.push_back(centroid_sad);
	}
	if(do_max_){
		GoodDiagram<dataType> centroid_max = diagramToCentroid(current_bidder_diagrams_max_[random_idx]);
		centroids_max_.push_back(centroid_max);
	}
	
	while((int) indexes_clusters.size()<k_){
		std::vector<dataType> min_distance_to_centroid(numberOfInputs_);
		std::vector<dataType> probabilities(numberOfInputs_);
		// Uncomment for a deterministic algorithm
		// dataType maximal_distance = 0;
		//int candidate_centroid = -1;
		for(int i=0; i<numberOfInputs_; i++){
			min_distance_to_centroid[i] = std::numeric_limits<dataType>::max();
			if(std::find(indexes_clusters.begin(), indexes_clusters.end(), i) != indexes_clusters.end()){
				min_distance_to_centroid[i] = 0;
			}
			else{
				for(unsigned int j=0; j<indexes_clusters.size(); ++j){
					dataType distance = 0;
					if(do_min_){
						GoodDiagram<dataType> centroid_min = centroidWithZeroPrices(centroids_min_[j]);
						distance += computeDistance(current_bidder_diagrams_min_[i], centroid_min);
					}
					if(do_sad_){
						GoodDiagram<dataType> centroid_saddle = centroidWithZeroPrices(centroids_saddle_[j]);
						distance += computeDistance(current_bidder_diagrams_saddle_[i], centroid_saddle);
					}
					if(do_max_){
						GoodDiagram<dataType> centroid_max = centroidWithZeroPrices(centroids_max_[j]);
						distance += computeDistance(current_bidder_diagrams_max_[i], centroid_max);
					}
					if(distance<min_distance_to_centroid[i]){
						min_distance_to_centroid[i] = distance;
					}
				}
			}
			probabilities[i] = pow(min_distance_to_centroid[i], 2);
			
			// The following block is useful in case of need for a deterministic algoritm
			/*if(min_distance_to_centroid[i]>maximal_distance){
				maximal_distance = min_distance_to_centroid[i];
				candidate_centroid = i;
			}*/
			
		}
		// Comment the following four lines to make it deterministic
		std::random_device rd;
		std::mt19937 gen(rd());
		std::discrete_distribution<int> distribution (probabilities.begin(),probabilities.end());
		int candidate_centroid = distribution(gen);
		
		indexes_clusters.push_back(candidate_centroid);
		if(do_min_){
			GoodDiagram<dataType> centroid_min = diagramToCentroid(current_bidder_diagrams_min_[candidate_centroid]);
			centroids_min_.push_back(centroid_min);
		}
		if(do_sad_){
			GoodDiagram<dataType> centroid_sad = diagramToCentroid(current_bidder_diagrams_saddle_[candidate_centroid]);
			centroids_saddle_.push_back(centroid_sad);
		}
		if(do_max_){
			GoodDiagram<dataType> centroid_max = diagramToCentroid(current_bidder_diagrams_max_[candidate_centroid]);
			centroids_max_.push_back(centroid_max);
		}
	}
}

template <typename dataType>
void PDClustering<dataType>::initializeAcceleratedKMeans(){
	// r_ is a vector stating for each diagram if its distance to its centroid is
	// up to date (false) or needs to be recomputed (true)
	r_ = std::vector<bool>(numberOfInputs_);
	// u_ is a vector of upper bounds of the distance of each diagram to its closest centroid
	u_ = std::vector<dataType>(numberOfInputs_);
	inv_clustering_ = std::vector<int>(numberOfInputs_);
	for(int i=0; i<numberOfInputs_; i++){
		r_[i] = true;
		u_[i] = std::numeric_limits<dataType>::max();
		inv_clustering_[i] = -1;
	}
	// l_ is the matrix of lower bounds for the distance from each diagram
	// to each centroid
	l_ = std::vector<std::vector<dataType>>(numberOfInputs_);
	for(int i=0; i<numberOfInputs_; ++i){
		l_[i] = std::vector<dataType>(k_);
		for(int c=0; c<k_; ++c){
			l_[i][c]=0;
		}
	}
	
	// And d_ is a (K x K) matrix storing the distances between each pair of centroids 
	d_ = std::vector<std::vector<dataType>>(k_);
	for(int i=0; i<k_; ++i){
		for(int c=0; c<k_; ++c){
			d_[i].push_back(0);
		}
	}
	return ;
}


template <typename dataType>
std::vector<std::vector<dataType>> PDClustering<dataType>::getDistanceMatrix(){
	std::vector<std::vector<dataType>> D(numberOfInputs_);
	
	for(int i=0; i<numberOfInputs_; ++i){
		BidderDiagram<dataType> D1_min, D1_sad, D1_max;
		if(do_min_){
			D1_min = diagramWithZeroPrices(current_bidder_diagrams_min_[i]);
		}
		if(do_sad_){
			D1_sad = diagramWithZeroPrices(current_bidder_diagrams_saddle_[i]);
		}
		if(do_max_){
			D1_max = diagramWithZeroPrices(current_bidder_diagrams_max_[i]);
		}
		for(int c=0; c<k_; ++c){
			GoodDiagram<dataType> D2_min, D2_sad, D2_max;
			dataType distance = 0;
			if(do_min_){
				D2_min = centroids_min_[c];
				distance += computeDistance(D1_min, D2_min, 0.01);
			}
			if(do_sad_){
				D2_sad = centroids_saddle_[c];
				distance += computeDistance(D1_sad, D2_sad, 0.01);
			}
			if(do_max_){
				D2_max = centroids_max_[c];
				distance += computeDistance(D1_max, D2_max, 0.01);
			}
			D[i].push_back(distance);
		}
	}
	return D;
}

template <typename dataType>
void PDClustering<dataType>::getCentroidDistanceMatrix(){
	for(int i=0; i<k_; ++i){
		GoodDiagram<dataType> D1_min, D1_sad, D1_max;
		if(do_min_){
			D1_min = centroidWithZeroPrices(centroids_min_[i]);
		}
		if(do_sad_){
			D1_sad = centroidWithZeroPrices(centroids_saddle_[i]);
		}
		if(do_max_){
			D1_max = centroidWithZeroPrices(centroids_max_[i]);
		}
		for(int j=i+1; j<k_; ++j){
			dataType distance = 0;
			GoodDiagram<dataType> D2_min, D2_sad, D2_max;
			if(do_min_){
				D2_min = centroidWithZeroPrices(centroids_min_[j]);
				distance += computeDistance(D1_min, D2_min, 0.001);
			}
			if(do_sad_){
				D2_sad = centroidWithZeroPrices(centroids_saddle_[j]);
				distance += computeDistance(D1_sad, D2_sad, 0.001);
			}
			if(do_max_){
				D2_max = centroidWithZeroPrices(centroids_max_[j]);
				distance += computeDistance(D1_max, D2_max, 0.001);
			}
			
			d_[i][j] = distance;
			d_[j][i] = distance;
		}
	}
	return;
}


template <typename dataType>
void PDClustering<dataType>::updateClusters(){
	std::vector<std::vector<dataType>> distance_matrix = getDistanceMatrix();
	old_clustering_ = clustering_;
	invertClusters();
	initializeEmptyClusters();
	
	for(int i=0; i<numberOfInputs_; ++i){	
		dataType min_distance_to_centroid = std::numeric_limits<dataType>::max();
		int cluster = -1;
		for(int c=0; c<k_; ++c){
			if(distance_matrix[i][c]<min_distance_to_centroid){
				min_distance_to_centroid = distance_matrix[i][c];
				cluster = c;
			}
		}
		
		clustering_[cluster].push_back(i);
		if(cluster!=inv_clustering_[i]){
			// New centroid attributed to this diagram
			if(do_min_){
				centroids_with_price_min_[i] = centroidWithZeroPrices(centroids_min_[cluster]);
			}
			if(do_sad_){
				centroids_with_price_saddle_[i] = centroidWithZeroPrices(centroids_saddle_[cluster]);
			}
			if(do_max_){
				centroids_with_price_max_[i] = centroidWithZeroPrices(centroids_max_[cluster]);
			}
			inv_clustering_[i] = cluster;
		}
	}
	for(int c=0; c<k_; ++c){
		if(clustering_[c].size()==0){
			clustering_[c].push_back( rand() % numberOfInputs_);
		}
	}
	return ;
}

template <typename dataType>
void PDClustering<dataType>::invertClusters(){
	/// Converts the clustering (vector of vector of diagram's id) into
	/// a vector of size numberOfInputs_ containg the cluster of each input diagram. 
	
	// Initializes clusters with -1
	inv_clustering_ = std::vector<int>(numberOfInputs_);
	for(int i=0; i<numberOfInputs_; ++i){
		inv_clustering_[i] = -1;
	}
	
	// Fill in the clusters
	for(int c=0; c<k_; ++c){
		for(unsigned int j=0; j<clustering_[c].size(); ++j){
			int idx = clustering_[c][j];
			inv_clustering_[idx] = c;
		}
	}
}

template <typename dataType>
void PDClustering<dataType>::invertInverseClusters(){
	clustering_ = std::vector<std::vector<int>>(k_);
	for(int i=0; i<numberOfInputs_; ++i){
		clustering_[inv_clustering_[i]].push_back(i);
	}
	
	// Check if a cluster was left without diagram
	for(int c=0; c<k_; ++c){
		if(clustering_[c].size() == 0){
			std::cout<< "Problem in invertInverseClusters()... \nCluster " << c << " was left with no diagram attached to it... " << std::endl;
		}
	}
}


template <typename dataType>
void PDClustering<dataType>::acceleratedUpdateClusters(){
	// Step 1
	getCentroidDistanceMatrix();
	old_clustering_ = clustering_;
	//self.old_clusters = copy.copy(self.clusters)
	invertClusters();
	initializeEmptyClusters();
	
	for(int i=0; i<numberOfInputs_; ++i){
		// Step 3 find potential changes of clusters
		BidderDiagram<dataType> D1_min, D1_sad, D1_max;
		if(do_min_){
			D1_min = diagramWithZeroPrices(current_bidder_diagrams_min_[i]);
		}
		if(do_sad_){
			D1_sad = diagramWithZeroPrices(current_bidder_diagrams_saddle_[i]);
		}
		if(do_max_){
			D1_max = diagramWithZeroPrices(current_bidder_diagrams_max_[i]);
		}
		
		for(int c=0; c<k_; ++c){
			if(inv_clustering_[i]==-1){
				// If not yet assigned, assign it first to a random cluster
				inv_clustering_[i] = rand() % (k_);
				r_[i] = true;
				if(do_min_){
					centroids_with_price_min_[i] = centroidWithZeroPrices(centroids_min_[inv_clustering_[i]]);
				}
				if(do_sad_){
					centroids_with_price_saddle_[i] = centroidWithZeroPrices(centroids_saddle_[inv_clustering_[i]]);
				}
				if(do_max_){
					centroids_with_price_max_[i] = centroidWithZeroPrices(centroids_max_[inv_clustering_[i]]);
				}
			}
			
			if(c!=inv_clustering_[i] && u_[i]>l_[i][c] && u_[i]>0.5*d_[inv_clustering_[i]][c]){
				// Step 3a, If necessary, recompute the distance to centroid
				if(r_[i]){
					dataType distance = 0;
					GoodDiagram<dataType> centroid_min, centroid_sad, centroid_max;
					if(do_min_){
						centroid_min = centroidWithZeroPrices(centroids_min_[inv_clustering_[i]]);
						distance += computeDistance(D1_min, centroid_min, 0.01);
					}
					if(do_sad_){
						centroid_sad = centroidWithZeroPrices(centroids_saddle_[inv_clustering_[i]]);
						distance += computeDistance(D1_sad, centroid_sad, 0.01);
					}
					if(do_max_){
						centroid_max = centroidWithZeroPrices(centroids_max_[inv_clustering_[i]]);
						distance += computeDistance(D1_max, centroid_max, 0.01);
					}
					r_[i] = false;
					u_[i] = distance;
					l_[i][inv_clustering_[i]] = distance;
				}
				// Step 3b, check if still potential change of clusters
				if(u_[i]>l_[i][c] || u_[i]>0.5*d_[inv_clustering_[i]][c]){
					BidderDiagram<dataType> diagram_min, diagram_sad, diagram_max;
					GoodDiagram<dataType> centroid_min, centroid_sad, centroid_max;
					dataType distance = 0;
					
					if(do_min_){
						centroid_min = centroidWithZeroPrices(centroids_min_[c]);
						diagram_min = diagramWithZeroPrices(current_bidder_diagrams_min_[i]);
						distance += computeDistance(diagram_min, centroid_min, 0.01);
					}
					if(do_sad_){
						centroid_sad = centroidWithZeroPrices(centroids_saddle_[c]);
						diagram_sad = diagramWithZeroPrices(current_bidder_diagrams_saddle_[i]);
						distance += computeDistance(diagram_sad, centroid_sad, 0.01);
					}
					if(do_max_){
						centroid_max = centroidWithZeroPrices(centroids_max_[c]);
						diagram_max = diagramWithZeroPrices(current_bidder_diagrams_max_[i]);
						distance += computeDistance(diagram_max, centroid_max, 0.01);
					}
					l_[i][c] = distance;
					// TODO Prices are lost here... If distance<self.u[i], we should keep the prices
					if(distance<u_[i]){
						// Changing cluster 
						u_[i] = distance;
						inv_clustering_[i] = c;
						
						if(do_min_){
							centroids_with_price_min_[i] = centroidWithZeroPrices(centroids_min_[c]);
						}
						if(do_sad_){
							centroids_with_price_saddle_[i] = centroidWithZeroPrices(centroids_saddle_[c]);
						}
						if(do_max_){
							centroids_with_price_max_[i] = centroidWithZeroPrices(centroids_max_[c]);
						}
					}
				}
			}
		}
	}
	invertInverseClusters();
	for(int c=0; c<k_; ++c){
		if(clustering_[c].size()==0){
			std::cout<< "Adding artificial centroid because a cluster was empty" <<std::endl;
			bool idx_acceptable = false;
			int idx;
			while(!idx_acceptable){
				idx = rand() % k_;
				if(inv_clustering_[idx]<k_ && inv_clustering_[idx]>=0  && clustering_[inv_clustering_[idx]].size()>1){
					idx_acceptable = true;
					int cluster_removal = inv_clustering_[idx];
					// Removing the index to remove
					clustering_[cluster_removal].erase(std::remove(clustering_[cluster_removal].begin(), clustering_[cluster_removal].end(), idx), clustering_[cluster_removal].end());
				}
				
			}
			
			clustering_[c].push_back(idx);
			inv_clustering_[idx] = c;
			
			if(do_min_){
				centroids_min_[c] = diagramToCentroid(current_bidder_diagrams_min_[idx]);
				centroids_with_price_min_[idx] = centroidWithZeroPrices(centroids_min_[c]);
			}
			if(do_sad_){
				centroids_saddle_[c] = diagramToCentroid(current_bidder_diagrams_saddle_[idx]);
				centroids_with_price_saddle_[idx] = centroidWithZeroPrices(centroids_saddle_[c]);
			}
			if(do_max_){
				centroids_max_[c] = diagramToCentroid(current_bidder_diagrams_max_[idx]);
				centroids_with_price_max_[idx] = centroidWithZeroPrices(centroids_max_[c]);
			}
		}
	}
	return ;
}

template <typename dataType>
dataType PDClustering<dataType>::updateCentroidsPosition(){	 
	dataType max_shift = 0;
	dataType max_wasserstein_shift = 0;
	cost_ = 0;
	for(int c=0; c<k_; ++c){
		std::vector<GoodDiagram<dataType> > centroids_with_price_min, centroids_with_price_sad, centroids_with_price_max;
		int count = 0;
		for(int idx : clustering_[c]){
			int number_of_points = 0;
			// Find the position of diagrams[idx] in old cluster c
			vector<int>::iterator i = std::find(old_clustering_[c].begin(), old_clustering_[c].end (), idx);
			int pos = (i==old_clustering_[c].end()) ? -1 : std::distance(old_clustering_[c].begin(), i);
			if(pos>=0){
				// Diagram was already linked to this centroid before
				if(do_min_){
					centroids_with_price_min.push_back(centroids_with_price_min_[idx]);
					number_of_points += centroids_with_price_min_[idx].size() + current_bidder_diagrams_min_[idx].size();
				}
				if(do_sad_){
					centroids_with_price_sad.push_back(centroids_with_price_saddle_[idx]);
					number_of_points += centroids_with_price_saddle_[idx].size() + current_bidder_diagrams_saddle_[idx].size();
				}
				if(do_max_){
					centroids_with_price_max.push_back(centroids_with_price_max_[idx]);
					number_of_points += centroids_with_price_max_[idx].size() + current_bidder_diagrams_max_[idx].size();
				}
			}
			else{
				// Otherwise, centroid is given 0 prices and the diagram is given 0 diagonal-prices
				if(do_min_){
					centroids_with_price_min.push_back(centroidWithZeroPrices( centroids_min_[c] ));
					current_bidder_diagrams_min_[idx] = diagramWithZeroPrices(current_bidder_diagrams_min_[idx]);
					number_of_points += centroids_with_price_min_[idx].size() + current_bidder_diagrams_min_[idx].size();
				}
				if(do_sad_){
					centroids_with_price_sad.push_back(centroidWithZeroPrices( centroids_saddle_[c] ));
					current_bidder_diagrams_saddle_[idx] = diagramWithZeroPrices(current_bidder_diagrams_saddle_[idx]);
					number_of_points += centroids_with_price_saddle_[idx].size() + current_bidder_diagrams_saddle_[idx].size();
				}
				if(do_max_){
					centroids_with_price_max.push_back(centroidWithZeroPrices( centroids_max_[c] ));
					current_bidder_diagrams_max_[idx] = diagramWithZeroPrices(current_bidder_diagrams_max_[idx]);
					number_of_points += centroids_with_price_max_[idx].size() + current_bidder_diagrams_max_[idx].size();
				}
				
				if(n_iterations_>1){
					// If diagram new to cluster and we're not at first iteration, 
					// precompute prices for the objects via compute_distance()
					number_of_points /= (int) do_min_ + (int) do_sad_ +  (int) do_max_;
					dataType d_estimated = pow(cost_/numberOfInputs_, 1./wasserstein_)+ 1e-7;
					dataType estimated_delta_lim = 2 * number_of_points * epsilon_ / d_estimated;
					
					if(estimated_delta_lim>1){
						estimated_delta_lim=1;
					}
					
					// We use pointer in the auction in order to keep the prices at the end
					if(do_min_){
						computeDistance(&(current_bidder_diagrams_min_[idx]), &(centroids_with_price_min[count]), estimated_delta_lim);
					}
					if(do_sad_){
						computeDistance(&(current_bidder_diagrams_saddle_[idx]), &(centroids_with_price_sad[count]), estimated_delta_lim);
					}
					if(do_max_){
						computeDistance(&(current_bidder_diagrams_max_[idx]), &(centroids_with_price_max[count]), estimated_delta_lim);
					}
					
				}
			}
			count++;
		}
		std::vector<BidderDiagram<dataType>> diagrams_c_min, diagrams_c_sad, diagrams_c_max;
		dataType total_cost = 0;
		dataType wasserstein_shift = 0;
		if(do_min_){
			
			for(int idx : clustering_[c]){
				diagrams_c_min.push_back(current_bidder_diagrams_min_[idx]);
			}
			
			std::vector<dataType> min_diag_price(diagrams_c_min.size());
			std::vector<dataType> min_price(diagrams_c_min.size());
			for(unsigned int i=0; i<diagrams_c_min.size(); i++){
				min_diag_price[i] = 0;
				min_price[i] = 0;
			}
			std::vector<int> sizes(diagrams_c_min.size());
			for(unsigned int i=0; i<diagrams_c_min.size(); i++){
				sizes[i] = diagrams_c_min[i].size();
			}
			
			PDBarycenter<dataType> barycenter_computer = PDBarycenter<dataType>();
			barycenter_computer.setThreadNumber(threadNumber_);
			barycenter_computer.setWasserstein(wasserstein_);
			barycenter_computer.setNumberOfInputs(diagrams_c_min.size());
			barycenter_computer.setDiagramType(0);
			barycenter_computer.setUseProgressive(false);
			barycenter_computer.setGeometricalFactor(geometrical_factor_);
			
			barycenter_computer.setCurrentBidders(diagrams_c_min);
			barycenter_computer.setCurrentBarycenter(centroids_with_price_min);
			std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> pair;
			bool use_kdt = false;
			if(centroids_with_price_min[0].size()>0){
				pair = barycenter_computer.getKDTree();
				use_kdt=true;
			}
			KDTree<dataType>* kdt = pair.first;
			std::vector<KDTree<dataType>*>& correspondance_kdt_map = pair.second;
			
			
			std::vector<std::vector<matchingTuple>> all_matchings(diagrams_c_min.size());
			
			barycenter_computer.runMatching(&total_cost, 
										epsilon_,
										sizes,
										kdt, 
										&correspondance_kdt_map, 
										&min_diag_price, 
										&min_price,
										&all_matchings,
										use_kdt);
			dataType max_shift_c_min = barycenter_computer.updateBarycenter(all_matchings);
			if(max_shift_c_min>max_shift){
				max_shift = max_shift_c_min;
			}
			
			// Now that barycenters and diagrams are updated in PDBarycenter class,
			// we import the results here.
			diagrams_c_min = barycenter_computer.getCurrentBidders();
			centroids_with_price_min = barycenter_computer.getCurrentBarycenter();
			int i = 0;
			for(int idx : clustering_[c]){
				current_bidder_diagrams_min_[idx] = diagrams_c_min[i];
				centroids_with_price_min_[idx] = centroids_with_price_min[i];
				i++;
			}
			
			GoodDiagram<dataType> old_centroid = centroids_min_[c];
			centroids_min_[c] = centroidWithZeroPrices(centroids_with_price_min_[clustering_[c][0]]);
			wasserstein_shift += computeDistance(old_centroid, centroids_min_[c], 0.01);
		}
		
		if(do_sad_){
			
			for(int idx : clustering_[c]){
				diagrams_c_sad.push_back(current_bidder_diagrams_saddle_[idx]);
			}
			
			std::vector<dataType> min_diag_price(diagrams_c_sad.size());
			std::vector<dataType> min_price(diagrams_c_sad.size());
			for(unsigned int i=0; i<diagrams_c_sad.size(); i++){
				min_diag_price[i] = 0;
				min_price[i] = 0;
			}
			std::vector<int> sizes(diagrams_c_sad.size());
			for(unsigned int i=0; i<diagrams_c_sad.size(); i++){
				sizes[i] = diagrams_c_sad[i].size();
			}
			
			PDBarycenter<dataType> barycenter_computer = PDBarycenter<dataType>();
			barycenter_computer.setThreadNumber(threadNumber_);
			barycenter_computer.setWasserstein(wasserstein_);
			barycenter_computer.setNumberOfInputs(diagrams_c_sad.size());
			barycenter_computer.setDiagramType(1);
			barycenter_computer.setUseProgressive(false);
			barycenter_computer.setGeometricalFactor(geometrical_factor_);
			
			barycenter_computer.setCurrentBidders(diagrams_c_sad);
			barycenter_computer.setCurrentBarycenter(centroids_with_price_sad);
			std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> pair;
			bool use_kdt = false;
			if(centroids_with_price_sad[0].size()>0){
				pair = barycenter_computer.getKDTree();
				use_kdt=true;
			}
			KDTree<dataType>* kdt = pair.first;
			std::vector<KDTree<dataType>*>& correspondance_kdt_map = pair.second;
			
			
			std::vector<std::vector<matchingTuple>> all_matchings(diagrams_c_sad.size());
			
			barycenter_computer.runMatching(&total_cost, 
										epsilon_,
										sizes,
										kdt, 
										&correspondance_kdt_map, 
										&min_diag_price, 
										&min_price,
										&all_matchings,
										use_kdt);
			dataType max_shift_c_sad = barycenter_computer.updateBarycenter(all_matchings);
			if(max_shift_c_sad>max_shift){
				max_shift = max_shift_c_sad;
			}
			
			// Now that barycenters and diagrams are updated in PDBarycenter class,
			// we import the results here.
			diagrams_c_sad = barycenter_computer.getCurrentBidders();
			centroids_with_price_sad = barycenter_computer.getCurrentBarycenter();
			int i = 0;
			for(int idx : clustering_[c]){
				current_bidder_diagrams_saddle_[idx] = diagrams_c_sad[i];
				centroids_with_price_saddle_[idx] = centroids_with_price_sad[i];
				i++;
			}
			GoodDiagram<dataType> old_centroid = centroids_saddle_[c];
			centroids_saddle_[c] = centroidWithZeroPrices(centroids_with_price_saddle_[clustering_[c][0]]);
			wasserstein_shift += computeDistance(old_centroid, centroids_saddle_[c], 0.01);
		}
		
		if(do_max_){
			
			for(int idx : clustering_[c]){
				diagrams_c_max.push_back(current_bidder_diagrams_max_[idx]);
			}
			std::vector<dataType> min_diag_price(diagrams_c_max.size());
			std::vector<dataType> min_price(diagrams_c_max.size());
			for(unsigned int i=0; i<diagrams_c_max.size(); i++){
				min_diag_price[i] = 0;
				min_price[i] = 0;
			}
			std::vector<int> sizes(diagrams_c_max.size());
			for(unsigned int i=0; i<diagrams_c_max.size(); i++){
				sizes[i] = diagrams_c_max[i].size();
			}
			
			PDBarycenter<dataType> barycenter_computer = PDBarycenter<dataType>();
			barycenter_computer.setThreadNumber(threadNumber_);
			barycenter_computer.setWasserstein(wasserstein_);
			barycenter_computer.setNumberOfInputs(diagrams_c_max.size());
			barycenter_computer.setDiagramType(2);
			barycenter_computer.setUseProgressive(false);
			barycenter_computer.setGeometricalFactor(geometrical_factor_);
			
			barycenter_computer.setCurrentBidders(diagrams_c_max);
			barycenter_computer.setCurrentBarycenter(centroids_with_price_max);
			std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> pair;
			bool use_kdt = false;
			if(centroids_with_price_max[0].size()>0){
				pair = barycenter_computer.getKDTree();
				use_kdt=true;
			}
			KDTree<dataType>* kdt = pair.first;
			std::vector<KDTree<dataType>*>& correspondance_kdt_map = pair.second;
			
			std::vector<std::vector<matchingTuple>> all_matchings(diagrams_c_max.size());
			
			barycenter_computer.runMatching(&total_cost, 
										epsilon_,
										sizes,
										kdt, 
										&correspondance_kdt_map, 
										&min_diag_price, 
										&min_price,
										&all_matchings,
										use_kdt);
			dataType max_shift_c_max = barycenter_computer.updateBarycenter(all_matchings);
			if(max_shift_c_max>max_shift){
				max_shift = max_shift_c_max;
			}
										
			// Now that barycenters and diagrams are updated in PDBarycenter class,
			// we import the results here.
			diagrams_c_max = barycenter_computer.getCurrentBidders();
			centroids_with_price_max = barycenter_computer.getCurrentBarycenter();
			int i = 0;
			for(int idx : clustering_[c]){
				current_bidder_diagrams_max_[idx] = diagrams_c_max[i];
				centroids_with_price_max_[idx] = centroids_with_price_max[i];
				i++;
			}
			GoodDiagram<dataType> old_centroid = centroids_max_[c];
			centroids_max_[c] = centroidWithZeroPrices(centroids_with_price_max_[clustering_[c][0]]);
			wasserstein_shift += computeDistance(old_centroid, centroids_max_[c], 0.01);
		}
		cost_+=total_cost;
		if(wasserstein_shift>max_wasserstein_shift){
			max_wasserstein_shift = wasserstein_shift;
		}
		if(use_accelerated_){
			for(int i=0; i< numberOfInputs_; ++i){
			// Step 5 of Accelerated KMeans: Update the lower bound on distance thanks to the triangular inequality
				l_[i][c] = pow( pow(l_[i][c], 1./wasserstein_) - pow(wasserstein_shift, 1./wasserstein_), wasserstein_);
				if(l_[i][c]<0){
					l_[i][c] = 0;
				}	
			}
			for(int idx : clustering_[c]){
				// Step 6, update the upper bound on the distance to the centroid thanks to the triangle inequality
				u_[idx] = pow( pow(u_[idx], 1./wasserstein_) + pow(wasserstein_shift, 1./wasserstein_), wasserstein_);
				r_[idx] = true;
			}
		}
	}
	// Normally return max_shift, but it seems there is a bug 
	// yielding max_shift > 100 * max_wasserstein_shift
	// which should logically not really happen...
	// This is supposed to be only a temporary patch...
	return std::min(max_shift, max_wasserstein_shift);
}


template <typename dataType>
void PDClustering<dataType>::setBidderDiagrams(){
  
	for(int i=0; i<numberOfInputs_; i++){
		if(do_min_){
			std::vector<diagramTuple> *CTDiagram = &((*inputDiagramsMin_)[i]);
			
			BidderDiagram<dataType> bidders;
			for(unsigned int j=0; j<CTDiagram->size(); j++){
				//Add bidder to bidders
				Bidder<dataType> b((*CTDiagram)[j], j);

				b.setPositionInAuction(bidders.size());
				bidders.addBidder(b);
				if(b.isDiagonal() || b.x_==b.y_){
					std::cout<<"Diagonal point in diagram !!!"<<std::endl;
				}
			}
			bidder_diagrams_min_.push_back(bidders);
			current_bidder_diagrams_min_.push_back(BidderDiagram<dataType>());
			centroids_with_price_min_.push_back(GoodDiagram<dataType>());
		}
		
		if(do_sad_){
			std::vector<diagramTuple> *CTDiagram = &((*inputDiagramsSaddle_)[i]);
			
			BidderDiagram<dataType> bidders;
			for(unsigned int j=0; j<CTDiagram->size(); j++){
				//Add bidder to bidders
				Bidder<dataType> b((*CTDiagram)[j], j);

				b.setPositionInAuction(bidders.size());
				bidders.addBidder(b);
				if(b.isDiagonal() || b.x_==b.y_){
					std::cout<<"Diagonal point in diagram !!!"<<std::endl;
				}
			}
			bidder_diagrams_saddle_.push_back(bidders);
			current_bidder_diagrams_saddle_.push_back(BidderDiagram<dataType>());
			centroids_with_price_saddle_.push_back(GoodDiagram<dataType>());
		}
		
		if(do_max_){
			std::vector<diagramTuple> *CTDiagram = &((*inputDiagramsMax_)[i]);
			
			BidderDiagram<dataType> bidders;
			for(unsigned int j=0; j<CTDiagram->size(); j++){
				//Add bidder to bidders
				Bidder<dataType> b((*CTDiagram)[j], j);

				b.setPositionInAuction(bidders.size());
				bidders.addBidder(b);
				if(b.isDiagonal() || b.x_==b.y_){
					std::cout<<"Diagonal point in diagram !!!"<<std::endl;
				}
			}
			bidder_diagrams_max_.push_back(bidders);
			current_bidder_diagrams_max_.push_back(BidderDiagram<dataType>());
			centroids_with_price_max_.push_back(GoodDiagram<dataType>());
		}
		
	}
	return;
}


template <typename dataType>
dataType PDClustering<dataType>::enrichCurrentBidderDiagrams(dataType previous_min_persistence, dataType min_persistence, std::vector<std::vector<dataType>> initial_diagonal_prices, int min_points_to_add){
	
  dataType new_min_persistence = min_persistence;
	
  // 1. Get size of the largest current diagram, deduce the maximal number of 
  // points to append
	int max_diagram_size=0;
	if(do_min_){
		for(int i=0; i<numberOfInputs_; i++){
			if(current_bidder_diagrams_min_[i].size()>max_diagram_size){
				max_diagram_size = current_bidder_diagrams_min_[i].size();
			}
		}
	}
	if(do_sad_){
		for(int i=0; i<numberOfInputs_; i++){
			if(current_bidder_diagrams_saddle_[i].size()>max_diagram_size){
				max_diagram_size = current_bidder_diagrams_saddle_[i].size();
			}
		}
	}
	if(do_max_){
		for(int i=0; i<numberOfInputs_; i++){
			if(current_bidder_diagrams_max_[i].size()>max_diagram_size){
				max_diagram_size = current_bidder_diagrams_max_[i].size();
			}
		}
	}
	int max_points_to_add = std::max(min_points_to_add, min_points_to_add + (int) (max_diagram_size/10));
	
	// 2. Get which points can be added, deduce the new minimal persistence
	std::vector<std::vector<int>> candidates_to_be_added_min(numberOfInputs_);
	std::vector<std::vector<int>> candidates_to_be_added_sad(numberOfInputs_);
	std::vector<std::vector<int>> candidates_to_be_added_max(numberOfInputs_);
	std::vector<std::vector<int>> idx_min(numberOfInputs_);
	std::vector<std::vector<int>> idx_sad(numberOfInputs_);
	std::vector<std::vector<int>> idx_max(numberOfInputs_);
	
	if(do_min_){
		for(int i=0; i<numberOfInputs_; i++){
			std::vector<dataType> persistences;
			for(int j=0; j<bidder_diagrams_min_[i].size(); j++){
				Bidder<dataType> b = bidder_diagrams_min_[i].get(j);
				dataType persistence = b.getPersistence();
				if(persistence>=min_persistence && persistence<previous_min_persistence){
					candidates_to_be_added_min[i].push_back(j);
					idx_min[i].push_back(idx_min[i].size());
					persistences.push_back(persistence);
				}
			}
			sort(idx_min[i].begin(), idx_min[i].end(), [&persistences](int& a, int& b){
				return ((persistences[a] > persistences[b])
				||((persistences[a] == persistences[b])&&(a > b)));
				});
			int size =  candidates_to_be_added_min[i].size();
			if(size>=max_points_to_add){
				dataType last_persistence_added = persistences[idx_min[i][max_points_to_add-1]];
				if(last_persistence_added>new_min_persistence){
					new_min_persistence = last_persistence_added;
				}
			}
		}
	}
	
	if(do_sad_){
		for(int i=0; i<numberOfInputs_; i++){
			std::vector<dataType> persistences;
			for(int j=0; j<bidder_diagrams_saddle_[i].size(); j++){
				Bidder<dataType> b = bidder_diagrams_saddle_[i].get(j);
				dataType persistence = b.getPersistence();
				if(persistence>=min_persistence && persistence<previous_min_persistence){
					candidates_to_be_added_sad[i].push_back(j);
					idx_sad[i].push_back(idx_sad[i].size());
					persistences.push_back(persistence);
				}
			}
			sort(idx_sad[i].begin(), idx_sad[i].end(), [&persistences](int& a, int& b){
				return ((persistences[a] > persistences[b])
				||((persistences[a] == persistences[b])&&(a > b)));
				});
			int size =  candidates_to_be_added_sad[i].size();
			if(size>=max_points_to_add){
				dataType last_persistence_added = persistences[idx_sad[i][max_points_to_add-1]];
				if(last_persistence_added>new_min_persistence){
					new_min_persistence = last_persistence_added;
				}
			}
		}
	}
	if(do_max_){
		for(int i=0; i<numberOfInputs_; i++){
			std::vector<dataType> persistences;
			for(int j=0; j<bidder_diagrams_max_[i].size(); j++){
				Bidder<dataType> b = bidder_diagrams_max_[i].get(j);
				dataType persistence = b.getPersistence();
				if(persistence>=min_persistence && persistence<previous_min_persistence){
					candidates_to_be_added_max[i].push_back(j);
					idx_max[i].push_back(idx_max[i].size());
					persistences.push_back(persistence);
				}
			}
			sort(idx_max[i].begin(), idx_max[i].end(), [&persistences](int& a, int& b){
				return ((persistences[a] > persistences[b])
				||((persistences[a] == persistences[b])&&(a > b)));
				});
			int size =  candidates_to_be_added_max[i].size();
			if(size>=max_points_to_add){
				dataType last_persistence_added = persistences[idx_max[i][max_points_to_add-1]];
				if(last_persistence_added>new_min_persistence){
					new_min_persistence = last_persistence_added;
				}
			}
		}
	}
	
	// 3. Add the points to the current diagrams
	if(do_min_){
		for(int i=0; i<numberOfInputs_; i++){
			int size =  candidates_to_be_added_min[i].size();
			for(int j=0; j<std::min(max_points_to_add, size); j++){
				Bidder<dataType> b = bidder_diagrams_min_[i].get(candidates_to_be_added_min[i][idx_min[i][j]]);
				dataType persistence = b.getPersistence();
				if(persistence>=new_min_persistence){
					b.id_ = current_bidder_diagrams_min_[i].size();
					b.setPositionInAuction(current_bidder_diagrams_min_[i].size());
					b.setDiagonalPrice(initial_diagonal_prices[0][i]);
					current_bidder_diagrams_min_[i].addBidder(b);
					
					if(use_accelerated_ && n_iterations_>0){
						for(int c=0; c<k_; ++c){
							// Step 5 of Accelerated KMeans: Update the lower bound on distance thanks to the triangular inequality
							l_[i][c] = pow( pow(l_[i][c], 1./wasserstein_) - persistence/pow(2, 0.5), wasserstein_);
							if(l_[i][c]<0){
								l_[i][c] = 0;
							}
						}
						// Step 6, update the upper bound on the distance to the centroid thanks to the triangle inequality
						u_[i] = pow( pow(u_[i], 1./wasserstein_) + persistence/pow(2, 0.5), wasserstein_);
						r_[i] = true;
					}
				}
			}
			std::cout<< " Diagram " << i << " size : " << current_bidder_diagrams_min_[i].size() << std::endl;
		}
	}
	if(do_sad_){
		for(int i=0; i<numberOfInputs_; i++){
			int size =  candidates_to_be_added_sad[i].size();
			for(int j=0; j<std::min(max_points_to_add, size); j++){
				Bidder<dataType> b = bidder_diagrams_saddle_[i].get(candidates_to_be_added_sad[i][idx_sad[i][j]]);
				dataType persistence = b.getPersistence();
				if(persistence>=new_min_persistence){
					b.id_ = current_bidder_diagrams_saddle_[i].size();
					b.setPositionInAuction(current_bidder_diagrams_saddle_[i].size());
					b.setDiagonalPrice(initial_diagonal_prices[1][i]);
					current_bidder_diagrams_saddle_[i].addBidder(b);
					
					if(use_accelerated_ && n_iterations_>0){
						for(int c=0; c<k_; ++c){
							// Step 5 of Accelerated KMeans: Update the lower bound on distance thanks to the triangular inequality
							l_[i][c] = pow( pow(l_[i][c], 1./wasserstein_) - persistence/pow(2, 0.5), wasserstein_);
							if(l_[i][c]<0){
								l_[i][c] = 0;
							}
						}
						// Step 6, update the upper bound on the distance to the centroid thanks to the triangle inequality
						u_[i] = pow( pow(u_[i], 1./wasserstein_) + persistence/pow(2, 0.5), wasserstein_);
						r_[i] = true;
					}
				}
			}
			std::cout<< " Diagram " << i << " size : " << current_bidder_diagrams_saddle_[i].size() << std::endl;
		}
	}
	if(do_max_){
		for(int i=0; i<numberOfInputs_; i++){
			int size =  candidates_to_be_added_max[i].size();
			for(int j=0; j<std::min(max_points_to_add, size); j++){
				Bidder<dataType> b = bidder_diagrams_max_[i].get(candidates_to_be_added_max[i][idx_max[i][j]]);
				dataType persistence = b.getPersistence();
				if(persistence>=new_min_persistence){
					b.id_ = current_bidder_diagrams_max_[i].size();
					b.setPositionInAuction(current_bidder_diagrams_max_[i].size());
					b.setDiagonalPrice(initial_diagonal_prices[2][i]);
					current_bidder_diagrams_max_[i].addBidder(b);
					
					if(use_accelerated_ && n_iterations_>0){
						for(int c=0; c<k_; ++c){
							// Step 5 of Accelerated KMeans: Update the lower bound on distance thanks to the triangular inequality
							l_[i][c] = pow( pow(l_[i][c], 1./wasserstein_) - persistence/pow(2, 0.5), wasserstein_);
							if(l_[i][c]<0){
								l_[i][c] = 0;
							}
						}
						// Step 6, update the upper bound on the distance to the centroid thanks to the triangle inequality
						u_[i] = pow( pow(u_[i], 1./wasserstein_) + persistence/pow(2, 0.5), wasserstein_);
						r_[i] = true;
					}
				}
			}
			std::cout<< " Diagram " << i << " size : " << current_bidder_diagrams_max_[i].size() << std::endl;
		}
	}
	
	return new_min_persistence;
}


#endif
