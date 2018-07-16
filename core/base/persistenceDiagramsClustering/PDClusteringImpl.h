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
	//dataType cost = std::numeric_limits<dataType>::max();
	setBidderDiagrams();
	cost_ = std::numeric_limits<dataType>::max();
	dataType min_cost = std::numeric_limits<dataType>::max();
	epsilon_ = pow(getMostPersistent(), 2)/8.;
	dataType epsilon0 = epsilon_;
	
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
	
	bool converged = false;
	bool diagrams_complete= true;
	int n_iterations_ = 0;
	
	double total_time = 0;
	
	while(!converged || !diagrams_complete){
		Timer t;{
			std::cout<< "Epsilon = "<< epsilon_<< std::endl;
			
			n_iterations_++;
			dataType max_shift = updateCentroidsPosition();
			
			if(max_shift==0 || epsilon_<5e-5){
				converged = true;
			}
			epsilon_ = std::max(std::min(max_shift/8., 10.*epsilon0/n_iterations_), epsilon_/5.);
			/*if(use_progressive_){
				continue;
				//TODO Implement it ! (cf Python code)
			}*/
			
			if(use_accelerated_){
				acceleratedUpdateClusters();
			}
			else{
				updateClusters();
			}
			
			if(cost_<min_cost && n_iterations_>1){
				min_cost=cost_;
			}
			else if(n_iterations_>2){
				// TODO Adapt the condition with progressive KMeans
				converged = true;
			}
			std::cout<< "Cost = "<< cost_<< std::endl;
			printClustering();
		}
		total_time +=t.getElapsedTime();
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
			GoodDiagram<dataType> centroid_min = diagramToCentroid(bidder_diagrams_min_[idx[c]]);
			centroids_min_.push_back(centroid_min);
		}
		if(do_sad_){
			GoodDiagram<dataType> centroid_sad = diagramToCentroid(bidder_diagrams_saddle_[idx[c]]);
			centroids_saddle_.push_back(centroid_sad);
		}
		if(do_max_){
			GoodDiagram<dataType> centroid_max = diagramToCentroid(bidder_diagrams_max_[idx[c]]);
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
		GoodDiagram<dataType> centroid_min = diagramToCentroid(bidder_diagrams_min_[random_idx]);
		centroids_min_.push_back(centroid_min);
	}
	if(do_sad_){
		GoodDiagram<dataType> centroid_sad = diagramToCentroid(bidder_diagrams_saddle_[random_idx]);
		centroids_saddle_.push_back(centroid_sad);
	}
	if(do_max_){
		GoodDiagram<dataType> centroid_max = diagramToCentroid(bidder_diagrams_max_[random_idx]);
		centroids_max_.push_back(centroid_max);
	}
	
	while((int) indexes_clusters.size()<k_){
		std::vector<dataType> min_distance_to_centroid(numberOfInputs_);
		std::vector<dataType> probabilities(numberOfInputs_);
		dataType maximal_distance = 0;
		int candidate_centroid = -1;
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
						distance += computeDistance(bidder_diagrams_min_[i], centroid_min);
					}
					if(do_sad_){
						GoodDiagram<dataType> centroid_saddle = centroidWithZeroPrices(centroids_saddle_[j]);
						distance += computeDistance(bidder_diagrams_saddle_[i], centroid_saddle);
					}
					if(do_max_){
						GoodDiagram<dataType> centroid_max = centroidWithZeroPrices(centroids_max_[j]);
						distance += computeDistance(bidder_diagrams_max_[i], centroid_max);
					}
					if(distance<min_distance_to_centroid[i]){
						min_distance_to_centroid[i] = distance;
					}
				}
			}
			probabilities[i] = pow(min_distance_to_centroid[i], 2);
			
			// The following block is usefull in case of need for a deterministic algoritm
			if(min_distance_to_centroid[i]>maximal_distance){
				maximal_distance = min_distance_to_centroid[i];
				candidate_centroid = i;
			}
			
		}
		// Comment the following two lines to make it deterministic
		std::random_device rd;
		std::mt19937 gen(rd());
		std::discrete_distribution<int> distribution (probabilities.begin(),probabilities.end());
		candidate_centroid = distribution(gen);
		
		indexes_clusters.push_back(candidate_centroid);
		if(do_min_){
			GoodDiagram<dataType> centroid_min = diagramToCentroid(bidder_diagrams_min_[candidate_centroid]);
			centroids_min_.push_back(centroid_min);
		}
		if(do_sad_){
			GoodDiagram<dataType> centroid_sad = diagramToCentroid(bidder_diagrams_saddle_[candidate_centroid]);
			centroids_saddle_.push_back(centroid_sad);
		}
		if(do_max_){
			GoodDiagram<dataType> centroid_max = diagramToCentroid(bidder_diagrams_max_[candidate_centroid]);
			centroids_max_.push_back(centroid_max);
		}
	}
}

template <typename dataType>
void PDClustering<dataType>::initializeAcceleratedKMeans(){
	r_ = std::vector<bool>(numberOfInputs_);
	u_ = std::vector<dataType>(numberOfInputs_);
	inv_clustering_ = std::vector<int>(numberOfInputs_);
	for(int i=0; i<numberOfInputs_; i++){
		r_[i] = true;
		u_[i] = std::numeric_limits<dataType>::max();
		inv_clustering_[i] = -1;
	}
	
	l_ = std::vector<std::vector<dataType>>(numberOfInputs_);
	for(int i=0; i<numberOfInputs_; ++i){
		l_[i] = std::vector<dataType>(k_);
		for(int c=0; c<k_; ++c){
			l_[i][c]=0;
		}
	}
	
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
			D1_min = diagramWithZeroPrices(bidder_diagrams_min_[i]);
		}
		if(do_sad_){
			D1_sad = diagramWithZeroPrices(bidder_diagrams_saddle_[i]);
		}
		if(do_max_){
			D1_max = diagramWithZeroPrices(bidder_diagrams_max_[i]);
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
				distance += computeDistance(D1_sad, D2_sad);
			}
			if(do_max_){
				D2_max = centroidWithZeroPrices(centroids_max_[j]);
				distance += computeDistance(D1_max, D2_max);
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
			std::cout<< " Problem in invertInverseClusters()... \n Cluster " << c << " was left with no diagram attached to it... " << std::endl;
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
			D1_min = diagramWithZeroPrices(bidder_diagrams_min_[i]);
		}
		if(do_sad_){
			D1_sad = diagramWithZeroPrices(bidder_diagrams_saddle_[i]);
		}
		if(do_max_){
			D1_max = diagramWithZeroPrices(bidder_diagrams_max_[i]);
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
				}
				// Step 3b, check if still potential change of clusters
				if(u_[i]>l_[i][c] || u_[i]>0.5*d_[inv_clustering_[i]][c]){
					BidderDiagram<dataType> diagram_min, diagram_sad, diagram_max;
					GoodDiagram<dataType> centroid_min, centroid_sad, centroid_max;
					dataType distance = 0;
					
					if(do_min_){
						centroid_min = centroidWithZeroPrices(centroids_min_[c]);
						diagram_min = diagramWithZeroPrices(bidder_diagrams_min_[i]);
						distance += computeDistance(diagram_min, centroid_min, 0.01);
					}
					if(do_sad_){
						centroid_sad = centroidWithZeroPrices(centroids_saddle_[c]);
						diagram_sad = diagramWithZeroPrices(bidder_diagrams_saddle_[i]);
						distance += computeDistance(diagram_sad, centroid_sad, 0.01);
					}
					if(do_max_){
						centroid_max = centroidWithZeroPrices(centroids_max_[c]);
						diagram_max = diagramWithZeroPrices(bidder_diagrams_max_[i]);
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
			int idx = rand() % k_;
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
	cost_ = 0;
	for(int c=0; c<k_; ++c){
		std::vector<GoodDiagram<dataType> > centroids_with_price_min, centroids_with_price_sad, centroids_with_price_max;
		int count = 0;
		for(int idx : clustering_[c]){
			int number_of_points = 0;
			// Find the position of diagrams[idx] in old cluster c
			vector<int>::iterator i = find (old_clustering_[c].begin(), old_clustering_[c].end (), idx);
			int pos = (i!=old_clustering_[c].end()) ? -1 : distance (old_clustering_[c].begin(), i);
			if(pos>=0){
				// Diagram was already linked to this centroid before
				if(do_min_){
					centroids_with_price_min.push_back(centroids_with_price_min_[idx]);
					number_of_points += centroids_with_price_min_[idx].size() + bidder_diagrams_min_[idx].size();
				}
				if(do_sad_){
					centroids_with_price_sad.push_back(centroids_with_price_saddle_[idx]);
					number_of_points += centroids_with_price_saddle_[idx].size() + bidder_diagrams_saddle_[idx].size();
				}
				if(do_max_){
					centroids_with_price_max.push_back(centroids_with_price_max_[idx]);
					number_of_points += centroids_with_price_max_[idx].size() + bidder_diagrams_max_[idx].size();
				}
			}
			else{
				// Otherwise, centroid is given 0 prices and the diagram is given 0 diagonal-prices
				if(do_min_){
					centroids_with_price_min.push_back(centroidWithZeroPrices( centroids_min_[c] ));
					bidder_diagrams_min_[idx] = diagramWithZeroPrices(bidder_diagrams_min_[idx]);
					number_of_points += centroids_with_price_min_[idx].size() + bidder_diagrams_min_[idx].size();
				}
				if(do_sad_){
					centroids_with_price_sad.push_back(centroidWithZeroPrices( centroids_saddle_[c] ));
					bidder_diagrams_saddle_[idx] = diagramWithZeroPrices(bidder_diagrams_saddle_[idx]);
					number_of_points += centroids_with_price_saddle_[idx].size() + bidder_diagrams_saddle_[idx].size();
				}
				if(do_max_){
					centroids_with_price_max.push_back(centroidWithZeroPrices( centroids_max_[c] ));
					bidder_diagrams_max_[idx] = diagramWithZeroPrices(bidder_diagrams_max_[idx]);
					number_of_points += centroids_with_price_max_[idx].size() + bidder_diagrams_max_[idx].size();
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
						computeDistance(&(bidder_diagrams_min_[idx]), &(centroids_with_price_min[count]), estimated_delta_lim);
					}
					if(do_sad_){
						computeDistance(&(bidder_diagrams_saddle_[idx]), &(centroids_with_price_sad[count]), estimated_delta_lim);
					}
					if(do_max_){
						computeDistance(&(bidder_diagrams_max_[idx]), &(centroids_with_price_max[count]), estimated_delta_lim);
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
				diagrams_c_min.push_back(bidder_diagrams_min_[idx]);
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
			std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> pair = barycenter_computer.getKDTree();
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
										&all_matchings);
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
				bidder_diagrams_min_[idx] = diagrams_c_min[i];
				centroids_with_price_min_[idx] = centroids_with_price_min[i];
				i++;
			}
			
			GoodDiagram<dataType> old_centroid = centroids_min_[c];
			centroids_min_[c] = centroidWithZeroPrices(centroids_with_price_min_[clustering_[c][0]]);
			wasserstein_shift += computeDistance(old_centroid, centroids_min_[c], 0.01);
		}
		
		if(do_sad_){
			
			for(int idx : clustering_[c]){
				diagrams_c_sad.push_back(bidder_diagrams_saddle_[idx]);
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
			std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> pair = barycenter_computer.getKDTree();
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
										&all_matchings);
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
				bidder_diagrams_saddle_[idx] = diagrams_c_sad[i];
				centroids_with_price_saddle_[idx] = centroids_with_price_sad[i];
				i++;
			}
			GoodDiagram<dataType> old_centroid = centroids_saddle_[c];
			centroids_saddle_[c] = centroidWithZeroPrices(centroids_with_price_saddle_[clustering_[c][0]]);
			wasserstein_shift += computeDistance(old_centroid, centroids_saddle_[c], 0.01);
		}
		
		if(do_max_){
			
			for(int idx : clustering_[c]){
				diagrams_c_max.push_back(bidder_diagrams_max_[idx]);
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
			std::pair<KDTree<dataType>*, std::vector<KDTree<dataType>*>> pair = barycenter_computer.getKDTree();
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
										&all_matchings);
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
				bidder_diagrams_max_[idx] = diagrams_c_max[i];
				centroids_with_price_max_[idx] = centroids_with_price_max[i];
				i++;
			}
			GoodDiagram<dataType> old_centroid = centroids_max_[c];
			centroids_max_[c] = centroidWithZeroPrices(centroids_with_price_max_[clustering_[c][0]]);
			wasserstein_shift += computeDistance(old_centroid, centroids_max_[c], 0.01);
		}
		cost_+=total_cost;
		
		if(use_accelerated_){
			for(int idx : clustering_[c]){
				// Step 5 of Accelerated KMeans: Update the lower bound on distance thanks to the triangular inequality
				l_[idx][c] = pow( pow(l_[idx][c], 1./wasserstein_) - pow(wasserstein_shift, 1./wasserstein_), wasserstein_);
				if(l_[idx][c]<0){
					l_[idx][c] = 0;
				}
				// Step 6, update the upper bound on the distance to the centroid thanks to the triangle inequality
				u_[idx] = pow( pow(u_[idx], 1./wasserstein_) + pow(wasserstein_shift, 1./wasserstein_), wasserstein_);
				r_[idx] = true;
			}
		}
	}
	return max_shift;
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

#endif
