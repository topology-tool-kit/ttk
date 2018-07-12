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
	
	
	} // End of timer
	return 0;
}

template <typename dataType>
dataType PDClustering<dataType>::getMostPersistent(int id_of_diagram, int diagram_type){
	//TODO
	return 0;
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
	//TODO
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
			if(min_distance_to_centroid[i]>maximal_distance){
				// TODO Make it probabilistic !
				maximal_distance = min_distance_to_centroid[i];
				candidate_centroid = i;
			}
		}
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
	
	// Check if a diagram was left without cluster
	for(int i=0; i<numberOfInputs_; ++i){
		if(inv_clustering_[i] == -1){
			std::cout<< " Problem in invertClusters()... \n Diagram " << i << " was left with no cluster attached to it... " << std::endl;
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
	//TODO
	// Step 1
	getCentroidDistanceMatrix();
	
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
		}
	}
	return ;
}

template <typename dataType>
void PDClustering<dataType>::updateCentroidsPosition(){
	//TODO
	
	
	
	
	
	return ;
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
		}
	}
	return;
}

#endif
