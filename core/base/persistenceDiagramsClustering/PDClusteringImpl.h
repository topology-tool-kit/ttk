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
int PDClustering<dataType>::execute(){
	return 0;
}

template <typename dataType>
dataType PDClustering<dataType>::getMostPersistent(int id_of_diagram, int diagram_type){
	//TODO
	return 0;
}

template <typename dataType>
dataType PDClustering<dataType>::computeDistance(BidderDiagram<dataType> D1, BidderDiagram<dataType> D2, dataType delta_lim=0.0001){
	//TODO
	return 0;
}

template <typename dataType>
dataType PDClustering<dataType>::computeDistance(BidderDiagram<dataType> D1, GoodDiagram<dataType> D2, dataType delta_lim=0.0001){
	//TODO
	return 0;
}

template <typename dataType>
dataType PDClustering<dataType>::computeDistance(GoodDiagram<dataType> D1, GoodDiagram<dataType> D2, dataType delta_lim=0.0001){
	//TODO
	return 0;
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
BidderDiagram<dataType> PDClustering<dataType>::centroidToDiagram(GoodDiagram<dataType> centroid){
	BidderDiagram<dataType> BD = BidderDiagram<dataType>();
	for(int i=0; i<centroid.size(); i++){
		Good<dataType> g = centroid.get(i);
		
		Bidder<dataType> b = Bidder<dataType>(g.x_, g.y_, g.is_diagonal_, BD.size());
		b.SetCriticalCoordinates(g.coords_x_, g.coords_y_, g.coords_z_);
		BD.addBidder(b);
	}
	return BD;
}

template <typename dataType>
GoodDiagram<dataType> PDClustering<dataType>::diagramToCentroid(BidderDiagram<dataType> diagram){
	GoodDiagram<dataType> GD = GoodDiagram<dataType>();
	for(int i=0; i<diagram.size(); i++){
		Bidder<dataType> b = diagram.get(i);
		
		Good<dataType> g = Good<dataType>(b.x_, b.y_, b.is_diagonal_, BD.size());
		g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
		BD.addGood(g);
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
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(idx.begin(), idx.end(), g);
	
	for(int c=0; c<k_; c++){
		if(do_min_;){
			GoodDiagram<dataType> centroid_min = diagramToCentroid(bidder_diagrams_min_[idx[c]]);
			centroids_min_.push_back(centroid_min);
		}
		if(do_sad_){
			GoodDiagram<dataType> centroid_sad = diagramToCentroid(bidder_diagrams_sad_[idx[c]]);
			centroids_sad_.push_back(centroid_sad);
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
	return ;
}

template <typename dataType>
void PDClustering<dataType>::initializeAcceleratedKMeans(){
	//TODO
	return ;
}


template <typename dataType>
std::vector<std::vector<dataType>> PDClustering<dataType>::getDistanceMatrix(){
	//TODO
	std::vector<std::vector<dataType>> D;
	return D;
}

template <typename dataType>
std::vector<std::vector<dataType>> PDClustering<dataType>::getCentroidDistanceMatrix(){
	//TODO
	std::vector<std::vector<dataType>> D;
	return D;
}


template <typename dataType>
void PDClustering<dataType>::updateClusters(){
	//TODO
	return ;
}

template <typename dataType>
void PDClustering<dataType>::invertClusters(){
	//TODO
	return ;
}

template <typename dataType>
void PDClustering<dataType>::invertInverseClusters(){
	//TODO
	return ;
}


template <typename dataType>
void PDClustering<dataType>::acceleratedUpdateClusters(){
	//TODO
	return ;
}

template <typename dataType>
void PDClustering<dataType>::updateCentroidsPosition(){
	//TODO
	return ;
}



#endif
