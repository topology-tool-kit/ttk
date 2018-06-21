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
	std::cout<<"kdt built"<<std::endl;
	KDTree<dataType>* kdt = pair.first;
	std::cout<<kdt->id_<<std::endl;
	std::cout<<kdt->left_->id_<<std::endl;
	std::vector<KDTree<dataType>*>& correspondance_kdt_map = pair.second;
	

	for(int i=0; i<numberOfInputs_; i++){
		std::cout<<"Building Auction..."<<std::endl;
		Auction<dataType> auction = Auction<dataType>(bidder_diagrams_[i], barycenter_goods_[i], wasserstein_, geometrical_factor_, 0.01, kdt, correspondance_kdt_map);
		int n_biddings = 0;
		std::cout<<"Building unassigned bidders"<<std::endl;
		auction.buildUnassignedBidders();
		std::cout<<"Initializing goods"<<std::endl;
		auction.reinitializeGoods();
		std::cout<<"Running Auction..."<<std::endl;
		auction.runAuctionRound(n_biddings, i);
		std::cout<<"Auction round finished in " << n_biddings <<" biddings"<<std::endl;
		std::vector<matchingTuple> matchings;
		dataType cost = auction.getMatchingsAndDistance(&matchings);
		std::cout<<"cost = " << cost <<std::endl;
	}
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
	int random_idx = rand() % numberOfInputs_;
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
