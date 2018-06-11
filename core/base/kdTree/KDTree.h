/// \ingroup base
/// \class ttk::KDTree 
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date June 2018
///
/// \brief TTK KD-Tree
///

#pragma once

#ifndef _KDTREE_H
#define _KDTREE_H

// base code includes
#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm>
#include <Debug.h>



namespace ttk{
  template<typename dataType>
  class KDTree : public Debug
  {
	
	public:
		
		KDTree(bool include_weights, int p){
			left_ = nullptr;
			right_ = nullptr;
			parent_ = nullptr;
			
			coords_number_ = 0;
			p_ = p;
			include_weights_ = include_weights;
		}
		
		
		
		KDTree(KDTree* father, int coords_number, bool is_left){
			left_ = nullptr;
			right_ = nullptr;
			parent_ = father;
			
			coords_number_ = coords_number;
			p_ = father->p_;
			include_weights_ = father->include_weights_;
		}
		
		~KDTree(){
			delete left_;
			delete right_;
		}
		
		void build(dataType* coordinates, const int& ptNumber, const int& dimension);
		void buildRecursive(dataType* coordinates, std::vector<int>& indexes, const int& ptNumber, const int& dimension, KDTree<dataType>* parent);
		void updateWeight(dataType new_weight);
		void updateMinSubweight(); 
		std::vector<int> getKClosest(int k, std::vector<dataType>& coordinates);
		
		bool isLeaf();
		bool isRoot();
		
	protected:
		KDTree* left_;		   // Lower half for the coordinate specified    
		KDTree* right_;		   // Higher half
		KDTree* parent_;
		int id_;			   // ID of the object saved here. The whole object is not kept in the KDTree
							   // Users should keep track of them in a table for instance
		bool is_left_;		   // Boolean indicating if the current node is a left node of its parent
		int coords_number_;    // Indicates according to which coordinate the tree splits its elements
		
		int p_;                // Power used for the computation of distances. p=2 yields euclidean distance
		
		bool include_weights_;   // Wether or not the KDTree should include weights that add up 
							   // to distance for the computation of nearest neighbours
		
		dataType weight_;
		dataType min_subweights_;
		std::vector<dataType> coordinates_;
		
	};
  
	template<typename dataType>
	void KDTree<dataType>::build(dataType* data, const int& ptNumber, const int& dimension){
		
		// First, perform a argsort on the data
		// initialize original index locations
		std::vector<int> idx(ptNumber * dimension);
		for(int i=0; i<ptNumber*dimension; i++){
			idx[i] = i;
		}
		// sort indexes based on comparing values in coordinates
		sort(idx.begin(), idx.end(), [&](int i1, int i2) {return *(data + dimension*i1+coords_number_) < *(data + dimension*i2+coords_number_);});
		
		int median_idx = idx[(int) (ptNumber-1)/2];
		for(int axis=0; axis<dimension; axis++){
			coordinates_.push_back(*(data + dimension*median_idx + axis) ); 
		}
		weight_ = 0;
		min_subweights_ = 0;
		id_ = median_idx;
		parent_ = nullptr;
		
		if(idx.size()>1){
			// Build left leaf
			int idx_left[median_idx];
			for(int i=0; i<median_idx; i++){
				idx_left[i] = idx[i];
			}
			
			KDTree left = KDTree(this, (coords_number_+1)%dimension, true);
			left.buildRecursive(data, idx_left, ptNumber, dimension, this);
			left_ = &left;
		}
		
		if(idx.size()>0){
			// Build right leaf
			int idx_right[ptNumber - median_idx - 1];
			for(int i=0; i<ptNumber - median_idx - 1; i++){
				idx_right[i] = idx[i + median_idx + 1];
			}
			KDTree right = KDTree(this, (coords_number_+1)%dimension, false);
			right.buildRecursive(data, idx_right, ptNumber, dimension, this);
			right_ = &right;
		}

		return ;
	}
  
	template<typename dataType>
	void KDTree<dataType>::buildRecursive(dataType* data, std::vector<int>& idx, const int& ptNumber, const int& dimension, KDTree<dataType>* parent){
		
		// First, perform a argsort on the data
		sort(idx.begin(), idx.end(), [&](int i1, int i2) {return *(data + dimension*idx[i1] + coords_number_) < *(data + dimension*idx[i2] + coords_number_);});
		int median_idx = idx[(int) (idx.size()+1)/2];
		
		for(int axis=0; axis<dimension; axis++){
			coordinates_.push_back(*(data + dimension*median_idx + axis) ); 
		}
		weight_ = 0;
		min_subweights_ = 0;
		id_ = median_idx;
		parent_ = parent;
		
		if(idx.size()>1){
			// Build left leaf
			int idx_left[median_idx];
			for(int i=0; i<median_idx; i++){
				idx_left[i] = idx[i];
			}
			
			KDTree left = KDTree(this, (coords_number_+1)%dimension, true);
			left.buildRecursive(data, idx_left, ptNumber, dimension, this);
			left_ = &left;
		}
		
		if(idx.size()>0){
			// Build right leaf
			int idx_right[ptNumber - median_idx - 1];
			for(int i=0; i<ptNumber - median_idx - 1; i++){
				idx_right[i] = idx[i + median_idx + 1];
			}
			KDTree right = KDTree(this, (coords_number_+1)%dimension, false);
			right.buildRecursive(data, idx_right, ptNumber, dimension, this);
			right_ = &right;
		}
		return;
	}
	
	template<typename dataType>
	void KDTree< dataType>::updateWeight(dataType new_weight){
		// TODO
		return;
	}
	
	template<typename dataType>
	void KDTree<dataType>::updateMinSubweight(){
		// TODO
		return;
	}
	
	
	template<typename dataType>
	std::vector<int> KDTree<dataType>::getKClosest(int k, std::vector<dataType>& coordinates){
		//TODO
		std::vector<int> neighbours;
		return neighbours;
	}
	
	
	
	template<typename dataType>
	bool KDTree<dataType>::isLeaf(){
		return left_==nullptr && right_==nullptr;
	}
	
	template<typename dataType>
	bool KDTree<dataType>::isRoot(){
		return parent_==nullptr;
	}
}

#endif
