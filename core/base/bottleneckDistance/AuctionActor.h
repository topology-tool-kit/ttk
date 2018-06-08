/// \ingroup base
/// \class ttk::AuctioActor 
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date 4th June 2018

#ifndef _AUCTIONACTOR_H
#define _AUCTIONACTOR_H

#include <cmath>
#include <limits>
#include <iostream>
#include <Debug.h>
#include <PersistenceDiagram.h>


namespace ttk{
  template<typename dataType> class Good;
  template<typename dataType> class Bidder;
  template<typename dataType> class GoodDiagram;
  template<typename dataType> class BidderDiagram;
	
	
 template<typename dataType>
  class AuctionActor : public Debug
  {
	public:
		dataType x_;
		dataType y_;
		int id_;
		
		AuctionActor() {
			id_ = 0;
			x_= 0;
			y_= 0;
			is_diagonal_ = false;
		}
		
		AuctionActor(dataType x, dataType y, bool is_diagonal, int id) {
			id_ = id;
			x_ = x;
			y_ = y;
			is_diagonal_ = is_diagonal;
		}
		~AuctionActor() {};
		
		
		void SetCoordinates(dataType x, dataType y);
		void projectOnDiagonal();
		int getId();
		dataType getPersistence();
		bool isDiagonal();
		
		template<typename type>
		inline static type abs(const type var) {
			return (var >= 0) ? var : -var;
		}
		
		dataType cost(AuctionActor& g, int wasserstein);
		
		inline dataType cost(AuctionActor* g, int wasserstein);
		
	protected:
		bool is_diagonal_;

      
  };
  
	template<typename dataType>
	int AuctionActor<dataType>::getId(){
		return id_;
	}
  
	template<typename dataType>
	void AuctionActor<dataType>::SetCoordinates(dataType x, dataType y){
		x_ = x;
		y_ = y;
	}
	
	template<typename dataType>
	void AuctionActor<dataType>::projectOnDiagonal(){
		x_ = (x_+y_)/2;
		y_ = x_;
		is_diagonal_ = true;
	}

	template<typename dataType>
	dataType AuctionActor<dataType>::getPersistence(){
		return y_-x_;
	}
	
	
	template<typename dataType>
	bool AuctionActor<dataType>::isDiagonal(){
		return (is_diagonal_ || (x_==y_));
	}
	
	template<typename dataType>
	dataType AuctionActor<dataType>::cost(AuctionActor& g, int wasserstein){
		if(is_diagonal_ && g.isDiagonal()){
			return 0;
		}
		else{
			return pow(abs<dataType>(x_ - g.x_) , wasserstein) + pow(abs<dataType>(y_ - g.y_), wasserstein);
		}
	}
	
	template<typename dataType>
	dataType AuctionActor<dataType>::cost(AuctionActor* g, int wasserstein){
		return this->cost(*g, wasserstein);
	}

	
  template<typename dataType>
  class Good : public AuctionActor<dataType>{
   
     public:
         Good() {};
         Good(dataType x, dataType y, bool is_diagonal, int id)  : AuctionActor<dataType>(x, y, is_diagonal, id)
		 {
			 owner_= -1;
			 price_ = 0;
         }
         
         Good(diagramTuple tuple, int id) : AuctionActor<dataType>(){
			AuctionActor<dataType>::id_ = id;
			
			dataType x = std::get<6>(tuple);
			dataType y = std::get<10>(tuple);
			this->SetCoordinates(x, y);
			if(x==y){AuctionActor<dataType>::is_diagonal_=true;}
			else{AuctionActor<dataType>::is_diagonal_=false;}
            price_ = 0;
			owner_ = -1;
         }
         
         ~Good() {};
         
         inline void setPrice(dataType price){
             price_ = price;
         }
             
		inline dataType getPrice(){
			return price_;
		}
		
		inline void assign(int b, dataType price){
			price_ = price;
			owner_ = b;
		}
		
		inline int getOwner(){
			return owner_;
		}
		
		inline void setOwner(int idx){
			owner_ = idx;
		}
         
     protected:
         dataType price_;
		 // Position in Auction.bidders_ of the owner of this good.
		 // If the good is not owned, owner_ is set to -1
		 int owner_;
      
  };
  
	
	template<typename dataType>
	class GoodDiagram{

		public:
			GoodDiagram() {};
			~GoodDiagram() {};
			
			void addGood(Good<dataType> g);
			Good<dataType>& get(int idx);
			void set(Good<dataType> g, int idx); 
			
			inline int size(){
				return goods_.size(); 
			}
			
		private:
			std::vector<Good<dataType>> goods_;
	};


	template<typename dataType>
	void GoodDiagram<dataType>::addGood(Good<dataType> g){
		goods_.push_back(g);
	}

	template<typename dataType>
	Good<dataType>& GoodDiagram<dataType>::get(int idx){
		return goods_[idx];
	}  
	
	template<typename dataType>
	void GoodDiagram<dataType>::set(Good<dataType> g, int idx){
		goods_[idx] = g;
	}  
  
  
  
  template<typename dataType>
  class Bidder : public AuctionActor<dataType>
  {
     public:
         Bidder() {
		};
         Bidder(dataType x, dataType y, bool is_diagonal, int id) : AuctionActor<dataType>(x, y, is_diagonal, id)
		 {
             price_paid_ = 0;
             diagonal_price_ = 0;
			 property_ = NULL;
         }
         
		Bidder(diagramTuple tuple, int id) : AuctionActor<dataType>()  {
			AuctionActor<dataType>::id_ = id;
			dataType x = std::get<6>(tuple);
			dataType y = std::get<10>(tuple);
			this->SetCoordinates(x, y);
			if(x==y){
				AuctionActor<dataType>::is_diagonal_=true;
			}
			else{
				AuctionActor<dataType>::is_diagonal_=false;
			}
            price_paid_ = 0;
            diagonal_price_ = 0;
			property_ = NULL;
         }
         
		~Bidder() {}

		Good<dataType>* getProperty();
		int runBidding(GoodDiagram<dataType>& goods, Good<dataType>& diagonalGood, int wasserstein, dataType epsilon);
		void setDiagonalPrice(dataType price);
		void setPricePaid(dataType price);
        void setProperty(Good<dataType>* g);
		void setTwin(Good<dataType>* g);
		void setPositionInAuction(int pos);
		int getPositionInAuction();
		
     protected:
		bool is_diagonal_;
		dataType price_paid_;
		dataType diagonal_price_;
		Good<dataType>* property_;
		
	private:
		// Attribute stating at which position in Auction.bidders_ this bidder can be found
		// In a single Auction, this attribute could be deducted from id_, but with the use of Barycenter
		// Goods will be added and deleted, which would add and delete diagonal bidders.
		int position_in_auction_;
      
  };
  
	template<typename dataType>
	Good<dataType>* Bidder<dataType>::getProperty(){
		return property_;
	}
	
	
	template<typename dataType>
	int Bidder<dataType>::runBidding(GoodDiagram<dataType>& goods, Good<dataType>& diagonalGood, int wasserstein, dataType epsilon){
		//TODO Adjust for goodDiagrams with only one point...
		dataType best_val = std::numeric_limits<dataType>::lowest();   //Higher than sqrt(2)*100 times the closest good
		dataType second_val = std::numeric_limits<dataType>::lowest();
		Good<dataType>* best_good = nullptr;
		for(int i=0; i<goods.size(); i++){
			Good<dataType>& g = goods.get(i);
			dataType val = -this->cost(g, wasserstein);
			val -= g.getPrice();
			if(val>best_val){
				second_val = best_val;
				best_val = val;
				best_good = &g;
			}
			else if(val>second_val){
				second_val=val;
			}
		}
		// And now check for the corresponding twin bidder
		Good<dataType>& g = diagonalGood;
		dataType val = -this->cost(g, wasserstein);
		val -= g.getPrice();
		if(val>best_val){
			second_val = best_val;
			best_val = val;
			best_good = &g;
		}
		else if(val>second_val){
			second_val=val;
		}

		dataType old_price = best_good->getPrice();
		dataType new_price = old_price + best_val-second_val + epsilon;
		// Assign bidder to best_good
		this->setProperty(best_good);
		this->setPricePaid(new_price);
		
		// Assign best_good to bidder and unassign the previous owner of best_good if need be
		int idx_reassigned = best_good->getOwner();
		best_good->assign(this->position_in_auction_, new_price);
		return idx_reassigned;
	}

	template<typename dataType>
	void Bidder<dataType>::setDiagonalPrice(dataType price){
		diagonal_price_ = price;
	}
	
	template<typename dataType>
	void Bidder<dataType>::setPricePaid(dataType price){
		price_paid_ = price;
	}

	
	template<typename dataType>
	void Bidder<dataType>::setProperty(Good<dataType>* g){
		property_ = g;
	}
	
	template<typename dataType>
	void Bidder<dataType>::setPositionInAuction(int pos){
		position_in_auction_ = pos;
	}
		
	
	template<typename dataType>
	int Bidder<dataType>::getPositionInAuction(){
		return position_in_auction_;
	}


  
    
  
  
  template<typename dataType>
  class BidderDiagram{
   
  public:
	BidderDiagram() {};
	~BidderDiagram() {};
      
	inline void addBidder(Bidder<dataType> b){
		bidders_.push_back(b);
	}
	  
	inline int size(){
		 return bidders_.size(); 
	}
	  
	inline Bidder<dataType>& get(int idx){
		return bidders_[idx];
	}  
	  
  private:
      std::vector<Bidder<dataType>> bidders_;
  };



}


#endif
