/// \ingroup base
/// \class ttk::PersistenceDiagramsBarycenter
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \date August 2016.
///
/// \brief TTK processing package that takes an input ensemble data set 
/// (represented by a list of scalar fields) and which computes various 
/// vertexwise statistics (PDF estimation, bounds, moments, etc.)
///
/// \sa ttkPersistenceDiagramsBarycenter.cpp %for a usage example.

#ifndef _PERSISTENCEDIAGRAMSBARYCENTER_H
#define _PERSISTENCEDIAGRAMSBARYCENTER_H

// base code includes
#include                  <Wrapper.h>


namespace ttk{

  template <class dataType>
  class PDFBounds : public Debug {
  public:

    PDFBounds() {
      numberOfVertices_ = 0;
    }

    ~PDFBounds() {
      flush();
    }

    int evaluateRealization(const void *voidPointer) {
      #ifdef TTK_ENABLE_KAMIKAZE
      if(!(numberOfVertices_>0)) {
        return -1; // Number of vertices not defined
      }
      #endif
      const dataType *inputData =
        reinterpret_cast<const dataType*>(voidPointer);
      unsigned int numberOfVertices =
        static_cast<unsigned int>(numberOfVertices_);
      /* Initialize if first call since a change */
      if (!(upperBound_.size()==numberOfVertices) || !(lowerBound_.size()==numberOfVertices)) {
        upperBound_.resize(numberOfVertices);
        lowerBound_.resize(numberOfVertices);
        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for num_threads(threadNumber_)
        #endif
        for (size_t i = 0 ; i < numberOfVertices ; i++) {
          upperBound_[i] = inputData[i];
          lowerBound_[i] = inputData[i];
        }
      } else { /* Update the two fields with the new input */
        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for num_threads(threadNumber_)
        #endif
        for (size_t i = 0 ; i < numberOfVertices ; i++) {
          // Upper Bound
          if (inputData[i] > upperBound_[i]) {
            upperBound_[i] = inputData[i];
          }
          // Lower Bound
          if (inputData[i] < lowerBound_[i]) {
            lowerBound_[i] = inputData[i];
          }
        }
      }
      return 0;
    }

    std::pair<dataType,dataType> getRange() const {
      std::pair<dataType,dataType> range;
      range.first = getRangeMin();
      range.second = getRangeMax();
      return range;
    }

    dataType getRangeMax() const {
      if(upperBound_.size()) {
        dataType maxValue = upperBound_[0];
        for (size_t i = 1; i < upperBound_.size(); i++) {
          if (upperBound_[i] > maxValue) {
            maxValue = upperBound_[i];
          }
        }
        return maxValue;
      } else {
        return 0;
      }
    }

    dataType getRangeMin() const {
      if(lowerBound_.size()) {
        dataType minValue = lowerBound_[0];
        for (size_t i = 1; i < lowerBound_.size(); i++) {
          if (lowerBound_[i] < minValue) {
            minValue = lowerBound_[i];
          }
        }
        return minValue;
      } else {
        return 0;
      }
    }

    inline int flush() {
      numberOfVertices_ = 0;
      upperBound_.clear();
      lowerBound_.clear();
      return 0;
    }

    inline dataType* getLowerBoundPointer() {
      return lowerBound_.data();
    }

    inline dataType* getUpperBoundPointer() {
      return upperBound_.data();
    }

    inline int setNumberOfVertices(const int number) {
      numberOfVertices_ = number;
      return 0;
    }

  protected:
    int numberOfVertices_;
    std::vector<dataType> upperBound_;
    std::vector<dataType> lowerBound_;
  };

  class PDFHistograms : public Debug {
  public:

    PDFHistograms() {
      numberOfBins_ = 0;
      numberOfInputs_ = 0;
      numberOfVertices_ = 0;
      rangeMax_ = 0;
      rangeMin_ = 0;
    }

    ~PDFHistograms() {
      flush();
    }

    template <class dataType>
    int evaluateRealization(const dataType *inputData) {
      #ifdef TTK_ENABLE_KAMIKAZE
      if(!(rangeMin_<rangeMax_)) {
        return -1; // Range error
      }
      if(!(numberOfBins_>0)) {
        return -2; // Number of bins not defined
      }
      if(!(numberOfVertices_>0)) {
        return -3; // Number of vertices not defined
      }
      #endif
      if(numberOfInputs_ == 0) {
        /* Initialize */
        probability_.resize(numberOfBins_);
        double dx = (rangeMax_-rangeMin_) / static_cast<double>(numberOfBins_);
        for (size_t i=0 ; i < numberOfBins_ ; i++) {
          probability_[i].resize(numberOfVertices_);
          binValue_[i] = rangeMin_+(dx/2.0)+(static_cast<double>(i)*dx);
        }
      }
      /* Add input datas */
      for(unsigned int i=0 ; i<numberOfVertices_ ; i++) {
        unsigned int bin = static_cast<unsigned int>(floor((inputData[i]-rangeMin_)*numberOfBins_/(rangeMax_-rangeMin_)));
        bin = (bin == numberOfBins_) ? numberOfBins_-1 : bin;
        probability_[bin][i] += 1.0;
      }
      numberOfInputs_++;
      return 0;
    }

    inline int flush() {
      binValue_.clear();
      numberOfBins_ = 0;
      numberOfInputs_ = 0;
      numberOfVertices_ = 0;
      rangeMax_ = 0;
      rangeMin_ = 0;
      // selection_.clear(); // TODO : selection support
      return 0;
    }

    inline double* getBinFieldPointer(const unsigned int binId) {
      if(binId < numberOfBins_) {
        return probability_[binId].data();
      } else {
        return NULL;
      }
    }

    void getVertexHistogram(const unsigned int vertexId, std::vector<double> &histogram) const {
      histogram.resize(numberOfBins_);
      if(vertexId < numberOfVertices_) {
        #ifdef TTK_ENABLE_OPENMP
        #ifdef _WIN32
        #pragma omp parallel for num_threads(threadNumber_)
        #else
        #pragma omp parallel for num_threads(threadNumber_) \
          schedule(static, numberOfBins_/threadNumber_)
        #endif
        #endif
        for(int i=0 ; i< (int) numberOfBins_ ; i++) {
          if(probability_[i].size()==numberOfVertices_) {
            histogram[i] = probability_[i][vertexId];
          } else {
            histogram[i] = 0.0;
          }
        }
      } else {
        fill(histogram.begin(), histogram.end(), 0.0);
      }
    }

    int normalize() {
      const double normalization = 1.0 / static_cast<double>(numberOfInputs_);
      #ifdef TTK_ENABLE_OPENMP
      #ifdef _WIN32
      #pragma omp parallel for num_threads(threadNumber_)
      #else
      #pragma omp parallel for num_threads(threadNumber_) collapse(2) \
        schedule(static, (numberOfBins_*numberOfVertices_)/threadNumber_)
      #endif
      #endif
      for(int i=0 ; i< (int) numberOfBins_ ; i++) {
        for(int j=0 ; j< (int) numberOfVertices_ ; j++) {
          probability_[i][j] *= normalization;
        }
      }
      return 0;
    }

    inline int setNumberOfBins(const unsigned int number) {
      numberOfBins_ = number;
      return 0;
    }

    inline int setNumberOfVertices(const unsigned int number) {
      numberOfVertices_ = number;
      return 0;
    }

    inline int setRange(const double min, const double max) {
      rangeMin_ = min;
      rangeMax_ = max;
      return 0;
    }

  protected:
    std::vector<double>            binValue_;
    std::vector<std::vector<double> >   probability_;
    unsigned int              numberOfBins_;
    unsigned int              numberOfInputs_;
    unsigned int              numberOfVertices_;
    double                    rangeMin_;
    double                    rangeMax_;
    // std::vector<int> selection_; // TODO : selection support
  };

  class PersistenceDiagramsBarycenter : public Debug{

    public:

      PersistenceDiagramsBarycenter();

      ~PersistenceDiagramsBarycenter();

      /// Execute the package.
      /// \return Returns 0 upon success, negative values otherwise.
      template <class dataType>
        int execute() const;

      inline int setInputDataPointer(int idx, void *data){
        if(idx < numberOfInputs_){
          inputData_[idx] = data;
        }
        else{
          return -1;
        }
        return 0;
      }



      inline int setVertexNumber(const int &vertexNumber){
        vertexNumber_ = vertexNumber;
        return 0;
      }


      inline int setNumberOfInputs(int numberOfInputs){
        numberOfInputs_ = numberOfInputs;
        if(inputData_)
          free(inputData_);
        inputData_ = (void **) malloc(numberOfInputs*sizeof(void *));
        for(int i=0 ; i<numberOfInputs ; i++){
          inputData_[i] = NULL;
        }
        return 0;
      }




    protected:

      int                   vertexNumber_;
      int                   numberOfInputs_;
      void                  **inputData_; //TODO : std::vector<void*>

  };
}

// if the package is a pure template class, uncomment the following line
// #include                  <PersistenceDiagramsBarycenter.cpp>

// template functions
template <class dataType> int ttk::PersistenceDiagramsBarycenter::execute() const{

	Timer t;
	
		
	for(int i=0; i<numberOfInputs_; i++){
		std::cout << i <<std::endl;
	}
	/*for(int i=0; i<numberOfInputs_; i++){
		dataType* x = static_cast<dataType*>(inputData_[i]);
		std::cout << x[0] <<std::endl;
	}*/
	
	{
	std::stringstream msg;
	msg << "[PersistenceDiagramsBarycenter] Data-set (" << vertexNumber_
		<< " points) processed in "
		<< t.getElapsedTime() << " s. (" << threadNumber_
		<< " thread(s))."
		<< std::endl;
	dMsg(std::cout, msg.str(), timeMsg);
	}

	return 0;
}

#endif // PERSISTENCEDIAGRAMSBARYCENTER_H
