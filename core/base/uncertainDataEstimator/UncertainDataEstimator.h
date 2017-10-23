/// \ingroup baseCode
/// \class ttk::UncertainDataEstimator
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \date August 2016.
///
/// \brief TTK processing package that takes an input ensemble data set 
/// (represented by a list of scalar fields) and which computes various 
/// vertexwise statistics (PDF estimation, bounds, moments, etc.)
///
/// \sa vtkUncertainDataEstimator.cpp %for a usage example.

#ifndef _UNCERTAINDATAESTIMATOR_H
#define _UNCERTAINDATAESTIMATOR_H

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
        #pragma omp parallel for num_threads(threadNumber_)
        for (size_t i = 0 ; i < numberOfVertices ; i++) {
          upperBound_[i] = inputData[i];
          lowerBound_[i] = inputData[i];
        }
      } else { /* Update the two fields with the new input */
        #pragma omp parallel for num_threads(threadNumber_)
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

    pair<dataType,dataType> getRange() const {
      pair<dataType,dataType> range;
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
    vector<dataType> upperBound_;
    vector<dataType> lowerBound_;
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

    void getVertexHistogram(const unsigned int vertexId, vector<double> &histogram) const {
      histogram.resize(numberOfBins_);
      if(vertexId < numberOfVertices_) {
        #ifdef _WIN32
        #pragma omp parallel for num_threads(threadNumber_)
        #else
        #pragma omp parallel for num_threads(threadNumber_) \
          schedule(static, numberOfBins_/threadNumber_)
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
      #ifdef _WIN32
      #pragma omp parallel for num_threads(threadNumber_)
      #else
      #pragma omp parallel for num_threads(threadNumber_) collapse(2) \
        schedule(static, (numberOfBins_*numberOfVertices_)/threadNumber_)
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
    vector<double>            binValue_;
    vector<vector<double> >   probability_;
    unsigned int              numberOfBins_;
    unsigned int              numberOfInputs_;
    unsigned int              numberOfVertices_;
    double                    rangeMin_;
    double                    rangeMax_;
    // vector<int> selection_; // TODO : selection support
  };

  class UncertainDataEstimator : public Debug{

    public:

      UncertainDataEstimator();

      ~UncertainDataEstimator();

      /// Execute the package.
      /// \return Returns 0 upon success, negative values otherwise.
      template <class dataType>
        int execute() const;

      /// Pass a pointer to an input array representing a scalarfield.
      /// The array is expected to be correctly allocated. idx in [0,numberOfInputs_[
      /// \param idx Index of the input scalar field.
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa setNumberOfInputs() and setVertexNumber().
      inline int setInputDataPointer(int idx, void *data){
        if(idx < numberOfInputs_){
          inputData_[idx] = data;
        }
        else{
          return -1;
        }
        return 0;
      }

      /// Pass a pointer to an output array representing the lower bound scalar field.
      /// The array is expected to be correctly allocated.
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa setVertexNumber()
      inline int setOutputLowerBoundField(void *data){
        outputLowerBoundField_ = data;
        return 0;
      }

      /// Pass a pointer to an output array representing the upper bound scalar field.
      /// The array is expected to be correctly allocated.
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa setVertexNumber()
      inline int setOutputUpperBoundField(void *data){
        outputUpperBoundField_ = data;
        return 0;
      }

      inline int setOutputProbability(int idx, void *data){
        if(idx < binCount_){
          outputProbability_[idx] = data;
        }

        return 0;
      }

      inline int setOutputMeanField(void *data) {
        outputMeanField_ = data;
        return 0;
      }


      inline int setComputeLowerBound(const bool &state){
        computeLowerBound_ = state;
        return 0;
      }


      inline int setComputeUpperBound(const bool &state){
        computeUpperBound_ = state;
        return 0;
      }


      /// Set the number of vertices in the scalar field.
      /// \param vertexNumber Number of vertices in the data-set.
      /// \return Returns 0 upon success, negative values otherwise.
      inline int setVertexNumber(const int &vertexNumber){
        vertexNumber_ = vertexNumber;
        return 0;
      }

      inline int setBinCount(const int &binCount){
        binCount_ = binCount;
        if(outputProbability_)
          free(outputProbability_);
        outputProbability_ = (void **) malloc(binCount*sizeof(void *));
        for(int b=0 ; b<binCount ; b++)
          outputProbability_[b] = NULL;

        if(binValues_)
          free(binValues_);
        binValues_ = (double *) malloc(binCount*sizeof(double));

        return 0;
      }


      /// Set the number of input scalar fields
      /// \param numberOfInputs Number of input scalar fields.
      /// \return Returns 0 upon success, negative values otherwise
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


      inline double getBinValue(int b){
        if(b<binCount_)
          return binValues_[b];
        return 0.0;
      }


    protected:

      int                   vertexNumber_;
      int                   numberOfInputs_;
      int                   binCount_;
      double                *binValues_; //TODO : vector<double>
      bool                  computeLowerBound_;
      bool                  computeUpperBound_;
      void                  **inputData_; //TODO : vector<void*>
      void                  *outputLowerBoundField_;
      void                  *outputUpperBoundField_;
      void                  **outputProbability_; //TODO : vector<void*>
      void                  *outputMeanField_;
  };
}

// if the package is a pure template class, uncomment the following line
// #include                  <UncertainDataEstimator.cpp>

// template functions
template <class dataType> int UncertainDataEstimator::execute() const{

  Timer t;

  // Check the consistency of the variables
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfInputs_)
    return -1;
  if(!vertexNumber_)
    return -2;
  if(!inputData_)
    return -3;

  for(int i=0 ; i<numberOfInputs_ ; i++){
      if(!inputData_[i])
        return -4;
    }
  if(!outputLowerBoundField_)
    return -5;
  if(!outputUpperBoundField_)
    return -6;
#endif

  int count = 0;

  // Pointers type casting
  dataType *outputLowerBoundField = (dataType *) outputLowerBoundField_;
  dataType *outputUpperBoundField = (dataType *) outputUpperBoundField_;
  double **outputProbability = (double **) outputProbability_;
  dataType **inputData = (dataType **) inputData_;
  double *outputMeanField = static_cast<double*>(outputMeanField_);


#ifdef TTK_ENABLE_OPENMP
  omp_lock_t writeLock;
  omp_init_lock(&writeLock);
#pragma omp parallel for num_threads(threadNumber_)
#endif

  for(int v = 0; v < (int) vertexNumber_; v++){

    // Avoid any processing if the abort signal is sent
    if((!wrapper_)||((wrapper_)&&(!wrapper_->needsToAbort()))){

      // For the lower bound scalar field
      if(computeLowerBound_){
        // Initialisation : values of the first input
        outputLowerBoundField[v] = inputData[0][v];
        // Loop over the inputs
        for(int inp=1 ; inp < numberOfInputs_ ; inp++){
          // Minimum value
          if(computeLowerBound_)
            if(inputData[inp][v] < outputLowerBoundField[v])
              outputLowerBoundField[v] = inputData[inp][v];
        }
      }

      // For the upper bound scalar field
      if(computeUpperBound_){
        // Initialisation : values of the first input
        outputUpperBoundField[v] = inputData[0][v];
        // Loop over the inputs
        for(int inp=1 ; inp < numberOfInputs_ ; inp++){
          // Maximum value
          if(computeUpperBound_)
            if(inputData[inp][v] > outputUpperBoundField[v])
              outputUpperBoundField[v] = inputData[inp][v];
        }
      }

      // Update the progress bar of the wrapping code -- to adapt
      if(debugLevel_ > advancedInfoMsg){
#ifdef TTK_ENABLE_OPENMP
        omp_set_lock(&writeLock);
#endif
        if((wrapper_)
          &&(!(count % ((vertexNumber_)/10)))){
          wrapper_->updateProgress((count + 1.0)
            /vertexNumber_);
        }

        count++;
#ifdef TTK_ENABLE_OPENMP
        omp_unset_lock(&writeLock);
#endif
      }
    }
  }


  // Histogram
  if(computeUpperBound_ && computeLowerBound_){
    // Range
    double range[2];
    range[0] = outputLowerBoundField[0];
    range[1] = outputUpperBoundField[0];

    for(int v=0 ; v<vertexNumber_ ; v++){
      if(outputLowerBoundField[v] < range[0])
        range[0] = outputLowerBoundField[v];
      if(outputUpperBoundField[v] > range[1])
        range[1] = outputUpperBoundField[v];
    }

    // Interval between bins
    double dx = (range[1]-range[0]) / (double)binCount_;

    // Bin values
    for(int b=0 ; b<binCount_ ; b++){
      binValues_[b] = range[0]+(dx/2.0) + (double)b * dx;
    }

    int idx;
    double increment = 1.0 / (double)numberOfInputs_;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for private(idx) num_threads(threadNumber_)
#endif
    for(int v=0 ; v<vertexNumber_ ; v++){
      for(int i=0 ; i<numberOfInputs_ ; i++){
        idx = (unsigned int) floor((inputData[i][v]-range[0])*binCount_/(range[1]-range[0]));
        idx = (idx==binCount_) ? binCount_-1 : idx;
        outputProbability[idx][v] += increment;
      }
    }
  }

  // Mean field
  for(int v=0 ; v<vertexNumber_ ; v++) {
    double sum = 0.0;
    for(int i=0 ; i<numberOfInputs_ ; i++) {
      sum += static_cast<double>(inputData[i][v]);
    }
    outputMeanField[v] = sum / static_cast<double>(numberOfInputs_);
  }



#ifdef TTK_ENABLE_OPENMP
  omp_destroy_lock(&writeLock);
#endif

  {
    stringstream msg;
    msg << "[UncertainDataEstimator] Data-set (" << vertexNumber_
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // UNCERTAINDATAESTIMATOR_H
