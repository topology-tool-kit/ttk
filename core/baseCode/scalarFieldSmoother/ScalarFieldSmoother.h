/// \ingroup baseCode
/// \class ttk::ScalarFieldSmoother
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2014.
///
/// \brief TTK processing package for scalar field smoothing.
///
/// This class is a dummy example for the development of TTK classes. It 
/// smooths an input scalar field by averaging the scalar values on the link
/// of each vertex.
///
/// \param dataType Data type of the input scalar field (char, float, 
/// etc.).
///
/// \sa vtkScalarFieldSmoother.cpp %for a usage example.

#ifndef _SCALAR_FIELD_SMOOTHER_H
#define _SCALAR_FIELD_SMOOTHER_H

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>

namespace ttk{

  class ScalarFieldSmoother : public Debug{

    public:
        
      ScalarFieldSmoother();
      
      ~ScalarFieldSmoother();

      int setDimensionNumber(const int &dimensionNumber){
        dimensionNumber_ = dimensionNumber;
        return 0;
      }
      
      int setInputDataPointer(void *data){
        inputData_ = data;
        return 0;
      }
      
      int setOutputDataPointer(void *data){
        outputData_ = data;
        return 0;
      }

      int setMaskDataPointer(void *mask){
        mask_ = (char*)mask;
        return 0;
      }
    
      inline int setupTriangulation(Triangulation *triangulation){
        
        triangulation_ = triangulation;
        
        // Pre-condition functions.
        if(triangulation_){
          triangulation_->preprocessVertexNeighbors();
        }
        
        return 0;
      }

      template <class dataType>
        int smooth(const int &numberOfIterations) const;
      
    protected:
    
      int                   dimensionNumber_;
      void                  *inputData_, *outputData_;
      char                  *mask_;
      Triangulation         *triangulation_;
  };

}

// template functions
template <class dataType> int ScalarFieldSmoother::smooth(
  const int &numberOfIterations) const{

  Timer t;
  
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
  if(!dimensionNumber_)
    return -2;
  if(!inputData_)
    return -3;
  if(!outputData_)
    return -4;
#endif

  int count = 0;
  
  int vertexNumber = triangulation_->getNumberOfVertices();
  
  vector<dataType> tmpData(vertexNumber*dimensionNumber_, 0);
  
  dataType *outputData = (dataType *) outputData_;
  dataType *inputData = (dataType *) inputData_;
  
  // init the output
  for(int i = 0; i < vertexNumber; i++){
    for(int j = 0; j < dimensionNumber_; j++){
      outputData[dimensionNumber_*i + j] = inputData[dimensionNumber_*i + j];
    }
  }

  for(int it = 0; it < numberOfIterations; it++){
#ifdef TTK_ENABLE_OPENMP
    omp_lock_t writeLock;
    omp_init_lock(&writeLock);
#pragma omp parallel for num_threads(threadNumber_) 
#endif
    for(int i = 0; i < vertexNumber; i++){
    
      // avoid to process masked vertices
      if(mask_ != nullptr && mask_[i] == 0) continue;

      // avoid any processing if the abort signal is sent
      if((!wrapper_)||((wrapper_)&&(!wrapper_->needsToAbort()))){
        
        for(int j = 0; j < dimensionNumber_; j++){
          tmpData[dimensionNumber_*i + j] = 0;
          
          int neighborNumber = triangulation_->getVertexNeighborNumber(i);
          for(int k = 0; k < (int) neighborNumber; k++){
            int neighborId = -1;
            triangulation_->getVertexNeighbor(i, k, neighborId);
            tmpData[dimensionNumber_*i + j]
              += outputData[dimensionNumber_*(neighborId) + j];
          }
          tmpData[dimensionNumber_*i + j]
            /= ((double) neighborNumber);
        }

        if(debugLevel_ > advancedInfoMsg){
#ifdef TTK_ENABLE_OPENMP
          omp_set_lock(&writeLock);
#endif
          // update the progress bar of the wrapping code
          if((wrapper_)
            &&(!(count % ((numberOfIterations*vertexNumber)/10)))){
            wrapper_->updateProgress((count + 1.0)
              /(numberOfIterations*vertexNumber));
          }
          count++;
#ifdef TTK_ENABLE_OPENMP
          omp_unset_lock(&writeLock);
#endif
        }
      }
    }
   
    if(numberOfIterations){
      // assign the tmpData back to the output
      for(int i = 0; i < vertexNumber; i++){
        for(int j = 0; j < dimensionNumber_; j++){
           // only set value for unmasked points
           if (mask_ == nullptr || mask_[i] != 0)
           {
              outputData[dimensionNumber_ * i + j] = tmpData[dimensionNumber_ * i + j];
           }
        }
      }
    }
    
#ifdef TTK_ENABLE_OPENMP
    omp_destroy_lock(&writeLock);
#endif
  }
  
  {
    stringstream msg;
    msg << "[ScalarFieldSmoother] Data-set (" << vertexNumber
      << " points) smoothed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return 0;
}

#endif // _SCALAR_FIELD_SMOOTHER_H
