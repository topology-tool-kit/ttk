/// \ingroup base
/// \class ttk::ManifoldLearning 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %manifoldLearning processing package.
///
/// %ManifoldLearning is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkManifoldLearning.cpp %for a usage example.

#pragma once

#include<Triangulation.h>
#include<Wrapper.h>

namespace ttk{

  class ManifoldLearning : public Debug{

    public:

      ManifoldLearning();
      ~ManifoldLearning();

      inline int setInputModulePath(const std::string& modulePath){
        modulePath_=modulePath;
        return 0;
      }

      inline int setInputModuleName(const std::string& moduleName){
        moduleName_=moduleName;
        return 0;
      }

      inline int setInputFunctionName(const std::string& functionName){
        functionName_=functionName;
        return 0;
      }

      inline int setInputMatrixDimension(const int matrixDimension){
        matrixDimension_=matrixDimension;
        return 0;
      }

      inline int setInputMatrix(void* data){
        matrix_=data;
        return 0;
      }

      inline int setInputNumberOfComponents(int numberOfComponents){
        numberOfComponents_=numberOfComponents;
        return 0;
      }

      inline int setInputNumberOfNeighbors(int numberOfNeighbors){
        numberOfNeighbors_=numberOfNeighbors;
        return 0;
      }

       inline int setOutputComponents(std::vector<void*>* data){
        components_=data;
        return 0;
       }

       int execute() const;

    protected:
      std::string modulePath_;
      std::string moduleName_;
      std::string functionName_;
      int matrixDimension_;
      int numberOfComponents_;
      int numberOfNeighbors_;
      void* matrix_;
      std::vector<void*>* components_;
      char majorVersion_;
  };
}
