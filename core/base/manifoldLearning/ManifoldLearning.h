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

      inline int setInputMatrixDimensions(int numberOfRows, int numberOfColumns){
        numberOfRows_=numberOfRows;
        numberOfColumns_=numberOfColumns;
        return 0;
      }

      inline int setInputMatrix(void* data){
        matrix_=data;
        return 0;
      }

      inline int setInputMethod(int method){
        method_=method;
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

       inline int setOutputComponents(std::vector<std::vector<double>>* data){
        embedding_=data;
        return 0;
       }

       int execute() const;

    protected:
      std::string modulePath_;
      std::string moduleName_;
      std::string functionName_;
      int numberOfRows_;
      int numberOfColumns_;
      int method_;
      int numberOfComponents_;
      int numberOfNeighbors_;
      void* matrix_;
      std::vector<std::vector<double>>* embedding_;
      char majorVersion_;
  };
}
