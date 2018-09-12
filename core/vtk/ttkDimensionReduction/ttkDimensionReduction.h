/// \ingroup vtk
/// \class ttkDimensionReduction
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the dimensionReduction processing package.
///
/// VTK wrapping code for the @DimensionReduction package.
/// 
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa ttk::DimensionReduction
#pragma once

#include<vtkCharArray.h>
#include<vtkDataArray.h>
#include<vtkDataSet.h>
#include<vtkDoubleArray.h>
#include<vtkFiltersCoreModule.h>
#include<vtkFloatArray.h>
#include<vtkInformation.h>
#include<vtkIntArray.h>
#include<vtkObjectFactory.h>
#include<vtkPointData.h>
#include<vtkSmartPointer.h>
#include<vtkTable.h>
#include<vtkTableAlgorithm.h>

#include<DimensionReduction.h>
#include<ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkDimensionReduction
#else
class ttkDimensionReduction
#endif
: public vtkTableAlgorithm, public ttk::Wrapper{
  public:

    enum Method{
      SpectralEmbedding=0,
      LocallyLinearEmbedding,
      MDS,
      TSNE,
      Isomap,
      PCA
    };

    static ttkDimensionReduction* New();
    vtkTypeMacro(ttkDimensionReduction, vtkTableAlgorithm)

      // default ttk setters
      vtkSetMacro(debugLevel_, int);

    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = ttk::OsCall::getNumberOfCores();
      }
      Modified();
    }

    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    void SetScalarFields(std::string s){
      ScalarFields.push_back(s);
      Modified();
    }

    void ClearScalarFields(){
      ScalarFields.clear();
      Modified();
    }

    vtkSetMacro(ModulePath, std::string);
    vtkGetMacro(ModulePath, std::string);

    vtkSetMacro(ModuleName, std::string);
    vtkGetMacro(ModuleName, std::string);

    vtkSetMacro(FunctionName, std::string);
    vtkGetMacro(FunctionName, std::string);

    vtkSetMacro(NumberOfComponents, int);
    vtkGetMacro(NumberOfComponents, int);

    vtkSetMacro(NumberOfNeighbors, int);
    vtkGetMacro(NumberOfNeighbors, int);

    vtkSetMacro(IsDeterministic, int);
    vtkGetMacro(IsDeterministic, int);

    vtkSetMacro(Method, int);
    vtkGetMacro(Method, int);

    vtkSetMacro(KeepAllDataArrays, int);
    vtkGetMacro(KeepAllDataArrays, int);

    int FillInputPortInformation(int port, vtkInformation *info) override {
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
          break;
      }

      return 1;
    }

    int FillOutputPortInformation(int port, vtkInformation *info) override {
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
          break;
      }

      return 1;
    }

  protected:

    ttkDimensionReduction(){
      UseAllCores = true;
      ThreadNumber = 1;
      debugLevel_ = 3;

      outputData_=new std::vector<std::vector<double>>;
    }
    ~ttkDimensionReduction(){
      delete outputData_;
    };

    int RequestData(vtkInformation *request,
        vtkInformationVector **inputVector, vtkInformationVector *outputVector);

  private:

    int doIt(vtkTable* input, vtkTable* output);
    bool needsToAbort();
    int updateProgress(const float &progress);

    std::string ModulePath;
    std::string ModuleName;
    std::string FunctionName;
    int NumberOfComponents;
    int NumberOfNeighbors;
    int Method;
    int IsDeterministic;
    bool KeepAllDataArrays;
    bool UseAllCores;
    ttk::ThreadId ThreadNumber;
    ttk::DimensionReduction dimensionReduction_;

    std::vector<std::string> ScalarFields;
    std::vector<std::vector<double>>* outputData_;

};
