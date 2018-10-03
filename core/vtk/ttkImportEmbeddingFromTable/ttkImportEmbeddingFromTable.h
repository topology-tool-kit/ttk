/// \ingroup vtk
/// \class ttkImportEmbeddingFromTable
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date September 2018
///
/// \brief TTK VTK-filter that embeds a vtkPointSet with the data of a vtkTable
/// 
/// \param Input Input point set (vtkPointSet)
/// \param Input Input table (vtkTable)
/// \param Output Output point set (vtkPointSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
#pragma once

// VTK includes
#include<vtkCharArray.h>
#include<vtkDataArray.h>
#include<vtkPointSet.h>
#include<vtkTable.h>
#include<vtkPointSetAlgorithm.h>
#include<vtkDoubleArray.h>
#include<vtkFiltersCoreModule.h>
#include<vtkFloatArray.h>
#include<vtkInformation.h>
#include<vtkIntArray.h>
#include<vtkObjectFactory.h>
#include<vtkPointData.h>
#include<vtkSmartPointer.h>

// ttk code includes
#include<Wrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkImportEmbeddingFromTable
#else
class ttkImportEmbeddingFromTable
#endif
: public vtkPointSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkImportEmbeddingFromTable* New();
    vtkTypeMacro(ttkImportEmbeddingFromTable, vtkPointSetAlgorithm)

      // default ttk setters
      vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = ttk::OsCall::getNumberOfCores();
      }
      Modified();
    }

    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
        SetThreads();
    }

    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    vtkSetMacro(XColumn, std::string);
    vtkGetMacro(XColumn, std::string);

    vtkSetMacro(YColumn, std::string);
    vtkGetMacro(YColumn, std::string);

    vtkSetMacro(ZColumn, std::string);
    vtkGetMacro(ZColumn, std::string);

    vtkSetMacro(Embedding2D, bool);
    vtkGetMacro(Embedding2D, bool);

    int FillInputPortInformation(int port, vtkInformation *info){
      if(port == 0)
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
      if(port == 1)
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");

      return 1;
    }

  protected:

    ttkImportEmbeddingFromTable(){
      UseAllCores = true;
      ThreadNumber = 1;
      debugLevel_ = 3;

      SetNumberOfInputPorts(2);
    }

    ~ttkImportEmbeddingFromTable(){};

    int RequestData(vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector);

  private:

    std::string XColumn;
    std::string YColumn;
    std::string ZColumn;
    bool Embedding2D;

    bool UseAllCores;
    int ThreadNumber;

    int doIt(vtkPointSet* inputDataSet, vtkTable* inputTable,  vtkPointSet* output);
    bool needsToAbort();
    int updateProgress(const float &progress);
};
