/// \ingroup vtk
/// \class ttkPointDataSelector
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date December 2017
///
/// \brief TTK VTK-filter that selects scalar fields on input with shallow copy.
/// 
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
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
#include<vtkDataSet.h>
#include<vtkDataSetAlgorithm.h>
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
class VTKFILTERSCORE_EXPORT ttkPointDataSelector
#else
class ttkPointDataSelector
#endif
: public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkPointDataSelector* New();
    vtkTypeMacro(ttkPointDataSelector, vtkDataSetAlgorithm)

      // default ttk setters
      vtkSetMacro(debugLevel_, int);
      vtkSetMacro(RegexpString, std::string);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = ttk::OsCall::getNumberOfCores();
      }
      Modified();
    }

    void SetThreadNumber(int threadNumber){\
      ThreadNumber = threadNumber;\
        SetThreads();\
    }\
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

    int FillInputPortInformation(int port, vtkInformation *info){
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
          break;
        default:
          break;
      }

      return 1;
    }

    int FillOutputPortInformation(int port, vtkInformation *info){

      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
          break;
        default:
          break;
      }

      return 1;
    }

  protected:

    ttkPointDataSelector(){
      UseAllCores = false;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }

    ~ttkPointDataSelector(){};

    int RequestData(vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector);

  private:

    bool UseAllCores;
    int ThreadNumber;
    std::vector<std::string> ScalarFields;
    std::string RegexpString;

    int doIt(vtkDataSet *input, vtkDataSet *output);
    bool needsToAbort();
    int updateProgress(const float &progress);
};
