/// \ingroup vtk
/// \class ttkManifoldLearning
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the manifoldLearning processing package.
///
/// VTK wrapping code for the @ManifoldLearning package.
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
/// \sa ttk::ManifoldLearning
#pragma once

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

#include<ManifoldLearning.h>
#include<ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkManifoldLearning
#else
class ttkManifoldLearning
#endif
: public vtkDataSetAlgorithm, public ttk::Wrapper{
  public:

    static ttkManifoldLearning* New();
    vtkTypeMacro(ttkManifoldLearning, vtkDataSetAlgorithm)

      // default ttk setters
      vtkSetMacro(debugLevel_, int);

    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    int FillInputPortInformation(int port, vtkInformation *info) override {
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
          break;
      }

      return 1;
    }

    int FillOutputPortInformation(int port, vtkInformation *info) override {
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
          break;
      }

      return 1;
    }

  protected:

    ttkManifoldLearning(){
      UseAllCores = true;
      ThreadNumber = 1;
      debugLevel_ = 3;
    }
    ~ttkManifoldLearning(){};

    TTK_SETUP();

  private:

    ttk::ManifoldLearning manifoldLearning_;

};
