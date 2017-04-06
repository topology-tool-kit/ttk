/// \ingroup vtkWrappers
/// \class vtkTriangulationRequest
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the triangulationRequest processing package.
///
/// VTK wrapping code for the @TriangulationRequest package.
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
/// \sa ttk::TriangulationRequest
#pragma once

// ttk code includes
#include                  <ttkWrapper.h>

// VTK includes
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkSmartPointer.h>

class VTKFILTERSCORE_EXPORT vtkTriangulationRequest 
: public vtkDataSetAlgorithm, public Wrapper{

  public:

    enum Simplex{
      Vertex=0,
      Edge,
      Triangle,
      Tetra
    };

    enum Request{
      ComputeSimplex=0,
      ComputeFacet,
      ComputeCofacet,
      ComputeStar,
      ComputeLink
    };

    static vtkTriangulationRequest* New();
    vtkTypeMacro(vtkTriangulationRequest, vtkDataSetAlgorithm)

      // default ttk setters
      vtkSetMacro(debugLevel_, int);

    void SetThreadNumber(int threadNumber){\
      ThreadNumber = threadNumber;\
        SetThreads();\
    }\
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    vtkSetMacro(SimplexType, int);
    vtkGetMacro(SimplexType, int);

    vtkSetMacro(SimplexId, int);
    vtkGetMacro(SimplexId, int);

    vtkSetMacro(RequestType, int);
    vtkGetMacro(RequestType, int);

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
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); 
          break;

        default:
          break;
      }

      return 1;
    }

  protected:

    vtkTriangulationRequest(){
      UseAllCores=false;
      SimplexType=0;
      SimplexId=0;
      RequestType=0;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }

    ~vtkTriangulationRequest(){};

    TTK_SETUP();

  private:

    int SimplexType;
    int SimplexId;
    int RequestType;

};
