/// \ingroup vtk
/// \class ttkDepthImageBasedGeometryApproximation
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.07.2018
///
/// \brief TTK VTK-filter that wraps the depthImageBasedGeometryApproximation processing package.
///
/// VTK wrapping code for the @DepthImageBasedGeometryApproximation package.
///
/// \param Input Depth Image (vtkImageData)
/// \param Output Unstructured Grid (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::DepthImageBasedGeometryApproximation
#pragma once

// VTK includes -- to adapt
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
#include                  <vtkCellData.h>
#include                  <vtkSmartPointer.h>

// ttk code includes
#include                  <DepthImageBasedGeometryApproximation.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkDepthImageBasedGeometryApproximation
#else
class ttkDepthImageBasedGeometryApproximation
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkDepthImageBasedGeometryApproximation* New();
    vtkTypeMacro(ttkDepthImageBasedGeometryApproximation, vtkDataSetAlgorithm)

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

    vtkSetMacro(Downsampling, int);
    vtkGetMacro(Downsampling, int);

    // Over-ride the input types.
    int FillInputPortInformation(int port, vtkInformation *info) override {

      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
          break;
        default:
          break;
      }

      return 1;
    }

    // Over-ride the output types.
    int FillOutputPortInformation(int port, vtkInformation *info) override {

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

    ttkDepthImageBasedGeometryApproximation(){

      Downsampling = 0;

      UseAllCores = false;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }

    ~ttkDepthImageBasedGeometryApproximation(){};

    TTK_SETUP();


  private:

    int                   Downsampling;
    ttk::DepthImageBasedGeometryApproximation            depthImageBasedGeometryApproximation_;

};
