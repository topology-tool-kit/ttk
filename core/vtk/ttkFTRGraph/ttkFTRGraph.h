/// \ingroup vtk
/// \class ttkFTRGraph
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the fTRGraph processing package.
///
/// VTK wrapping code for the @FTRGraph package.
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
/// \sa ttk::FTRGraph
#pragma once

// VTK includes -- to adapt
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

// ttk code includes
#include <ttkWrapper.h>
#include <FTRGraph.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkFTRGraph
#else
class ttkFTRGraph
#endif
    : public vtkDataSetAlgorithm,
      public Wrapper
{
  private:
   string        ScalarField;
   vtkDataArray* outputScalarField_;
   ftr::FTRGraph fTRGraph_;

  public:
   static ttkFTRGraph* New();
   vtkTypeMacro(ttkFTRGraph, vtkDataSetAlgorithm)

       // default ttk setters
       vtkSetMacro(debugLevel_, int);

   void SetThreadNumber(int threadNumber)
   {
      ThreadNumber = threadNumber;
      SetThreads();
   }
   void SetUseAllCores(bool onOff)
   {
      UseAllCores = onOff;
      SetThreads();
   }
   // end of default ttk setters

   vtkSetMacro(ScalarField, string);
   vtkGetMacro(ScalarField, string);
   // end of TODO-4

   // TODO-2
   // Over-ride the input types.
   // By default, this filter has one input and one output, of the same type.
   // Here, you can re-define the input types, on a per input basis.
   // In this example, the first input type is forced to vtkUnstructuredGrid.
   // The second input type is forced to vtkImageData.
   //     int FillInputPortInformation(int port, vtkInformation *info){
   //
   //       switch(port){
   //         case 0:
   //           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
   //           break;
   //         case 1:
   //           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
   //           break;
   //         default:
   //           break;
   //       }
   //
   //       return 1;
   //     }
   // end of TODO-2

   // TODO-3
   // Over-ride the output types.
   // By default, this filter has one input and one output, of the same type.
   // Here, you can re-define the output types, on a per output basis.
   // In this example, the first output type is forced to vtkUnstructuredGrid.
   // The second output type is forced to vtkImageData.
   //     int FillOutputPortInformation(int port, vtkInformation *info){
   //
   //       switch(port){
   //         case 0:
   //           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
   //           break;
   //         case 1:
   //           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
   //           break;
   //         default:
   //           break;
   //       }
   //
   //       return 1;
   //     }
   // end of TODO-3

  protected:
   ttkFTRGraph()
   {
      // init
      outputScalarField_  = nullptr;
      UseAllCores = false;

      // TODO-1
      // Specify the number of input and output ports.
      // By default, this filter has one input and one output.
      // In this example, we define 2 inputs and 2 outputs.
      //       SetNumberOfInputPorts(2);
      //       SetNumberOfOutputPorts(2);
      // end of TODO-1
   }

   ~ttkFTRGraph(){};

   TTK_SETUP();
};
