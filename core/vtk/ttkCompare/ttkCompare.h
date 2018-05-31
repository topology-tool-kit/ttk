/// \ingroup vtk
/// \class ttkCompare
/// \author Charles Gueunet
/// \date 2018-05-31
///
/// \brief TTK VTK-filter that wraps the compare processing package.
///
/// VTK wrapping code for the @Compare package.
/// This packages compare two data sets and highlight
/// their difference for both mesh and scalar values.
///
/// \param Input two scalar fields (vtkDataSet)
/// \param Output the first scalar fild with difference highlighted (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \sa ttk::Compare
#pragma once

// VTK includes -- to adapt
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>

// ttk code includes
#include <Compare.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCompare
#else
class ttkCompare
#endif
    : public vtkDataSetAlgorithm,
      public ttk::Wrapper
{
  public:
   static ttkCompare *New();
   vtkTypeMacro(ttkCompare, vtkDataSetAlgorithm)

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

   vtkSetMacro(meshOnly, bool);
   vtkGetMacro(meshOnly, bool);

   // TODO-2
   // Over-ride the input types.
   // By default, this filter has one input and one output, of the same type.
   // Here, you can re-define the input types, on a per input basis.
   // In this example, the first input type is forced to vtkUnstructuredGrid.
   // The second input type is forced to vtkImageData.
   //     int FillInputPortInformation(int port, vtkInformation *info) override {
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

  protected:
   ttkCompare()
   {
      // init
      meshOnly    = false;
      UseAllCores = false;

      SetNumberOfInputPorts(2);
   }

   ~ttkCompare(){};

   TTK_SETUP();

   // Tools

   template <typename VTKArrayType>
   vtkSmartPointer<VTKArrayType> createVTKArray(const char *      name,
                                                const std::size_t size,
                                                const char        nbEl = 1);

  private:
   // if true, do not check point/cell data
   bool         meshOnly;
   ttk::Compare compare_;
};

template <typename VTKArrayType>
vtkSmartPointer<VTKArrayType> ttkCompare::createVTKArray(const char *      name,
                                                         const std::size_t size,
                                                         const char        nbEl)
{
   vtkSmartPointer<VTKArrayType> resArray = vtkSmartPointer<VTKArrayType>::New();
   resArray->SetName(name);
   resArray->SetNumberOfComponents(nbEl);
   resArray->SetNumberOfTuples(size);

   return resArray;
}
