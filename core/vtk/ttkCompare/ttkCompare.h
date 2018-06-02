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
  private:
   // if true, do not check point/cell data
   bool                meshOnly;
   int                 diffBFlag_;
   ttk::Compare        compare_;

  public:
   static ttkCompare* New();
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

   // return the bit flag corresponding to the diff
   // 0 : same,
   // 1 : difference.
   //
   // first bit: vertices
   // second bit: cell
   // third bit: vertices scalars
   // fourth bit: cells scalars
   int getDiffBFlag(void) const
   {
      return diffBFlag_;
   }

  protected:
   ttkCompare()
   {
      // init
      meshOnly    = false;
      UseAllCores = false;
      diffBFlag_  = 0;

      SetNumberOfInputPorts(2);
   }

   ~ttkCompare(){};

   TTK_SETUP();

  private:
   void computeMeshDiff(vtkDataSet* input1, vtkDataSet* input2, vtkDataSet* output);

   void computeVertsDiff(vtkDataSet* input1, vtkDataSet* input2, vtkDataSet* output);

   void computeCellsDiff(vtkDataSet* input1, vtkDataSet* input2, vtkDataSet* output);

   void addFlagFieldData(vtkDataSet* output);

   // Tools

   template <typename VTKArrayType>
   vtkSmartPointer<VTKArrayType> createVTKArray(const char*       name,
                                                const std::size_t size,
                                                const char        nbEl = 1);
};

template <typename VTKArrayType>
vtkSmartPointer<VTKArrayType> ttkCompare::createVTKArray(const char*       name,
                                                         const std::size_t size,
                                                         const char        nbEl)
{
   vtkSmartPointer<VTKArrayType> resArray = vtkSmartPointer<VTKArrayType>::New();
   resArray->SetName(name);
   resArray->SetNumberOfComponents(nbEl);
   resArray->SetNumberOfTuples(size);

   return resArray;
}
